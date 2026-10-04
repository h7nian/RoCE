"""Paired component and calibration ablations; all outputs live on scratch."""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import time
import numpy as np
from experiments import binary_data
from fitted_strata import fitted_stratum_data
from calibration_remainders import calibration_remainder_budgets
from repairs import infer_repair, variance_certificate
from run_repairs import save_case


def make_plan(profile):
    counts = [32, 512, 2048] if profile == "binary" else [8, 32, 128]
    patterns = ["all_valid", "weak_boundary", "strong"] if profile == "binary" else ["all_valid", "weak_boundary"]
    rows = []
    for count in counts:
        for pattern in patterns:
            for floor in [0., .3]:
                rows.append(dict(cell_id=len(rows) + 1, panel=profile, source_count=count,
                    pattern=pattern, shared_scale=floor, n_per_site=1000,
                    source_evaluation_size=750 if profile == "binary" else 250,
                    target_evaluation_size=1000 if profile == "binary" else 500,
                    pilot_size=250, valid_fraction=.75, seed_offset=9100000))
    return rows


def run_cell(task):
    setting, repeats, output = task
    started = time.time()
    fitted = setting["panel"] == "fitted"
    data = fitted_stratum_data(setting, repeats) if fitted else binary_data(setting, repeats)
    lower, upper = variance_certificate(data["pilot_variance"], data["score_ranges"],
        250, setting["source_evaluation_size"], .005)
    common = dict(target_calibration="empirical", bounded=True,
        target_alpha=.02 if fitted else .025, variance_lower=lower, variance_upper=upper,
        variance_failure=.005, bias_failure=.005 if fitted else 0.)
    if fitted:
        common.update(source_bias_budget=data["shared_source_budget"], target_bias_budget=data["shared_target_budget"])
    results = {}
    for mode in ["dispersion", "fourier", "hybrid"]:
        results[mode + "_legacy"] = infer_repair(data, source_mode=mode, **common)
        if mode != "fourier":
            results[mode + "_direct"] = infer_repair(data, source_mode=mode,
                dispersion_calibration="bounded_exponential", **common)
    if fitted:
        valid = math.ceil(.75 * setting["source_count"])
        budgets = calibration_remainder_budgets(data, valid)
        aggregate = dict(common, source_bias_budget=budgets["source_mean_budget"],
            source_squared_bias_budget=budgets["source_squared_budget"], target_bias_budget=budgets["target_budget"])
        for mode in ["dispersion", "hybrid"]:
            results[mode + "_aggregate_budget"] = infer_repair(data, source_mode=mode,
                dispersion_calibration="bounded_exponential", **aggregate)
        # Oracle validation is recorded separately, never supplied to the fits.
        drift = data["nuisance_drift"]
        sorted_drift = np.sort(drift, axis=1)
        worst_average = np.maximum(abs(sorted_drift[:, :valid].mean(axis=1)),
                                   abs(sorted_drift[:, -valid:].mean(axis=1)))
        worst_squared = np.sort(drift**2, axis=1)[:, -valid:].sum(axis=1) / setting["source_count"]
        diagnostics = dict(
            old_source_budget=data["shared_source_budget"].mean(axis=0).tolist(),
            direct_source_budget=budgets["source_mean_budget"].mean(axis=0).tolist(),
            oracle_worst_subset_mean=worst_average.mean(axis=0).tolist(),
            old_target_budget=data["shared_target_budget"].mean(axis=0).tolist(),
            direct_target_budget=budgets["target_budget"].mean(axis=0).tolist(),
            certificate_failure=float(np.mean(np.any(worst_average > budgets["source_mean_budget"] + 1e-12, axis=1)
                | np.any(worst_squared > budgets["source_squared_budget"] + 1e-12, axis=1)
                | np.any(data["oracle_target_budget"] > budgets["target_budget"] + 1e-12, axis=1))))
        np.savez_compressed(str(Path(output) / "diagnostics" / ("budget_{:03d}.npz".format(setting["cell_id"]))),
            mean_budget=budgets["source_mean_budget"], squared_budget=budgets["source_squared_budget"],
            target_budget=budgets["target_budget"], oracle_worst_subset_mean=worst_average,
            oracle_worst_subset_squared=worst_squared)
        (Path(output) / "diagnostics" / ("cell_{:03d}.json".format(setting["cell_id"]))).write_text(json.dumps(diagnostics, indent=2) + "\n")
    return save_case(output, setting, data, results, time.time() - started)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("output")
    parser.add_argument("--profile", choices=["binary", "fitted"], required=True)
    parser.add_argument("--repeats", type=int, default=1000)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    output = Path(args.output)
    if args.repeats < 1 or not 1 <= args.workers <= 4:
        raise ValueError("Positive repeats and one to four workers required")
    if not str(output).startswith("/scratch.global/zhan9381/FACE-HD/") or output.exists():
        raise ValueError("Use a new FACE-HD scratch directory")
    plan = make_plan(args.profile)
    if args.smoke:
        plan = plan[:1]
    output.mkdir(parents=True)
    (output / "cells").mkdir()
    (output / "diagnostics").mkdir()
    hashes = {file.name: hashlib.sha256(file.read_bytes()).hexdigest() for file in Path(__file__).parent.glob("*.py")}
    (output / "configuration.json").write_text(json.dumps(dict(settings=plan, repeats=args.repeats,
        code_sha256=hashes, scope="Paired bounded-score ablation; all methods reserve .005 for variance certificates, even when unused; no high-dimensional fitted-RoCE claim"), indent=2) + "\n")
    records = []
    with multiprocessing.Pool(args.workers) as pool:
        for metrics in pool.imap_unordered(run_cell, [(row, args.repeats, str(output)) for row in plan]):
            records.extend(metrics)
            print("Completed", args.profile, "cell", metrics[0]["cell_id"], flush=True)
    with (output / "metrics.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=sorted(set(key for row in records for key in row)))
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda row: (row["cell_id"], row["method"])))
    (output / "COMPLETE").write_text("All paired repetitions retained, including failures and fallback outcomes.\n")


if __name__ == "__main__":
    main()
