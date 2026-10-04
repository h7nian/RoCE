#!/usr/bin/env python3
"""Summarize paired diagnostic seeds, explicitly recording all missing completions."""
import argparse
import csv
import json
import math
import statistics
from pathlib import Path

from summarize_production import wilson_interval


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("expected_seeds", type=int)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    expected = set(range(1, args.expected_seeds + 1))
    summary = []
    status = {}
    for configuration in ("C1", "C2", "C3"):
        root = args.input_root / configuration
        if not root.is_dir():
            continue
        completed = {}
        baseline_identity_checks = 0
        for directory in sorted(root.glob("seed_*")):
            if not (directory / "COMPLETED").exists():
                continue
            with (directory / "results.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            seed = int(rows[0]["seed"])
            lambda_sensitivity = "rule" in rows[0]
            expected_variants = 2 if lambda_sensitivity else 6
            if seed not in expected or seed in completed or len(rows) != expected_variants:
                raise ValueError(f"Unexpected or duplicate seed/variant count in {directory}")
            if seed != int(directory.name[5:]) or any(
                    row["configuration"] != configuration or int(row["seed"]) != seed
                    for row in rows):
                raise ValueError(f"Filename, configuration and row identifiers disagree in {directory}")
            if any(int(row["n_site"]) != 1000 or int(row["p"]) != 100 or
                   int(row["n_folds"]) != 10 for row in rows):
                raise ValueError(f"Wrong scientific setting in {directory}")
            if lambda_sensitivity:
                if {row["rule"] for row in rows} != {"min", "1se"}:
                    raise ValueError(f"Missing or duplicate penalty rule in {directory}")
                reference_path = args.input_root.parent / "target_native" / configuration / directory.name / "results.csv"
                if reference_path.exists():
                    with reference_path.open() as handle:
                        reference = next(row for row in csv.DictReader(handle) if row["variant"] == "lasso")
                    minimum = next(row for row in rows if row["rule"] == "min")
                    if any(abs(float(minimum[field]) - float(reference["target_" + field])) > 1e-12
                           for field in ("estimate", "bias", "se")):
                        raise ValueError(f"Baseline differs between independent diagnostic runs: {directory}")
                    baseline_identity_checks += 1
                for row in rows:
                    row.update(variant="lasso" if row["rule"] == "min" else "lasso_1se",
                               target_bias=row["bias"], target_se=row["se"],
                               full_bias="NA", target_remainder="nan")
            elif {row["variant"] for row in rows} != {
                    "lasso", "oracle_outcome", "oracle_propensity", "oracle_both",
                    "rcal_propensity", "rcal_weighted"}:
                raise ValueError(f"Missing or duplicate nuisance variant in {directory}")
            completed[seed] = {row["variant"]: row for row in rows}
        status[configuration] = dict(expected=len(expected), completed=len(completed),
                                     missing=sorted(expected - set(completed)))
        if baseline_identity_checks:
            status[configuration]["independent_baseline_identity_checks"] = baseline_identity_checks
        if len(completed) < 2:
            continue
        seeds = sorted(completed)
        for scope in ("target", "full"):
            if completed[seeds[0]]["lasso"][f"{scope}_bias"] == "NA":
                continue
            for variant in completed[seeds[0]]:
                biases = [float(completed[seed][variant][f"{scope}_bias"]) for seed in seeds]
                standard_errors = [float(completed[seed][variant][f"{scope}_se"]) for seed in seeds]
                baseline = [float(completed[seed]["lasso"][f"{scope}_bias"]) for seed in seeds]
                if not all(math.isfinite(value) for value in biases + standard_errors + baseline):
                    raise ValueError("Nonfinite successful-replicate output")
                deltas = [value - reference for value, reference in zip(biases, baseline)]
                mse_deltas = [value**2 - reference**2 for value, reference in zip(biases, baseline)]
                remainders = [float(completed[seed][variant]["target_remainder"]) for seed in seeds]
                count = len(seeds)
                covered = sum(abs(value) <= 1.959963984540054 * se
                              for value, se in zip(biases, standard_errors))
                lower, upper = wilson_interval(covered, count)
                empirical_sd = statistics.stdev(biases)
                remainder_mean = statistics.mean(remainders) if all(map(math.isfinite, remainders)) else math.nan
                remainder_mcse = statistics.stdev(remainders) / math.sqrt(count) if all(map(math.isfinite, remainders)) else math.nan
                summary.append(dict(configuration=configuration, scope=scope, variant=variant,
                                    n=count, expected_n=len(expected), bias=statistics.mean(biases),
                                    bias_mcse=empirical_sd / math.sqrt(count), sd=empirical_sd,
                                    mean_se=statistics.mean(standard_errors),
                                    se_over_sd=statistics.mean(standard_errors) / empirical_sd,
                                    rmse=math.sqrt(statistics.mean(value**2 for value in biases)),
                                    coverage=covered / count, coverage_lower=lower, coverage_upper=upper,
                                    paired_shift=statistics.mean(deltas),
                                    paired_shift_mcse=statistics.stdev(deltas) / math.sqrt(count),
                                    paired_mse_difference=statistics.mean(mse_deltas),
                                    paired_mse_mcse=statistics.stdev(mse_deltas) / math.sqrt(count),
                                    target_remainder=remainder_mean,
                                    target_remainder_mcse=remainder_mcse))
    (args.output_root / "completion.json").write_text(json.dumps(status, indent=2))
    if summary:
        with (args.output_root / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=summary[0].keys())
            writer.writeheader()
            writer.writerows(summary)
    print(json.dumps(status))
    for row in summary:
        print("{configuration} {scope} {variant:18s} n={n} bias={bias:+.5f} RMSE={rmse:.5f} "
              "cov={coverage:.3f} remainder={target_remainder:+.5f} paired_shift={paired_shift:+.5f} "
              "+/-{paired_shift_mcse:.5f}".format(**row))


if __name__ == "__main__":
    main()
