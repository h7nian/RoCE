#!/usr/bin/env python3
"""Review a complete, paired Two-layer calibration ablation at 200 seeds/cell."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import math
from pathlib import Path
import shutil
import statistics
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import data_identity, read_campaign
from review_completed_highdim import completed_result
from review_layer_comparison import metrics, paired_variance_comparison
from summarize_repeat_pilot import write_csv

MODES = {"common_tate": "one_round_crossfit_ate",
         "separate_arms": "one_round_crossfit_ate_separate_arms",
         "joint_tate": "one_round_crossfit_ate_joint_tate"}
QUANTITIES = ("mu1", "mu0", "tate")
CONTROLS = ("n_per_site", "dgp_type", "dgp_control", "n_folds", "nlambda",
            "nuisance_tol", "source_radius", "target_score_radius",
            "aggregation_lambda", "recipes", "nuisance_cv_certificate", "library")


def validate_configurations(full, ablation, component):
    for configuration in (full, ablation):
        if (configuration.get("crossfit_layers") != 2 or
                configuration.get("source_validation_method") != "outer_fit"):
            raise ValueError("This review requires Two-layer outer_fit scores")
    for field in CONTROLS:
        if full.get(field) != ablation.get(field):
            raise ValueError("Paired control differs: " + field)
    expected = (("calibrated", "hou_calibrated"),
                ("standard", "hou_calibrated") if component == "source" else
                ("calibrated", "lasso"))
    for configuration, pair in zip((full, ablation), expected):
        actual = (configuration.get("source_nuisance_method", "calibrated"),
                  configuration.get("target_nuisance_method", "hou_calibrated"))
        if actual != pair:
            raise ValueError("Unexpected source/target calibration combination")


def index_tasks(tasks, configuration, dimension):
    selected = {}
    seeds = defaultdict(set)
    for task in tasks:
        if (int(task["p"]) != dimension or int(task["K"]) != 2 or
                float(task["rho"]) != 0 or task["config"] not in ("C1", "C2", "C3")):
            continue
        if task["protocol"] != "one_round":
            raise ValueError("The current calibration panel requires one_round")
        key = data_identity(task, configuration)
        if key in selected:
            raise ValueError("Duplicate calibration data setting")
        selected[key] = task
        seeds[task["config"]].add(int(task["sim_id"]))
    if set(seeds) != {"C1", "C2", "C3"} or any(
            value != set(range(1, 201)) for value in seeds.values()):
        raise ValueError("C1-C3 must each contain exactly the prescribed seeds1-200")
    return selected


def require_close(first, second, label):
    if not math.isclose(float(first), float(second), rel_tol=0, abs_tol=1e-12):
        raise ValueError("Paired invariant differs: " + label)


def validate_references(full, ablation, full_hash, ablation_hash, component):
    if full_hash != ablation_hash:
        raise ValueError("Paired generated-data hashes differ")
    for result in (full, ablation):
        if not (set(MODES.values()) | {"target_only_ate"}).issubset(result):
            raise ValueError("A requested estimator is missing")
    for field in ("estimate", "se", "truth"):
        require_close(full["target_only_ate"][field], ablation["target_only_ate"][field],
                      "ordinary-target " + field)
        if component == "source":
            require_close(full["target_anchor_ate"][field], ablation["target_anchor_ate"][field],
                          "calibrated-target " + field)
    if component == "target":
        for method in (MODES["separate_arms"], MODES["joint_tate"]):
            for arm in ("mu1", "mu0"):
                for source in ("s1", "s2"):
                    field = arm + "_source_" + source + "_source_estimate"
                    require_close(full[method][field], ablation[method][field], field)


def read_layers(directory, results):
    with (directory / "layer_estimates.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    indexed = {(row["aggregation_mode"], row["quantity"]): row for row in rows}
    expected = {(mode, quantity) for mode in MODES for quantity in QUANTITIES}
    if len(rows) != len(expected) or set(indexed) != expected:
        raise ValueError("Missing or duplicated arm/TATE summaries")
    for row in rows:
        if row["recipe"] != "score_derivative" or int(row["crossfit_layers"]) != 2:
            raise ValueError("Unexpected layer-summary configuration")
        for field in ("estimate", "truth", "variance", "variance_fixed_weights",
                      "covariance_mu1_mu0", "covariance_mu1_mu0_fixed_weights"):
            if not math.isfinite(float(row[field])):
                raise ValueError("Nonfinite arm/TATE summary")
        if min(float(row["variance"]), float(row["variance_fixed_weights"])) <= 0:
            raise ValueError("Nonpositive arm/TATE variance")
    for mode, method in MODES.items():
        treated, control, tate = [indexed[mode, quantity] for quantity in QUANTITIES]
        for field in ("estimate", "truth"):
            require_close(float(treated[field])-float(control[field]), tate[field],
                          "TATE " + field + " decomposition")
            require_close(tate[field], results[method][field], "result/layer " + field)
        require_close(tate["variance"], float(results[method]["se"])**2, "result/layer variance")
        for variance, covariance in (("variance", "covariance_mu1_mu0"),
                                     ("variance_fixed_weights", "covariance_mu1_mu0_fixed_weights")):
            require_close(treated[covariance], control[covariance], "arm covariance")
            require_close(float(treated[variance])+float(control[variance])-2*float(treated[covariance]),
                          tate[variance], "TATE variance decomposition")
    return indexed


def paired_summary(full, ablation):
    if len(full) != len(ablation) or len(full) < 3:
        raise ValueError("At least three paired repeats are required")
    for first, second in zip(full, ablation):
        require_close(first["truth"], second["truth"], "paired truth")
    first_errors = [float(row["estimate"])-float(row["truth"]) for row in full]
    second_errors = [float(row["estimate"])-float(row["truth"]) for row in ablation]
    differences = [x*x-y*y for x, y in zip(first_errors, second_errors)]
    coverage_difference = [int(abs(x) <= 1.959963984540054*math.sqrt(float(a["variance"]))) -
                           int(abs(y) <= 1.959963984540054*math.sqrt(float(b["variance"])))
                           for x, y, a, b in zip(first_errors, second_errors, full, ablation)]
    count = len(full)
    mcse = statistics.stdev(differences)/math.sqrt(count)
    variance = paired_variance_comparison(first_errors, second_errors,
        [float(row["variance"]) for row in full], [float(row["variance"]) for row in ablation])
    return dict(repeats=count, mse_difference=statistics.mean(differences),
        mse_difference_mcse=mcse,
        mse_difference_z=statistics.mean(differences)/mcse if mcse else None,
        relative_rmse_change=math.sqrt(statistics.mean(x*x for x in first_errors)/
                                      statistics.mean(y*y for y in second_errors))-1,
        coverage_difference=statistics.mean(coverage_difference),
        coverage_difference_mcse=statistics.stdev(coverage_difference)/math.sqrt(count),
        empirical_variance_difference=variance["empirical_variance_difference_two_minus_three"],
        reported_variance_difference=variance["mean_reported_variance_difference_two_minus_three"],
        variance_difference_discrepancy=variance["variance_difference_discrepancy"],
        variance_discrepancy_jackknife_se=variance["variance_discrepancy_jackknife_se"])


def review(full_root, ablation_root, output, component, dimension):
    full_root, ablation_root, output = map(pilot.scratch_path, (full_root, ablation_root, output))
    configurations = [pilot.verify(root) for root in (full_root, ablation_root)]
    validate_configurations(*configurations, component)
    tasks = [index_tasks(read_campaign(root)[1], configuration, dimension)
             for root, configuration in zip((full_root, ablation_root), configurations)]
    if set(tasks[0]) != set(tasks[1]):
        raise ValueError("Full and ablated data specifications differ")

    def read_pair(key):
        entries = []
        for root, indexed in zip((full_root, ablation_root), tasks):
            task = indexed[key]
            results, data_hash = completed_result(root, task)
            entries.append((results, data_hash, read_layers(root / "tasks" / task["task_id"], results)))
        validate_references(entries[0][0], entries[1][0], entries[0][1], entries[1][1], component)
        for layer_key in entries[0][2]:
            require_close(entries[0][2][layer_key]["truth"], entries[1][2][layer_key]["truth"], "arm truth")
        return tasks[0][key], entries

    with ThreadPoolExecutor(max_workers=8) as pool:
        pairs = list(pool.map(read_pair, sorted(tasks[0])))
    groups = defaultdict(list)
    rows, target_records = [], []
    for task, entries in pairs:
        for variant, entry in zip(("full", "ablation"), entries):
            for (mode, quantity), row in entry[2].items():
                groups[task["config"], mode, quantity, variant].append(row)
                rows.append(dict(config=task["config"], p=dimension, K=2, rho=0,
                    sim_id=int(task["sim_id"]), variant=variant, data_sha256=entry[1], **row))
        for method in ("target_only_ate", "target_anchor_ate"):
            row = entries[0][0][method]
            target_records.append(dict(config=task["config"], p=dimension, K=2, rho=0,
                sim_id=int(task["sim_id"]), method=method, estimate=row["estimate"],
                truth=row["truth"], variance=float(row["se"])**2,
                variance_fixed_weights=float(row["se"])**2))
    summaries = [dict(config=config, p=dimension, K=2, rho=0, aggregation_mode=mode,
        quantity=quantity, variant=variant, **metrics(values, 200))
        for (config, mode, quantity, variant), values in sorted(groups.items())]
    comparisons = []
    for config in ("C1", "C2", "C3"):
        for mode in MODES:
            for quantity in QUANTITIES:
                comparisons.append(dict(config=config, p=dimension, K=2, rho=0,
                    aggregation_mode=mode, quantity=quantity,
                    **paired_summary(groups[config, mode, quantity, "full"],
                                     groups[config, mode, quantity, "ablation"])))
    target_groups = defaultdict(list)
    for row in target_records:
        target_groups[row["config"], row["method"]].append(row)
    target_metrics = [dict(config=config, p=dimension, K=2, rho=0, method=method,
        **metrics(values, 200)) for (config, method), values in sorted(target_groups.items())]
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "paired_estimates.csv", rows)
    write_csv(output / "metrics.csv", summaries)
    write_csv(output / "paired_comparisons.csv", comparisons)
    write_csv(output / "target_references.csv", target_metrics)
    pilot.write_json(output / "pairing_checks.json", dict(pairs=len(pairs), repeats_per_cell=200,
        component=component, dimension=dimension, crossfit_layers=2,
        data_hash_matches=len(pairs), ordinary_target_matches=len(pairs),
        calibrated_target_matches=len(pairs) if component == "source" else 0,
        source_candidate_matches=len(pairs) if component == "target" else 0,
        configurations=configurations, roots=list(map(str, (full_root, ablation_root))),
        hashes={str(root / name): pilot.digest(root / name)
                for root in (full_root, ablation_root) for name in ("configuration.json", "manifest.csv", "provenance.json")}))
    workflow = output / "workflow"
    workflow.mkdir()
    for name in ("review_current_calibration.py", "review_layer_comparison.py",
                 "review_completed_highdim.py", "prepare_highdim_validation.py",
                 "submit_repeat_pilot.py", "summarize_repeat_pilot.py", "advance_validation_plan.py"):
        shutil.copy2(Path(__file__).with_name(name), workflow / name)
    pilot.write_json(output / "workflow_hashes.json", {p.name: pilot.digest(p) for p in workflow.iterdir()})
    lines = ["# Complete Two-layer calibration ablation", "",
        "Component: {}. p={}. All600 pairs use identical data and200 prescribed seeds per C1-C3 cell.".format(component, dimension),
        "The full method retains both calibrations. Source ablation retains initial balancing and omits final source calibration; target ablation uses an ordinary lasso target anchor with calibrated sources.",
        "Data, fold configuration, cutoff2, nuisance grid100, fitting library and aggregation program are matched. Each variant relearns its own eta; identical eta values are not imposed.",
        "Every contrast is full minus ablation; negative MSE/RMSE changes favor full calibration. Paired MCSEs use independent repeats, not folds or correlated scenarios.", "",
        "| Scenario | Quantity | Full RMSE | Ablated RMSE | RMSE change | Paired MSE / MCSE | Full coverage | Ablated coverage |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    lookup = {(row["config"], row["aggregation_mode"], row["quantity"], row["variant"]): row for row in summaries}
    for row in comparisons:
        if row["aggregation_mode"] != "joint_tate":
            continue
        key = row["config"], "joint_tate", row["quantity"]
        full, ablated = lookup[key+("full",)], lookup[key+("ablation",)]
        lines.append("|{}|{}|{:.5f}|{:.5f}|{:+.2%}|{:.2f}|{:.1%}|{:.1%}|".format(
            row["config"], row["quantity"], full["rmse"], ablated["rmse"], row["relative_rmse_change"],
            row["mse_difference_z"] or 0, full["coverage"], ablated["coverage"]))
    lines += ["", "metrics.csv includes both arms, TATE, reported/empirical variance, Wilson coverage intervals and fixed-weight variance controls for all three aggregation modes.",
        "paired_comparisons.csv includes paired MSE and coverage differences and a paired jackknife SE for the difference in reported-minus-empirical variance discrepancies.",
        "The nuisance fits and active sets are held fixed by the reported eta-sensitivity variance; these are not full nuisance-refit variances. This study does not establish uniform improvement or growing-K validity.", ""]
    (output / "report.md").write_text("\n".join(lines))
    print("Verified {} complete {}-calibration pairs; {} comparisons.".format(len(pairs), component, len(comparisons)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("full_root", type=Path)
    parser.add_argument("ablation_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--component", choices=("source", "target"), required=True)
    parser.add_argument("--dimension", type=int, choices=(100, 200), required=True)
    args = parser.parse_args()
    review(args.full_root, args.ablation_root, args.output, args.component, args.dimension)
