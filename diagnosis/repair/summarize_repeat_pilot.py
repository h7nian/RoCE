#!/usr/bin/env python3
"""Summarize every pilot task, retaining failures and paired recipe differences.

This is a development-pilot report, not a Monte Carlo coverage summary.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
from collections import Counter, defaultdict

CELL_FIELDS = ("config", "K", "p", "rho", "protocol", "deviation_mechanism", "n_deviated_sites",
               "source_nuisance_method", "target_nuisance_method", "source_validation_method",
               "crossfit_layers")


def wilson_interval(successes, count):
    if count == 0:
        return (None, None)
    z = 1.959963984540054
    proportion = successes / count
    denominator = 1 + z * z / count
    center = (proportion + z * z / (2 * count)) / denominator
    radius = z * math.sqrt(proportion * (1 - proportion) / count + z * z / (4 * count * count)) / denominator
    return center - radius, center + radius


def summarize_metrics(manifest, recipes, estimates, output, requested_methods=None):
    planned = Counter(tuple(row.get(key, "") for key in CELL_FIELDS) for row in manifest)
    groups = defaultdict(list)
    lookup = {}
    for row in estimates:
        if row["state"] != "COMPLETE":
            continue
        cell = tuple(row.get(key, "") for key in CELL_FIELDS)
        groups[(cell, row["recipe"], row["method"])].append(row)
        lookup[(row["task_id"], row["recipe"], row["method"])] = row
    metrics, paired_mse = [], []
    for cell, expected_count in sorted(planned.items()):
        protocol = cell[CELL_FIELDS.index("protocol")]
        prefix = protocol + "_crossfit_ate"
        methods = requested_methods or [prefix, prefix + "_separate_arms", prefix + "_joint_tate",
                                        prefix + "_armwise", "target_only_ate", "target_anchor_ate"]
        for recipe in recipes:
            for method in methods:
                rows = groups[(cell, recipe, method)]
                record = dict(zip(CELL_FIELDS, cell))
                record.update(recipe=recipe, method=method, n_planned=expected_count,
                              n_success=len(rows), n_unavailable=expected_count-len(rows))
                if not rows:
                    record.update(bias=None, bias_mcse=None, rmse=None, empirical_sd=None,
                        mean_se=None, se_to_sd_ratio=None, coverage=None, coverage_mcse=None,
                        coverage_lower=None, coverage_upper=None, mean_ci_width=None)
                else:
                    truths = [float(row["truth"]) for row in rows]
                    if max(truths) - min(truths) > 1e-12:
                        raise ValueError("Population truth changes within a simulation cell")
                    errors = [float(row["estimate"]) - float(row["truth"]) for row in rows]
                    n = len(errors)
                    sd = statistics.stdev(errors) if n > 1 else None
                    standard_errors = [float(row["se"]) for row in rows]
                    covered = sum(abs(error) <= 1.959963984540054 * se
                                  for error, se in zip(errors, standard_errors))
                    lower, upper = wilson_interval(covered, n)
                    mean_se = statistics.mean(standard_errors)
                    record.update(bias=statistics.mean(errors), bias_mcse=sd/math.sqrt(n) if sd is not None else None,
                        rmse=math.sqrt(statistics.mean(error*error for error in errors)), empirical_sd=sd,
                        mean_se=mean_se, se_to_sd_ratio=mean_se/sd if sd else None,
                        coverage=covered/n, coverage_mcse=math.sqrt((covered/n)*(1-covered/n)/n),
                        coverage_lower=lower, coverage_upper=upper,
                        mean_ci_width=2*1.959963984540054*mean_se)
                    differences = []
                    for row, error in zip(rows, errors):
                        reference = lookup.get((row["task_id"], recipe, "target_only_ate"))
                        if reference is None or not math.isclose(float(reference["truth"]), float(row["truth"]), rel_tol=0, abs_tol=1e-12):
                            raise ValueError("Unpaired target benchmark for task " + row["task_id"])
                        reference_error = float(reference["estimate"]) - float(reference["truth"])
                        differences.append(error*error-reference_error*reference_error)
                    paired_mse.append(dict(zip(CELL_FIELDS, cell), recipe=recipe, method=method,
                        n_pairs=n, mean_squared_error_difference=statistics.mean(differences),
                        mcse=statistics.stdev(differences)/math.sqrt(n) if n > 1 else None))
                metrics.append(record)
    write_csv(output / "metrics.csv", metrics)
    write_csv(output / "paired_mse_vs_target.csv", paired_mse)


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(root, output):
    root = root.resolve()
    output = output.resolve()
    scratch = Path("/scratch.global/zhan9381/FACE-HD")
    if scratch not in root.parents or scratch not in output.parents:
        raise ValueError("Pilot and report must be under the FACE-HD scratch root")
    output.mkdir(parents=True, exist_ok=True)
    with (root / "manifest.csv").open() as stream:
        manifest = list(csv.DictReader(stream))
    configuration = json.loads((root / "configuration.json").read_text())
    recipes = configuration["recipes"]
    source_method = ("not_applicable" if configuration.get("baseline_methods") else
                     configuration.get("source_nuisance_method", "calibrated"))
    target_method = ("not_applicable" if configuration.get("baseline_methods") else
                     configuration.get("target_nuisance_method", "hou_calibrated"))
    validation_method = configuration.get("source_validation_method", "calibrated")
    layers = configuration.get("crossfit_layers")
    if configuration.get("baseline_methods"):
        validation_method = layers = "not_applicable"
    elif layers == 2 or validation_method == "outer_fit":
        validation_method, layers = "outer_fit", 2
    elif layers is None:
        layers = 3 if (target_method == "hou_calibrated" or
                       (source_method == "calibrated" and validation_method in ("calibrated", "complete"))) else 2
    for task in manifest:
        task.setdefault("p", str(configuration.get("p", "")))
        task.setdefault("n_deviated_sites", str(configuration.get("n_deviated_sites", 1)))
        task.setdefault("source_nuisance_method", source_method)
        task.setdefault("target_nuisance_method", target_method)
        task.setdefault("source_validation_method", validation_method)
        task.setdefault("crossfit_layers", str(layers))
    scheduler = {}
    if (root / "scheduler_status.csv").exists():
        with (root / "scheduler_status.csv").open() as stream:
            scheduler = {row["task_id"]: row["scheduler_state"] for row in csv.DictReader(stream)}
    statuses, estimates, paired = [], [], []
    for task in manifest:
        directory = root / "tasks" / task["task_id"]
        state_file = directory / "status.txt"
        state = state_file.read_text().strip() if state_file.exists() else "NOT_STARTED"
        scheduler_state = scheduler.get(task["task_id"], "")
        if state != "COMPLETE" and scheduler_state in {
                "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL"}:
            state = scheduler_state
        if state == "COMPLETE" and not (directory / "COMPLETE").exists():
            state = "INCOMPLETE_COMMIT"
        result_path = directory / "results.csv"
        if state == "COMPLETE" and not result_path.exists():
            state = "MISSING_RESULTS"
        failures = [path.read_text().strip() for path in sorted(directory.glob("*_FAILED.txt"))]
        status = dict(task, state=state, scheduler_state=scheduler_state, errors="; ".join(failures))
        statuses.append(status)
        if not result_path.exists():
            continue
        with result_path.open() as stream:
            rows = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
        try:
            valid = bool(rows) and all(all(math.isfinite(float(row[key])) for key in ("estimate", "se", "truth"))
                                      and float(row["se"]) > 0 for row in rows)
        except (ValueError, TypeError):
            valid = False
        if not valid:
            status["state"] = state = "INVALID_RESULTS"
            status["errors"] += "; nonfinite estimate/truth, nonpositive SE, or missing TATE rows"
        methods = {}
        for row in rows:
            recipe = row["pilot_recipe"]
            method = row["method"]
            key = (recipe, method)
            if key in methods or recipe not in recipes:
                raise ValueError("Duplicate or unexpected method/recipe in task " + task["task_id"])
            methods[key] = row
            estimates.append(dict(task, state=state, recipe=recipe, method=method,
                                  estimate=row["estimate"], se=row["se"], truth=row["truth"]))
        if state != "COMPLETE" or len(recipes) < 2:
            continue
        reference = recipes[0]
        reference_methods = {method for recipe, method in methods if recipe == reference}
        for recipe in recipes:
            if {method for candidate, method in methods if candidate == recipe} != reference_methods:
                raise ValueError("Incomplete method pairing in task " + task["task_id"])
        for (recipe, method), row in sorted(methods.items()):
            if recipe == reference:
                continue
            base = methods.get((reference, method))
            if base is None or float(base["truth"]) != float(row["truth"]):
                raise ValueError("Unpaired method/truth in task " + task["task_id"])
            paired.append(dict(task, reference_recipe=reference, candidate_recipe=recipe,
                method=method, estimate_difference=float(row["estimate"]) - float(base["estimate"]),
                se_difference=float(row["se"]) - float(base["se"])))
    write_csv(output / "task_status.csv", statuses)
    write_csv(output / "tate_estimates.csv", estimates)
    write_csv(output / "paired_differences.csv", paired)
    requested_methods = ([method + "_ate" for method in configuration["baseline_methods"]] +
                         ["target_only_ate"]) if configuration.get("baseline_methods") else None
    if requested_methods is None and configuration.get("target_nuisance_method") == "lasso":
        prefix = manifest[0]["protocol"] + "_crossfit_ate"
        requested_methods = [prefix, prefix + "_separate_arms", prefix + "_joint_tate",
                             prefix + "_armwise", "target_only_ate"]
    summarize_metrics(manifest, recipes, estimates, output, requested_methods)
    complete = sum(row["state"] == "COMPLETE" for row in statuses)
    report = ["# Simulation development pilot", "",
              "Completed tasks: {}/{}.".format(complete, len(statuses)), "",
              "Every manifest task is retained in task_status.csv, including failures.",
              "tate_estimates.csv includes task status so partial results cannot be mistaken for complete repeats.",
              "paired_differences.csv uses completed tasks with matched recipes, methods and truth.", "",
              "metrics.csv reports replicate metrics, Monte Carlo uncertainty and unavailable repeats.",
              "paired_mse_vs_target.csv uses the ordinary target-only result from the identical repeat.",
              "Small or incomplete cells are preliminary; a single repeat does not establish coverage or RMSE.", "",
              "| Task | Scenario | Sources | Seed | State |", "| --- | --- | --- | --- | --- |"]
    report += ["| {task_id} | {config} | {K} | {sim_id} | {state} |".format(**row) for row in statuses]
    (output / "report.md").write_text("\n".join(report) + "\n")
    print("Completed {}/{} tasks; report: {}".format(complete, len(statuses), output / "report.md"))
    return complete == len(statuses)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if summarize(args.pilot, args.output) else 1)
