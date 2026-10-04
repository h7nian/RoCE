#!/usr/bin/env python3
"""Pair completed high-dimensional method and baseline runs by data and seed."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import summarize_repeat_pilot as review

def identity(row):
    return (row["config"], int(row["K"]), int(row["p"]), float(row["rho"]),
            int(row["sim_id"]), row["deviation_mechanism"])


def manifest(root):
    with (root / "manifest.csv").open() as stream:
        return list(csv.DictReader(stream))


def completed_result(root, task):
    directory = root / "tasks" / task["task_id"]
    if not (directory / "COMPLETE").is_file() or (directory / "status.txt").read_text().strip() != "COMPLETE":
        raise ValueError("This report requires committed results: " + str(directory))
    with (directory / "results.csv").open() as stream:
        records = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
    methods = {row["method"]: row for row in records}
    if len(methods) != len(records):
        raise ValueError("Multiple recipes or duplicated method rows: " + str(directory))
    for row in records:
        if not all(math.isfinite(float(row[key])) for key in ("estimate", "se", "truth")) or float(row["se"]) <= 0:
            raise ValueError("Invalid committed result: " + str(directory))
    return methods, (directory / "data_sha256.txt").read_text().strip()


def summarize_paired_estimates(records):
    groups = {}
    for row in records:
        key = (row["config"], int(row["K"]), int(row["p"]), float(row["rho"]),
               row["n_deviated_sites"], row["method"])
        groups.setdefault(key, []).append(row)
    metrics = []
    for (scenario, sources, dimension, rho, deviated, method), rows in sorted(groups.items()):
        truths = [row["truth"] for row in rows]
        if max(truths) - min(truths) > 1e-12:
            raise ValueError("Population truth changes within a cell")
        errors = [row["estimate"] - row["truth"] for row in rows]
        count = len(rows)
        sd = statistics.stdev(errors)
        standard_errors = [row["se"] for row in rows]
        covered = sum(abs(error) <= 1.959963984540054 * se for error, se in zip(errors, standard_errors))
        lower, upper = review.wilson_interval(covered, count)
        mse_differences = [error**2 - row["target_error"]**2 for error, row in zip(errors, rows)]
        target_mse = statistics.mean(row["target_error"]**2 for row in rows)
        mse = statistics.mean(error**2 for error in errors)
        metrics.append(dict(config=scenario, K=sources, p=dimension, rho=rho, n_deviated_sites=deviated,
            method=method, repeats=count,
            bias=statistics.mean(errors), bias_mcse=sd/math.sqrt(count), rmse=math.sqrt(mse),
            empirical_sd=sd, mean_se=statistics.mean(standard_errors),
            se_to_sd_ratio=statistics.mean(standard_errors)/sd,
            coverage=covered/count, coverage_lower=lower, coverage_upper=upper,
            mean_ci_width=2*1.959963984540054*statistics.mean(standard_errors),
            rmse_ratio_to_target=math.sqrt(mse/target_mse),
            paired_mse_difference=statistics.mean(mse_differences),
            paired_mse_difference_mcse=statistics.stdev(mse_differences)/math.sqrt(count)))
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method_root", type=Path)
    parser.add_argument("baseline_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    method_root, baseline_root, output = map(pilot.scratch_path,
        (args.method_root, args.baseline_root, args.output))
    method_configuration = pilot.verify(method_root)
    baseline_configuration = pilot.verify(baseline_root)
    for field in ("n_per_site", "dgp_type", "dgp_control", "n_folds"):
        if method_configuration[field] != baseline_configuration[field]:
            raise ValueError("Method/baseline design differs: " + field)
    tasks = [row for row in manifest(method_root) if int(row["p"]) in (100, 200)]
    deviation_count = int(method_configuration.get("n_deviated_sites", 1))
    if (any(float(row["rho"]) != 0 for row in tasks) and
            deviation_count != int(baseline_configuration.get("n_deviated_sites", 1))):
        raise ValueError("Method/baseline deviation counts differ")
    baseline_rows = manifest(baseline_root)
    baseline_tasks = {identity(row): row for row in baseline_rows}
    if len(baseline_tasks) != len(baseline_rows):
        raise ValueError("Duplicated baseline seeds")
    if len({identity(row) for row in tasks}) != len(tasks):
        raise ValueError("Duplicated method seeds")

    def pair(task):
        key = identity(task)
        if key not in baseline_tasks:
            raise ValueError("Missing baseline seed: " + str(key))
        method, method_hash = completed_result(method_root, task)
        baseline, baseline_hash = completed_result(baseline_root, baseline_tasks[key])
        if method_hash != baseline_hash:
            raise ValueError("Different generated data: " + str(key))
        for field in ("estimate", "se", "truth"):
            if abs(float(method["target_only_ate"][field]) - float(baseline["target_only_ate"][field])) > 1e-12:
                raise ValueError("Ordinary target reference mismatch: " + str((key, field)))
        combined = dict(method)
        combined.update(baseline)
        target_error = float(combined["target_only_ate"]["estimate"]) - float(combined["target_only_ate"]["truth"])
        result = []
        for name, row in combined.items():
            if abs(float(row["truth"]) - float(combined["target_only_ate"]["truth"])) > 1e-12:
                raise ValueError("Population truth mismatch: " + str(key))
            result.append(dict(task, method=name, estimate=float(row["estimate"]), se=float(row["se"]),
                               n_deviated_sites=deviation_count, truth=float(row["truth"]),
                               target_error=target_error, data_sha256=method_hash))
        return result

    # Only independent file reads run concurrently; map retains manifest order.
    with ThreadPoolExecutor(max_workers=8) as pool:
        records = [record for paired in pool.map(pair, tasks) for record in paired]
    metrics = summarize_paired_estimates(records)
    output.mkdir(parents=True, exist_ok=True)
    review.write_csv(output / "paired_estimates.csv", records)
    review.write_csv(output / "metrics.csv", metrics)
    pilot.replace_json(output / "pairing_checks.json", dict(repeats=len(tasks),
        data_hash_matches=len(tasks), target_reference_matches=len(tasks),
        method_root=str(method_root), baseline_root=str(baseline_root)))
    lines = ["# Completed high-dimensional method/baseline comparison", "",
             "All {} repeats passed generated-data hash and ordinary-target estimate/SE/truth pairing checks.".format(len(tasks)),
             "All methods within each cell use identical seeds. This report does not pool correlated scenarios.", "",
             "| K | p | Scenario | rho | Deviated count | Method | N | Bias | RMSE | Coverage | SE/SD |",
             "|---:|---:|---|---:|---:|---|---:|---:|---:|---:|---:|"]
    for row in metrics:
        lines.append("| {K} | {p} | {config} | {rho} | {n_deviated_sites} | {method} | {repeats} | {bias:.5f} | {rmse:.5f} | {coverage:.1%} | {se_to_sd_ratio:.3f} |".format(**row))
    lines += ["", "Coverage uncertainty and paired MSE differences are in metrics.csv.",
              "SS/IVW average site effects; transportability of those averages requires its own conditions.",
              "The frozen federated/pooled baselines use their documented density models and inference conventions."]
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print("Verified {} matched high-dimensional repeats; saved {}".format(len(tasks), output))


if __name__ == "__main__":
    main()
