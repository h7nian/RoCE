#!/usr/bin/env python3
"""Compare final-calibration factorial cells on identical data and seeds."""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import prepare_highdim_validation as design
import submit_repeat_pilot as pilot
import summarize_repeat_pilot as summary


def read_variant(roots, label):
    records = {}
    for root in roots:
        configuration = pilot.verify(root)
        _, tasks = design.read_campaign(root)
        for task in tasks:
            if int(task["p"]) not in (100, 200) or float(task["rho"]) != 0 or int(task["sim_id"]) > 50:
                continue
            directory = root / "tasks" / task["task_id"]
            if not (directory / "COMPLETE").is_file():
                continue
            if (directory / "status.txt").read_text().strip() != "COMPLETE":
                raise ValueError("Inconsistent completion: " + str(directory))
            with (directory / "results.csv").open() as stream:
                rows = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
            methods = {row["method"]: row for row in rows}
            if len(methods) != len(rows):
                raise ValueError("Duplicate estimator in " + str(directory))
            if any(not math.isfinite(float(row[field])) for row in rows for field in ("estimate", "se", "truth")):
                raise ValueError("Nonfinite completed result in " + str(directory))
            if any(float(row["se"]) <= 0 for row in rows):
                raise ValueError("Nonpositive completed SE in " + str(directory))
            key = design.data_identity(task, configuration)
            if key in records:
                raise ValueError("Duplicate completed seed for " + label)
            records[key] = dict(task=task, methods=methods,
                                data_hash=(directory / "data_sha256.txt").read_text().strip())
    return records


def compare(left, right, label):
    groups = {}
    for key in sorted(set(left) & set(right)):
        a, b = left[key], right[key]
        if a["data_hash"] != b["data_hash"]:
            raise ValueError("Different paired data: " + str((label, key)))
        for field in ("estimate", "se", "truth"):
            if not math.isclose(float(a["methods"]["target_only_ate"][field]),
                                float(b["methods"]["target_only_ate"][field]), rel_tol=0, abs_tol=1e-12):
                raise ValueError("Different paired target reference: " + str((label, key)))
        for method in ("one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
                       "one_round_crossfit_ate_joint_tate"):
            first, second = a["methods"][method], b["methods"][method]
            if abs(float(first["truth"]) - float(second["truth"])) > 1e-12:
                raise ValueError("Different paired truths")
            errors = [float(row["estimate"]) - float(row["truth"]) for row in (first, second)]
            covered = [abs(error) <= 1.959963984540054 * float(row["se"])
                       for error, row in zip(errors, (first, second))]
            cell = (a["task"]["config"], int(a["task"]["p"]), method)
            groups.setdefault(cell, []).append((errors, covered))
    result = []
    for (scenario, dimension, method), pairs in sorted(groups.items()):
        differences = [errors[0] ** 2 - errors[1] ** 2 for errors, _ in pairs]
        result.append(dict(comparison=label, config=scenario, p=dimension, method=method,
            n_planned=50, n_pairs=len(pairs), mse_difference=statistics.mean(differences),
            mse_difference_mcse=statistics.stdev(differences)/math.sqrt(len(pairs)) if len(pairs)>1 else None,
            left_coverage=statistics.mean(covered[0] for _, covered in pairs),
            right_coverage=statistics.mean(covered[1] for _, covered in pairs),
            left_rmse=math.sqrt(statistics.mean(errors[0] ** 2 for errors, _ in pairs)),
            right_rmse=math.sqrt(statistics.mean(errors[1] ** 2 for errors, _ in pairs))))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("existing_root", type=Path)
    parser.add_argument("ablation_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    base, ablation, output = [pilot.scratch_path(path) for path in
                              (args.existing_root, args.ablation_root, args.output)]
    output.mkdir(parents=True, exist_ok=False)
    roots = dict(full=[base / "validation_highdim_mc200_cutoff2_v1"],
        target_lasso=[base / "extended_validation_v1/target_lasso", base / "highdim_validation_v2/target_lasso_p200"],
        standard_source=[ablation / ("standard_source_p" + str(p)) for p in (100, 200)],
        standard_both=[ablation / ("standard_source_target_lasso_p" + str(p)) for p in (100, 200)])
    records = {label: read_variant(paths, label) for label, paths in roots.items()}
    comparisons = []
    for left, right in (("standard_source", "full"), ("standard_both", "target_lasso"),
                        ("target_lasso", "full"), ("standard_both", "standard_source")):
        comparisons.extend(compare(records[left], records[right], left + " minus " + right))
    summary.write_csv(output / "paired_comparisons.csv", comparisons)
    status = dict(completed={label: len(rows) for label, rows in records.items()},
                  planned_per_variant=300, planned_per_cell=50,
                  roots={label: list(map(str, paths)) for label, paths in roots.items()})
    (output / "status.json").write_text(json.dumps(status, indent=2) + "\n")
    (output / "README.md").write_text(
        "# Calibration ablation review\n\n"
        "Only committed, data-hash-matched seeds are compared. Each cell has 50 prescribed repeats.\n"
        "The sign is left minus right; negative paired MSE difference favors the left variant.\n"
        "Missing pairs remain missing and are not failed-coverage events. Report Monte Carlo uncertainty.\n"
        "The standard source retains initial balancing and omits final score calibration.\n")
    print(json.dumps(status["completed"], sort_keys=True))


if __name__ == "__main__":
    main()
