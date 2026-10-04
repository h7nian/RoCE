#!/usr/bin/env python3
"""Report completed ablations and compare identical data/seed pairs only."""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import summarize_repeat_pilot as review


def records(root):
    with (root / "manifest.csv").open() as stream:
        manifest = list(csv.DictReader(stream))
    result = {}
    for task in manifest:
        directory = root / "tasks" / task["task_id"]
        if not (directory / "COMPLETE").is_file():
            continue
        if (directory / "status.txt").read_text().strip() != "COMPLETE":
            raise ValueError("Inconsistent completion: " + str(directory))
        with (directory / "results.csv").open() as stream:
            rows = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
        methods = {row["method"]: row for row in rows}
        if len(methods) != len(rows):
            raise ValueError("Expected one recipe and unique methods: " + str(directory))
        joint = methods[task["protocol"] + "_crossfit_ate_joint_tate"]
        target = methods["target_only_ate"]
        key = (task["config"], task["p"], task["K"], float(task["rho"]), task["sim_id"])
        result[key] = dict(task=task, joint=joint, target=target,
                           data_hash=(directory / "data_sha256.txt").read_text().strip())
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    extension = args.run_root / "extended_validation_v1"
    roots = dict(main=args.run_root / "validation_highdim_mc200_cutoff2_v1")
    roots.update({name: extension / name for name in (
        "outcome_shift", "covariate_shift_half", "covariate_shift_double",
        "target_lasso", "initial_validation", "two_round")})
    summaries = {}
    for name, root in roots.items():
        review.summarize(root, args.output / name)
        with (args.output / name / "metrics.csv").open() as stream:
            summaries[name] = [row for row in csv.DictReader(stream)
                               if row["method"].endswith("_joint_tate")]
    reference = records(roots["main"])
    reference.update(records(roots["outcome_shift"]))
    comparisons = []
    for name in ("target_lasso", "initial_validation", "two_round"):
        groups = defaultdict(list)
        for key, candidate in records(roots[name]).items():
            baseline = reference.get(key)
            if baseline is None:
                continue
            if candidate["data_hash"] != baseline["data_hash"]:
                raise ValueError("Paired data hashes differ: " + str((name, key)))
            for field in ("estimate", "se", "truth"):
                if not math.isclose(float(candidate["target"][field]),
                                    float(baseline["target"][field]), rel_tol=0, abs_tol=1e-12):
                    raise ValueError("Paired target references differ: " + str((name, key)))
            if abs(float(candidate["joint"]["truth"]) - float(baseline["joint"]["truth"])) > 1e-12:
                raise ValueError("Paired population truths differ")
            groups[key[:4]].append((candidate["joint"], baseline["joint"]))
        for cell, pairs in sorted(groups.items()):
            errors = [[float(row["estimate"]) - float(row["truth"]) for row in pair] for pair in pairs]
            differences = [a*a-b*b for a, b in errors]
            record = dict(ablation=name, config=cell[0], p=cell[1], K=cell[2], rho=cell[3],
                          n_pairs=len(pairs), mse_difference=statistics.mean(differences),
                          mse_difference_mcse=statistics.stdev(differences)/math.sqrt(len(pairs))
                          if len(pairs) > 1 else None)
            for index, label in enumerate(("ablation", "reference")):
                record[label + "_rmse"] = math.sqrt(statistics.mean(e[index]**2 for e in errors))
                record[label + "_coverage"] = statistics.mean(
                    abs(e[index]) <= 1.959963984540054 * float(pair[index]["se"])
                    for e, pair in zip(errors, pairs))
            comparisons.append(record)
    review.write_csv(args.output / "paired_ablations.csv", comparisons)
    (args.output / "joint_metrics.json").write_text(json.dumps(summaries, indent=2) + "\n")
    complete = all(int(row["n_success"]) == int(row["n_planned"])
                   for rows in summaries.values() for row in rows)
    title = "# Completed ablation review" if complete else "# Interim ablation review"
    lines = [title, "", "Only committed repeats are included. Planned cells have 50 repeats for extensions.",
             "All comparisons use joint_tate, n=1000 per site, p=100, cutoff=2.",
             "Use the matched seeds in paired_ablations.csv when comparing variants.", "",
             "| Variant | Scenario | rho | Completed | Bias | RMSE | Coverage |", "|---|---|---:|---:|---:|---:|---:|"]
    for name, rows in summaries.items():
        for row in rows:
            if row["p"] != "100":
                continue
            lines.append("| {} | {} | {} | {} | {:.5f} | {:.5f} | {:.1%} |".format(
                name, row["config"], row["rho"], row["n_success"],
                float(row["bias"]), float(row["rmse"]), float(row["coverage"])))
    lines += ["", "Paired comparisons, data-hash and ordinary-target-reference checks: paired_ablations.csv.",
              "A negative paired MSE difference favors the ablation; report its Monte Carlo uncertainty."]
    (args.output / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
