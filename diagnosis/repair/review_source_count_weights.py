#!/usr/bin/env python3
"""Review the first source's actual joint-TATE arm weights in completed runs."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from summarize_repeat_pilot import write_csv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root, output = map(pilot.scratch_path, (args.campaign, args.output))
    configuration = pilot.verify(root)
    if int(configuration.get("n_deviated_sites", 1)) != 1:
        raise ValueError("This review expects only source s1 to deviate")
    with (root / "manifest.csv").open() as stream:
        tasks = list(csv.DictReader(stream))

    def read_task(task):
        directory = root / "tasks" / task["task_id"]
        if not (directory / "COMPLETE").is_file():
            raise ValueError("Require all prescribed results to be committed")
        with (directory / "results.csv").open() as stream:
            records = [row for row in csv.DictReader(stream) if row["method"].endswith("_ate_joint_tate")]
        if len(records) != 1:
            raise ValueError("Expected one joint-TATE row per task")
        row = records[0]
        result = dict(task)
        for arm in ("mu1", "mu0"):
            weights = [float(row[arm+"_source_s1_fold_"+str(fold)+"_weight"])
                       for fold in range(1, int(configuration["n_folds"])+1)]
            if any(not math.isfinite(value) for value in weights):
                raise ValueError("Nonfinite weight")
            if abs(statistics.mean(weights)-float(row[arm+"_source_s1_fold_weight_mean"])) > 1e-12:
                raise ValueError("Fold means do not reproduce recorded weights")
            result[arm+"_mean_fold_weight"] = statistics.mean(weights)
            result[arm+"_fraction_zero_folds"] = statistics.mean(abs(value)<=1e-10 for value in weights)
            result[arm+"_all_folds_zero"] = all(abs(value)<=1e-10 for value in weights)
        return result

    with ThreadPoolExecutor(max_workers=8) as pool:
        records = list(pool.map(read_task, tasks))
    groups = defaultdict(list)
    for row in records:
        groups[(row["config"], row["K"], row["p"], row["rho"])].append(row)
    summary = []
    for (scenario, sources, dimension, rho), values in sorted(groups.items()):
        result = dict(config=scenario,K=sources,p=dimension,rho=rho,repeats=len(values))
        for arm in ("mu1", "mu0"):
            weights = [row[arm+"_mean_fold_weight"] for row in values]
            result[arm+"_mean_weight"] = statistics.mean(weights)
            result[arm+"_weight_mcse"] = statistics.stdev(weights)/math.sqrt(len(weights)) if len(weights)>1 else None
            result[arm+"_mean_fraction_zero_folds"] = statistics.mean(row[arm+"_fraction_zero_folds"] for row in values)
            result[arm+"_fraction_repeats_all_zero"] = statistics.mean(row[arm+"_all_folds_zero"] for row in values)
        summary.append(result)
    output.mkdir(parents=True,exist_ok=False)
    write_csv(output / "repeat_weights.csv", records)
    write_csv(output / "weight_summary.csv", summary)
    (output / "scope.json").write_text(json.dumps(dict(campaign=str(root),completed=len(records),
        method="joint_tate",source="s1",zero_tolerance=1e-10,
        note="s1 is invalid in the treated arm only when rho>0. Cross-K datasets and total sample sizes differ."),indent=2)+"\n")
    print("Verified source-s1 fold weights in {} completed repeats".format(len(records)))


if __name__ == "__main__":
    main()
