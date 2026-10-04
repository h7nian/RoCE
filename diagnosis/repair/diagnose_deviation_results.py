#!/usr/bin/env python3
"""Summarize shifted-source weights alongside bias and SE calibration.

Read committed repeats only. The legacy source-one summary is retained;
all_source_weight_diagnostics.csv reports every source and its deviation flag.
"""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import summarize_repeat_pilot as review


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root, output = pilot.scratch_path(args.campaign), pilot.scratch_path(args.output)
    configuration = pilot.verify(root)
    review.summarize(root, output)
    groups = defaultdict(list)
    with (root / "manifest.csv").open() as stream:
        tasks = list(csv.DictReader(stream))
    for task in tasks:
        directory = root / "tasks" / task["task_id"]
        if not (directory / "COMPLETE").exists():
            continue
        if (directory / "status.txt").read_text().strip() != "COMPLETE":
            raise ValueError("Inconsistent completed task: " + str(directory))
        with (directory / "results.csv").open() as stream:
            rows = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
        anchor = next(row for row in rows if row["method"] == "target_anchor_ate")
        for row in rows:
            if not row["method"].endswith(("_joint_tate", "_separate_arms")):
                continue
            for source in range(1, int(task["K"]) + 1):
                for arm in (0, 1):
                    error = float(row["mu{}_source_s{}_wald_identity_error_max".format(arm, source)])
                    if not math.isfinite(error) or error > 1e-8:
                        raise ValueError("Wald discrepancy identity failed: " + str(directory))
            groups[(task["config"], int(task["K"]), task["p"], task["rho"],
                    int(configuration.get("n_deviated_sites", 1)), row["method"])].append((row, anchor))
    diagnostics = []
    all_sources = []
    fields = ["mu{}_source_s1_{}".format(arm, suffix) for arm in (0, 1)
              for suffix in ("weight", "wald_mean", "penalty_activation_fraction",
                             "screening_discrepancy_mean", "evaluation_discrepancy_mean")]
    for (config, sources, p, rho, deviated, method), pairs in sorted(groups.items()):
        errors = [float(row["estimate"]) - float(row["truth"]) for row, anchor in pairs]
        anchor_errors = [float(anchor["estimate"]) - float(anchor["truth"]) for row, anchor in pairs]
        adjustment = [error - anchor_error for error, anchor_error in zip(errors, anchor_errors)]
        standard_errors = [float(row["se"]) for row, anchor in pairs]
        count = len(pairs)
        sd = statistics.stdev(errors) if count > 1 else None
        covered = sum(abs(error) <= 1.959963984540054 * se for error, se in zip(errors, standard_errors))
        lower, upper = review.wilson_interval(covered, count)
        record = dict(config=config, K=sources, p=p, rho=rho, n_deviated_sites=deviated,
                      method=method, n=count,
                      bias=statistics.mean(errors), bias_mcse=sd/math.sqrt(count) if sd is not None else None,
                      rmse=math.sqrt(statistics.mean(error*error for error in errors)),
                      empirical_sd=sd, mean_se=statistics.mean(standard_errors),
                      coverage=covered/count, coverage_lower=lower, coverage_upper=upper,
                      anchor_bias=statistics.mean(anchor_errors),
                      mean_tate_adjustment=statistics.mean(adjustment),
                      tate_adjustment_mcse=statistics.stdev(adjustment)/math.sqrt(count) if count > 1 else None)
        for field in fields:
            record[field] = statistics.mean(float(row[field]) for row, anchor in pairs)
        diagnostics.append(record)
        for source in range(1, sources + 1):
            item = dict(config=config, K=sources, p=p, rho=rho, n_deviated_sites=deviated,
                        method=method, source="s" + str(source),
                        outcome_deviated=float(rho) != 0 and source <= deviated, n=count)
            for arm in (0, 1):
                for suffix in ("weight", "wald_mean", "penalty_activation_fraction",
                               "screening_discrepancy_mean", "evaluation_discrepancy_mean"):
                    field = "mu{}_source_s{}_{}".format(arm, source, suffix)
                    item["mu{}_{}".format(arm, suffix)] = statistics.mean(float(row[field]) for row, anchor in pairs)
            all_sources.append(item)
    review.write_csv(output / "source_weight_diagnostics.csv", diagnostics)
    review.write_csv(output / "all_source_weight_diagnostics.csv", all_sources)
    print("Saved diagnostics for {} completed cells/methods to {}".format(len(diagnostics), output))


if __name__ == "__main__":
    main()
