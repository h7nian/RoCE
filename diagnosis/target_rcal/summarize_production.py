#!/usr/bin/env python3
"""Independently summarize existing v4 raw TATE rows, retaining cell counts.

This validates row consistency and recorded package provenance. It does not
replace the production checkpoint/checksum gates or call an incomplete cell final.
"""
import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def wilson_interval(successes, count):
    proportion = successes / count
    z = 1.959963984540054
    denominator = 1 + z * z / count
    center = (proportion + z * z / (2 * count)) / denominator
    radius = z * math.sqrt(proportion * (1 - proportion) / count + z * z / (4 * count * count)) / denominator
    return center - radius, center + radius


def write_paired_method_comparisons(groups, output_path):
    comparisons = []
    for key, direct in sorted(groups.items()):
        family, configuration, source_count, rho, method = key
        if method != "one_round_crossfit_ate":
            continue
        armwise = groups[(family, configuration, source_count, rho,
                          "one_round_crossfit_ate_armwise")]
        seeds = sorted(set(direct) & set(armwise))
        for metric, value in (("estimate", lambda row: row[0]),
                              ("squared_error", lambda row: row[0] ** 2),
                              ("coverage", lambda row: int(row[2]))):
            differences = [value(direct[seed]) - value(armwise[seed]) for seed in seeds]
            comparisons.append(dict(family=family, configuration=configuration, K=source_count,
                                    rho=rho, metric=metric, n=len(seeds),
                                    direct_minus_armwise=statistics.mean(differences),
                                    mcse=statistics.stdev(differences) / math.sqrt(len(seeds))))
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparisons[0].keys())
        writer.writeheader()
        writer.writerows(comparisons)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--compare-only", action="store_true",
                        help="Read the existing thin replicate table and refresh paired comparisons only.")
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    groups = defaultdict(dict)
    if args.compare_only:
        with (args.output_root / "replicate_metrics.csv").open() as handle:
            for row in csv.DictReader(handle):
                key = (row["family"], row["configuration"], int(row["K"]),
                       float(row["rho"]), row["method"])
                groups[key][int(row["seed"])] = (float(row["bias"]), float(row["se"]),
                                                 row["covered"] == "True", float(row["truth"]))
        write_paired_method_comparisons(groups, args.output_root / "paired_direct_armwise.csv")
        return
    fingerprints = set()
    errors = []
    file_count = 0
    for family, root in (("negative_transfer", args.production_root),
                         ("shared_shift", args.production_root / "shared_shift")):
        for path in sorted((root / "raw").glob("task_*.csv")):
            file_count += 1
            try:
                with path.open(newline="") as handle:
                    rows = list(csv.DictReader(handle))
                for row in rows:
                    if row["estimand_scope"] != "tate":
                        continue
                    n_total = int(row["n_total"])
                    source_count = int(row["K"])
                    if n_total != 1000 * (source_count + 1) or int(row["p"]) != 100 or int(row["n_folds"]) != 10:
                        raise ValueError("unexpected scientific settings")
                    bias, se, estimate = (float(row[key]) for key in ("bias", "se", "estimate"))
                    if not all(math.isfinite(value) for value in (bias, se, estimate)) or se <= 0:
                        raise ValueError("nonfinite estimate/bias or invalid SE")
                    covered = abs(bias) <= 1.959963984540054 * se
                    if covered != (row["coverage"] == "TRUE"):
                        raise ValueError("coverage disagrees with bias and SE")
                    fingerprints.add(row["package_fingerprint"])
                    key = (family, row["config"], source_count, float(row["rho"]), row["method"])
                    seed = int(row["sim_id"])
                    if seed in groups[key]:
                        raise ValueError("duplicate replicate in cell")
                    groups[key][seed] = (bias, se, covered, estimate - bias)
            except (KeyError, ValueError, OSError) as error:
                errors.append({"file": str(path), "error": str(error)})
    if errors:
        (args.output_root / "errors.json").write_text(json.dumps(errors, indent=2))
        raise RuntimeError(f"{len(errors)} raw-file consistency failures; inspect errors.json")
    if len(fingerprints) != 1:
        raise RuntimeError(f"Expected a single frozen package fingerprint, found {fingerprints}")
    summary = []
    for key, replicates in sorted(groups.items()):
        family, configuration, source_count, rho, method = key
        records = list(replicates.values())
        count = len(records)
        biases = [record[0] for record in records]
        truth_values = [record[3] for record in records]
        if max(truth_values) - min(truth_values) > 1e-10:
            raise ValueError(f"inconsistent truth in {key}")
        bias = statistics.mean(biases)
        sd = statistics.stdev(biases)
        mean_se = statistics.mean(record[1] for record in records)
        coverage_count = sum(record[2] for record in records)
        centered_coverage = sum(abs(record[0] - bias) <= 1.959963984540054 * record[1]
                                for record in records) / count
        lower, upper = wilson_interval(coverage_count, count)
        target = groups[(family, configuration, source_count, rho, "target_only_ate")]
        paired_seeds = sorted(set(target) & set(replicates))
        mse_differences = [replicates[seed][0] ** 2 - target[seed][0] ** 2 for seed in paired_seeds]
        paired_target_mse = statistics.mean(target[seed][0] ** 2 for seed in paired_seeds)
        paired_method_mse = statistics.mean(replicates[seed][0] ** 2 for seed in paired_seeds)
        summary.append(dict(family=family, configuration=configuration, K=source_count, rho=rho,
                            method=method, n=count, bias=bias, bias_mcse=sd / math.sqrt(count),
                            sd=sd, mean_se=mean_se, se_over_sd=mean_se / sd,
                            rmse=math.sqrt(statistics.mean(value * value for value in biases)),
                            coverage=coverage_count / count, coverage_lower=lower, coverage_upper=upper,
                            coverage_after_mc_recentering=centered_coverage,
                            paired_n=len(paired_seeds),
                            rmse_ratio_target=math.sqrt(paired_method_mse / paired_target_mse),
                            paired_mse_difference=statistics.mean(mse_differences),
                            paired_mse_mcse=statistics.stdev(mse_differences) / math.sqrt(len(paired_seeds))))
    with (args.output_root / "cell_metrics.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    with (args.output_root / "replicate_metrics.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("family", "configuration", "K", "rho", "method", "seed",
                         "bias", "se", "covered", "truth"))
        for key, replicates in sorted(groups.items()):
            for seed, record in sorted(replicates.items()):
                writer.writerow((*key, seed, *record))
    (args.output_root / "audit.json").write_text(json.dumps(dict(
        raw_files=file_count, cells=len(summary), package_fingerprint=next(iter(fingerprints)),
        invalid_files=0, checkpoint_gates_repeated=False), indent=2))
    write_paired_method_comparisons(groups, args.output_root / "paired_direct_armwise.csv")
    primary = [row for row in summary if row["method"] == "one_round_crossfit_ate"]
    print(f"Checked {file_count} files and {len(summary)} TATE cells; {len(primary)} primary cells.")
    for row in primary:
        print("{family:18s} {configuration} K={K} rho={rho:g} n={n} coverage={coverage:.3f} "
              "bias={bias:+.5f} RMSE={rmse:.5f} RMSE/target={rmse_ratio_target:.3f} SE/SD={se_over_sd:.3f}".format(**row))


if __name__ == "__main__":
    main()
