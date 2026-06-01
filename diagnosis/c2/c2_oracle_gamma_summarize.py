#!/usr/bin/env python3
"""Summarize C2 oracle-gamma convention probe outputs."""

import csv
import glob
import math
import os
from collections import defaultdict


OUTPUT_ROOT = os.environ.get(
    "C2_ORACLE_GAMMA_OUTPUT_ROOT",
    os.path.join("diagnosis", "c2", "oracle_gamma_probe"),
)
RESULT_DIR = os.path.join(OUTPUT_ROOT, "results")
SUMMARY_DIR = os.path.join(OUTPUT_ROOT, "summary")


def read_csvs(pattern):
    rows = []
    for path in sorted(glob.glob(os.path.join(RESULT_DIR, pattern))):
        with open(path, newline="") as handle:
            for row in csv.DictReader(handle):
                row["_file"] = os.path.basename(path)
                rows.append(row)
    return rows


def f(row, key):
    value = row.get(key, "")
    if value in ("", "NA", "NaN", "nan", "Inf", "-Inf"):
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def finite(values):
    return [value for value in values if math.isfinite(value)]


def mean(values):
    values = finite(values)
    return sum(values) / len(values) if values else math.nan


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def summarize_oracle_runs(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.get("method", "")].append(row)

    out = []
    for method, sub in sorted(grouped.items()):
        bias = [f(row, "bias") for row in sub]
        covered = [1.0 if row.get("covered", "") == "TRUE" else 0.0 for row in sub]
        out.append({
            "method": method,
            "n": len(sub),
            "mean_bias": mean(bias),
            "mean_abs_bias": mean([abs(value) for value in bias]),
            "min_bias": min(finite(bias)) if finite(bias) else math.nan,
            "max_bias": max(finite(bias)) if finite(bias) else math.nan,
            "coverage_rate": mean(covered),
            "mean_se": mean([f(row, "se") for row in sub]),
            "mean_weight_sum": mean([f(row, "weight_sum") for row in sub]),
        })
    return out


def compare_oracle_rows(rows):
    grouped = defaultdict(dict)
    for row in rows:
        grouped[row.get("tag", "")][row.get("method", "")] = row

    out = []
    for tag, methods in sorted(grouped.items()):
        raw = methods.get("oracle_raw_joint_gamma", {})
        cal = methods.get("oracle_source_conditional_gamma", {})
        target = methods.get("target_only", {})
        out.append({
            "tag": tag,
            "target_bias": f(target, "bias"),
            "raw_oracle_bias": f(raw, "bias"),
            "calibrated_oracle_bias": f(cal, "bias"),
            "raw_minus_calibrated_bias": f(raw, "bias") - f(cal, "bias"),
            "raw_oracle_covered": raw.get("covered", ""),
            "calibrated_oracle_covered": cal.get("covered", ""),
        })
    return out


def summarize_moments(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[row.get("gamma_label", "")].append(row)

    out = []
    for label, sub in sorted(grouped.items()):
        pop_gap = [f(row, "population_intercept") - 1.0 for row in sub]
        sample_gap = [f(row, "sample_intercept") - 1.0 for row in sub]
        out.append({
            "gamma_label": label,
            "n": len(sub),
            "mean_population_intercept_gap": mean(pop_gap),
            "mean_abs_population_intercept_gap": mean([abs(value) for value in pop_gap]),
            "mean_sample_intercept_gap": mean(sample_gap),
            "mean_abs_sample_intercept_gap": mean([abs(value) for value in sample_gap]),
            "mean_sample_weight_max": mean([f(row, "sample_weight_max") for row in sub]),
        })
    return out


def main():
    oracle_rows = read_csvs("*_oracle_runs.csv")
    moment_rows = read_csvs("*_source_gamma_moments.csv")
    if not oracle_rows:
        raise SystemExit("No oracle run CSVs found.")

    write_csv(
        os.path.join(SUMMARY_DIR, "oracle_bias_by_method.csv"),
        summarize_oracle_runs(oracle_rows),
        [
            "method", "n", "mean_bias", "mean_abs_bias", "min_bias", "max_bias",
            "coverage_rate", "mean_se", "mean_weight_sum",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "oracle_bias_by_setting.csv"),
        compare_oracle_rows(oracle_rows),
        [
            "tag", "target_bias", "raw_oracle_bias", "calibrated_oracle_bias",
            "raw_minus_calibrated_bias", "raw_oracle_covered",
            "calibrated_oracle_covered",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "source_gamma_moments.csv"),
        summarize_moments(moment_rows),
        [
            "gamma_label", "n", "mean_population_intercept_gap",
            "mean_abs_population_intercept_gap", "mean_sample_intercept_gap",
            "mean_abs_sample_intercept_gap", "mean_sample_weight_max",
        ],
    )
    print("[oracle-gamma-summary] oracle rows: %d" % len(oracle_rows))
    print("[oracle-gamma-summary] moment rows: %d" % len(moment_rows))
    print("[oracle-gamma-summary] wrote %s" % SUMMARY_DIR)


if __name__ == "__main__":
    main()
