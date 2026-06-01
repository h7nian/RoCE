#!/usr/bin/env python3
"""Summarize C2 core intermediate estimator biases from test outputs."""

import csv
import glob
import math
import os
import re
from collections import defaultdict


INPUT_DIR = os.environ.get(
    "C2_CORE_BIAS_INPUT_DIR",
    os.path.join("tests", "testthat", "c2_core_intermediate_output"),
)
OUTPUT_ROOT = os.environ.get(
    "C2_CORE_BIAS_OUTPUT_ROOT",
    os.path.join("diagnosis", "c2", "core_bias_summary"),
)
SUMMARY_DIR = os.path.join(OUTPUT_ROOT, "summary")
SINCE_PREFIX = os.environ.get("C2_CORE_BIAS_SINCE_PREFIX", "20260528")
SEEDS = {
    int(value) for value in os.environ.get(
        "C2_CORE_BIAS_SEEDS",
        "40,120,220,340,460,364,464,564",
    ).split(",") if value
}


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


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


def run_id_from_path(path, suffix):
    base = os.path.basename(path)
    if not base.endswith(suffix):
        raise ValueError(base)
    return base[: -len(suffix)]


def tuning_label(run_id):
    match = re.search(r"_seed\d+_(.+)$", run_id)
    if not match:
        return "unknown"
    label = match.group(1)
    for suffix in (
        "_method_runs", "_crossfit_sources", "_crossfit_fold_weights",
        "_crossfit_fold_components", "_crossfit_runs",
    ):
        if label.endswith(suffix):
            label = label[: -len(suffix)]
    return label


def keep_method_row(row, run_id):
    if run_id < SINCE_PREFIX:
        return False
    if row.get("config") != "C2":
        return False
    if int(f(row, "n_total")) != 5000:
        return False
    if int(f(row, "p")) != 10:
        return False
    if int(f(row, "n_folds")) != 10:
        return False
    if int(f(row, "nlambda_init")) != 100:
        return False
    return int(f(row, "seed")) in SEEDS


def latest_method_rows():
    latest = {}
    for path in glob.glob(os.path.join(INPUT_DIR, "*_method_runs.csv")):
        run_id = run_id_from_path(path, "_method_runs.csv")
        for row in read_csv(path):
            if not keep_method_row(row, run_id):
                continue
            label = tuning_label(run_id)
            key = (
                int(f(row, "seed")),
                int(f(row, "K")),
                int(f(row, "p")),
                label,
                row.get("method", ""),
            )
            row = dict(row)
            row["run_id"] = run_id
            row["tuning_label"] = label
            if key not in latest or run_id > latest[key]["run_id"]:
                latest[key] = row
    return list(latest.values())


def companion_rows(method_rows, suffix):
    run_ids = {row["run_id"] for row in method_rows}
    rows = []
    for run_id in sorted(run_ids):
        path = os.path.join(INPUT_DIR, run_id + suffix)
        if not os.path.exists(path):
            continue
        label = tuning_label(run_id)
        for row in read_csv(path):
            row = dict(row)
            row["run_id"] = run_id
            row["tuning_label"] = label
            rows.append(row)
    return rows


def method_by_setting(rows):
    out = []
    for row in sorted(rows, key=lambda r: (
        r["tuning_label"], int(f(r, "K")), int(f(r, "seed")), r.get("method", "")
    )):
        out.append({
            "tuning_label": row["tuning_label"],
            "seed": int(f(row, "seed")),
            "K": int(f(row, "K")),
            "p": int(f(row, "p")),
            "method": row.get("method", ""),
            "truth": f(row, "truth"),
            "estimate": f(row, "estimate"),
            "bias": f(row, "bias"),
            "se": f(row, "se"),
            "covered": row.get("covered", ""),
            "run_id": row["run_id"],
        })
    return out


def method_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tuning_label"], row.get("method", ""))].append(row)

    out = []
    for (label, method), sub in sorted(grouped.items()):
        bias = [f(row, "bias") for row in sub]
        covered = [1.0 if row.get("covered", "") == "TRUE" else 0.0 for row in sub]
        out.append({
            "tuning_label": label,
            "method": method,
            "n": len(sub),
            "mean_bias": mean(bias),
            "mean_abs_bias": mean([abs(value) for value in bias]),
            "min_bias": min(finite(bias)) if finite(bias) else math.nan,
            "max_bias": max(finite(bias)) if finite(bias) else math.nan,
            "coverage_rate": mean(covered),
            "mean_se": mean([f(row, "se") for row in sub]),
        })
    return out


def source_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tuning_label"], row.get("method", ""))].append(row)

    out = []
    for (label, method), sub in sorted(grouped.items()):
        bias = [f(row, "source_bias") for row in sub]
        out.append({
            "tuning_label": label,
            "method": method,
            "n": len(sub),
            "mean_source_bias": mean(bias),
            "mean_abs_source_bias": mean([abs(value) for value in bias]),
            "min_source_bias": min(finite(bias)) if finite(bias) else math.nan,
            "max_source_bias": max(finite(bias)) if finite(bias) else math.nan,
            "mean_avg_weight": mean([f(row, "avg_weight") for row in sub]),
        })
    return out


def fold_component_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tuning_label"], row.get("method", ""))].append(row)

    out = []
    for (label, method), sub in sorted(grouped.items()):
        fold_bias = [f(row, "fold_source_bias") for row in sub]
        out.append({
            "tuning_label": label,
            "method": method,
            "n": len(sub),
            "mean_fold_source_bias": mean(fold_bias),
            "mean_abs_fold_source_bias": mean([abs(value) for value in fold_bias]),
            "mean_delta_ts": mean([f(row, "delta_ts") for row in sub]),
            "mean_C_ot": mean([f(row, "C_ot") for row in sub]),
            "mean_V_s": mean([f(row, "V_s") for row in sub]),
        })
    return out


def weight_summary(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tuning_label"], row.get("method", ""))].append(row)

    out = []
    for (label, method), sub in sorted(grouped.items()):
        weights = [f(row, "weight") for row in sub]
        out.append({
            "tuning_label": label,
            "method": method,
            "n": len(sub),
            "mean_weight": mean(weights),
            "mean_abs_weight": mean([abs(value) for value in weights]),
            "min_weight": min(finite(weights)) if finite(weights) else math.nan,
            "max_weight": max(finite(weights)) if finite(weights) else math.nan,
            "mean_fold_lambda": mean([f(row, "fold_lambda") for row in sub]),
        })
    return out


def main():
    method_rows = latest_method_rows()
    if not method_rows:
        raise SystemExit("No C2 core method rows matched the requested filters.")

    source_rows = companion_rows(method_rows, "_crossfit_sources.csv")
    fold_rows = companion_rows(method_rows, "_crossfit_fold_components.csv")
    weight_rows = companion_rows(method_rows, "_crossfit_fold_weights.csv")

    write_csv(
        os.path.join(SUMMARY_DIR, "method_bias_by_setting.csv"),
        method_by_setting(method_rows),
        [
            "tuning_label", "seed", "K", "p", "method", "truth",
            "estimate", "bias", "se", "covered", "run_id",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "method_bias_by_method.csv"),
        method_summary(method_rows),
        [
            "tuning_label", "method", "n", "mean_bias", "mean_abs_bias",
            "min_bias", "max_bias", "coverage_rate", "mean_se",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "source_bias_by_method.csv"),
        source_summary(source_rows),
        [
            "tuning_label", "method", "n", "mean_source_bias",
            "mean_abs_source_bias", "min_source_bias", "max_source_bias",
            "mean_avg_weight",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "fold_component_bias_by_method.csv"),
        fold_component_summary(fold_rows),
        [
            "tuning_label", "method", "n", "mean_fold_source_bias",
            "mean_abs_fold_source_bias", "mean_delta_ts", "mean_C_ot",
            "mean_V_s",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "aggregation_weight_by_method.csv"),
        weight_summary(weight_rows),
        [
            "tuning_label", "method", "n", "mean_weight", "mean_abs_weight",
            "min_weight", "max_weight", "mean_fold_lambda",
        ],
    )

    print("[core-bias-summary] method rows: %d" % len(method_rows))
    print("[core-bias-summary] source rows: %d" % len(source_rows))
    print("[core-bias-summary] fold rows: %d" % len(fold_rows))
    print("[core-bias-summary] weight rows: %d" % len(weight_rows))
    print("[core-bias-summary] wrote %s" % SUMMARY_DIR)


if __name__ == "__main__":
    main()
