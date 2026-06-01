#!/usr/bin/env python3
"""Summarize C2 nuisance-refit probe outputs.

This script is intended to run as an MSI Slurm job after the refit array.  It
uses only the Python standard library so it does not depend on local packages.
"""

import csv
import glob
import math
import os
import statistics
from collections import defaultdict


RESULT_DIR = os.path.join("diagnosis", "c2", "nuisance_refit_probe", "results")
SUMMARY_DIR = os.path.join("diagnosis", "c2", "nuisance_refit_probe", "summary")


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
    return [x for x in values if math.isfinite(x)]


def mean(values):
    vals = finite(values)
    return sum(vals) / len(vals) if vals else math.nan


def sd(values):
    vals = finite(values)
    return statistics.stdev(vals) if len(vals) >= 2 else math.nan


def median(values):
    vals = finite(values)
    return statistics.median(vals) if vals else math.nan


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def summarize_eval(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[
            (
                row.get("tag", ""),
                row.get("mode", ""),
                row.get("gamma_label", ""),
                row.get("alpha_label", ""),
            )
        ].append(row)

    out = []
    for (tag, mode, gamma_label, alpha_label), sub in sorted(grouped.items()):
        bias = [f(row, "bias_vs_fold_truth") for row in sub]
        mu_bias = [f(row, "mu_pred_bias") for row in sub]
        delta = [f(row, "delta") for row in sub]
        ess = [f(row, "treated_weight_ess") for row in sub]
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "gamma_label": gamma_label,
                "alpha_label": alpha_label,
                "n": len(sub),
                "mean_bias": mean(bias),
                "mean_abs_bias": mean([abs(x) for x in bias]),
                "median_abs_bias": median([abs(x) for x in bias]),
                "sd_bias": sd(bias),
                "mean_mu_pred_bias": mean(mu_bias),
                "mean_delta": mean(delta),
                "mean_treated_ess": mean(ess),
            }
        )
    return out


def summarize_current_vs_best(summary_rows):
    grouped = defaultdict(list)
    for row in summary_rows:
        grouped[(str(row["tag"]), str(row["mode"]))].append(row)

    out = []
    for (tag, mode), sub in sorted(grouped.items()):
        current = [
            row
            for row in sub
            if row["gamma_label"] == "current" and row["alpha_label"] == "current"
        ]
        best = min(sub, key=lambda row: abs(float(row["mean_bias"])))
        cur = current[0] if current else {}
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "current_mean_bias": cur.get("mean_bias", math.nan),
                "current_mean_abs_bias": cur.get("mean_abs_bias", math.nan),
                "best_gamma_label": best["gamma_label"],
                "best_alpha_label": best["alpha_label"],
                "best_mean_bias": best["mean_bias"],
                "best_mean_abs_bias": best["mean_abs_bias"],
                "best_mean_delta": best["mean_delta"],
                "best_mean_mu_pred_bias": best["mean_mu_pred_bias"],
            }
        )
    return out


def summarize_oracle(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row.get("tag", ""), row.get("mode", ""), row.get("variant", ""))].append(row)

    out = []
    for (tag, mode, variant), sub in sorted(grouped.items()):
        bias = [f(row, "bias_vs_fold_truth") for row in sub]
        delta = [f(row, "delta") for row in sub]
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "variant": variant,
                "n": len(sub),
                "mean_bias": mean(bias),
                "mean_abs_bias": mean([abs(x) for x in bias]),
                "median_abs_bias": median([abs(x) for x in bias]),
                "mean_delta": mean(delta),
            }
        )
    return out


def summarize_runs(rows):
    out = []
    for row in sorted(rows, key=lambda r: (r.get("tag", ""), r.get("mode", ""))):
        out.append(
            {
                "tag": row.get("tag", ""),
                "mode": row.get("mode", ""),
                "estimate": f(row, "estimate"),
                "bias_vs_superpopulation": f(row, "bias_vs_superpopulation"),
                "se": f(row, "se"),
                "covered": row.get("covered", ""),
                "target_bias_vs_superpopulation": f(row, "target_bias_vs_superpopulation"),
                "weight_sum": f(row, "weight_sum"),
            }
        )
    return out


def summarize_refit_failures(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[
            (
                row.get("tag", ""),
                row.get("mode", ""),
                row.get("nuisance", ""),
                row.get("label", ""),
                row.get("ok", ""),
            )
        ].append(row)

    out = []
    for (tag, mode, nuisance, label, ok), sub in sorted(grouped.items()):
        lambdas = [f(row, "lambda") for row in sub]
        supports = [f(row, "support") for row in sub]
        errors = sorted({row.get("error", "") for row in sub if row.get("error", "")})
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "nuisance": nuisance,
                "label": label,
                "ok": ok,
                "n": len(sub),
                "mean_lambda": mean(lambdas),
                "mean_support": mean(supports),
                "errors": " | ".join(errors[:3]),
            }
        )
    return out


def main():
    eval_rows = read_csvs("*_source_eval_grid.csv")
    oracle_rows = read_csvs("*_oracle_checks.csv")
    run_rows = read_csvs("*_run_summary.csv")
    refit_rows = read_csvs("*_nuisance_refits.csv")

    if not eval_rows:
        raise SystemExit(f"No source_eval_grid CSVs found in {RESULT_DIR}")

    eval_summary = summarize_eval(eval_rows)
    current_vs_best = summarize_current_vs_best(eval_summary)
    oracle_summary = summarize_oracle(oracle_rows)
    run_summary = summarize_runs(run_rows)
    failure_summary = summarize_refit_failures(refit_rows)

    write_csv(
        os.path.join(SUMMARY_DIR, "source_eval_by_label.csv"),
        eval_summary,
        [
            "tag", "mode", "gamma_label", "alpha_label", "n", "mean_bias",
            "mean_abs_bias", "median_abs_bias", "sd_bias",
            "mean_mu_pred_bias", "mean_delta", "mean_treated_ess",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "current_vs_best.csv"),
        current_vs_best,
        [
            "tag", "mode", "current_mean_bias", "current_mean_abs_bias",
            "best_gamma_label", "best_alpha_label", "best_mean_bias",
            "best_mean_abs_bias", "best_mean_delta", "best_mean_mu_pred_bias",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "oracle_by_variant.csv"),
        oracle_summary,
        ["tag", "mode", "variant", "n", "mean_bias", "mean_abs_bias", "median_abs_bias", "mean_delta"],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "run_summary.csv"),
        run_summary,
        [
            "tag", "mode", "estimate", "bias_vs_superpopulation", "se",
            "covered", "target_bias_vs_superpopulation", "weight_sum",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "nuisance_refit_status.csv"),
        failure_summary,
        ["tag", "mode", "nuisance", "label", "ok", "n", "mean_lambda", "mean_support", "errors"],
    )

    print(f"[nuisance-refit-summary] eval rows: {len(eval_rows)}")
    print(f"[nuisance-refit-summary] oracle rows: {len(oracle_rows)}")
    print(f"[nuisance-refit-summary] run rows: {len(run_rows)}")
    print(f"[nuisance-refit-summary] refit rows: {len(refit_rows)}")
    print(f"[nuisance-refit-summary] wrote {SUMMARY_DIR}")


if __name__ == "__main__":
    main()
