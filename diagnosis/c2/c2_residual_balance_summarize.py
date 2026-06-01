#!/usr/bin/env python3
"""Summarize C2 residual-balance probe outputs on MSI."""

import csv
import glob
import math
import os
import statistics
from collections import defaultdict


OUTPUT_ROOT = os.environ.get(
    "C2_RBAL_OUTPUT_ROOT",
    os.path.join("diagnosis", "c2", "residual_balance_probe"),
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
    return [x for x in values if math.isfinite(x)]


def mean(values):
    vals = finite(values)
    return sum(vals) / len(vals) if vals else math.nan


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


def summarize_balance(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[
            (
                row.get("tag", ""),
                row.get("mode", ""),
                row.get("gamma_label", ""),
                row.get("alpha_label", ""),
                row.get("gamma_basis", ""),
                row.get("alpha_basis", ""),
            )
        ].append(row)

    out = []
    for key, sub in sorted(grouped.items()):
        tag, mode, gamma_label, alpha_label, gamma_basis, alpha_basis = key
        source_bias = [f(row, "source_bias_observed") for row in sub]
        model_bias_a = [f(row, "source_bias_model_observed_A") for row in sub]
        model_bias_pi = [f(row, "source_bias_model_true_pi") for row in sub]
        required = [f(row, "required_delta") for row in sub]
        observed = [f(row, "observed_delta") for row in sub]
        model_a = [f(row, "model_delta_observed_A") for row in sub]
        model_pi = [f(row, "model_delta_true_pi") for row in sub]
        noise = [f(row, "outcome_noise_delta") for row in sub]
        intercept_a = [f(row, "intercept_gap_observed_A") for row in sub]
        intercept_pi = [f(row, "intercept_gap_true_pi") for row in sub]
        raw_z_a = [f(row, "raw_z_gap_observed_A") for row in sub]
        raw_z_pi = [f(row, "raw_z_gap_true_pi") for row in sub]
        norm_z_a = [f(row, "norm_z_gap_observed_A") for row in sub]
        norm_z_pi = [f(row, "norm_z_gap_true_pi") for row in sub]
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "gamma_label": gamma_label,
                "alpha_label": alpha_label,
                "gamma_basis": gamma_basis,
                "alpha_basis": alpha_basis,
                "n": len(sub),
                "mean_source_bias_observed": mean(source_bias),
                "mean_abs_source_bias_observed": mean([abs(x) for x in source_bias]),
                "median_abs_source_bias_observed": median([abs(x) for x in source_bias]),
                "mean_source_bias_model_observed_A": mean(model_bias_a),
                "mean_source_bias_model_true_pi": mean(model_bias_pi),
                "mean_required_delta": mean(required),
                "mean_observed_delta": mean(observed),
                "mean_model_delta_observed_A": mean(model_a),
                "mean_model_delta_true_pi": mean(model_pi),
                "mean_outcome_noise_delta": mean(noise),
                "mean_intercept_gap_observed_A": mean(intercept_a),
                "mean_intercept_gap_true_pi": mean(intercept_pi),
                "mean_raw_z_gap_observed_A": mean(raw_z_a),
                "mean_raw_z_gap_true_pi": mean(raw_z_pi),
                "mean_norm_z_gap_observed_A": mean(norm_z_a),
                "mean_norm_z_gap_true_pi": mean(norm_z_pi),
            }
        )
    return out


def summarize_current_vs_best(summary_rows):
    grouped = defaultdict(list)
    for row in summary_rows:
        grouped[(row["tag"], row["mode"])].append(row)

    out = []
    for key, sub in sorted(grouped.items()):
        tag, mode = key
        current = [
            row for row in sub
            if row["gamma_label"] == "current"
            and row["alpha_label"] == "current"
            and row["gamma_basis"] == "fitted"
            and row["alpha_basis"] == "fitted"
        ]
        best = min(sub, key=lambda row: abs(float(row["mean_source_bias_observed"])))
        cur = current[0] if current else {}
        out.append(
            {
                "tag": tag,
                "mode": mode,
                "current_source_bias": cur.get("mean_source_bias_observed", math.nan),
                "current_required_delta": cur.get("mean_required_delta", math.nan),
                "current_observed_delta": cur.get("mean_observed_delta", math.nan),
                "current_model_delta_true_pi": cur.get("mean_model_delta_true_pi", math.nan),
                "current_intercept_gap_true_pi": cur.get("mean_intercept_gap_true_pi", math.nan),
                "current_raw_z_gap_true_pi": cur.get("mean_raw_z_gap_true_pi", math.nan),
                "best_gamma_label": best["gamma_label"],
                "best_alpha_label": best["alpha_label"],
                "best_gamma_basis": best["gamma_basis"],
                "best_alpha_basis": best["alpha_basis"],
                "best_source_bias": best["mean_source_bias_observed"],
                "best_required_delta": best["mean_required_delta"],
                "best_observed_delta": best["mean_observed_delta"],
                "best_model_delta_true_pi": best["mean_model_delta_true_pi"],
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
                "target_bias_vs_superpopulation": f(row, "target_bias_vs_superpopulation"),
                "weight_sum": f(row, "weight_sum"),
            }
        )
    return out


def main():
    balance_rows = read_csvs("*_residual_balance.csv")
    run_rows = read_csvs("*_run_summary.csv")
    error_rows = read_csvs("*_refit_errors.csv")
    if not balance_rows:
        raise SystemExit("No residual_balance CSVs found.")

    balance_summary = summarize_balance(balance_rows)
    current_vs_best = summarize_current_vs_best(balance_summary)
    run_summary = summarize_runs(run_rows)

    write_csv(
        os.path.join(SUMMARY_DIR, "residual_balance_by_label.csv"),
        balance_summary,
        [
            "tag", "mode", "gamma_label", "alpha_label", "gamma_basis",
            "alpha_basis", "n", "mean_source_bias_observed",
            "mean_abs_source_bias_observed", "median_abs_source_bias_observed",
            "mean_source_bias_model_observed_A",
            "mean_source_bias_model_true_pi", "mean_required_delta",
            "mean_observed_delta", "mean_model_delta_observed_A",
            "mean_model_delta_true_pi", "mean_outcome_noise_delta",
            "mean_intercept_gap_observed_A", "mean_intercept_gap_true_pi",
            "mean_raw_z_gap_observed_A", "mean_raw_z_gap_true_pi",
            "mean_norm_z_gap_observed_A", "mean_norm_z_gap_true_pi",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "current_vs_best.csv"),
        current_vs_best,
        [
            "tag", "mode", "current_source_bias", "current_required_delta",
            "current_observed_delta", "current_model_delta_true_pi",
            "current_intercept_gap_true_pi", "current_raw_z_gap_true_pi",
            "best_gamma_label", "best_alpha_label", "best_gamma_basis",
            "best_alpha_basis", "best_source_bias", "best_required_delta",
            "best_observed_delta", "best_model_delta_true_pi",
        ],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "run_summary.csv"),
        run_summary,
        ["tag", "mode", "estimate", "bias_vs_superpopulation", "se",
         "target_bias_vs_superpopulation", "weight_sum"],
    )
    write_csv(
        os.path.join(SUMMARY_DIR, "refit_errors.csv"),
        error_rows,
        ["tag", "mode", "fold", "source", "nuisance", "label", "error"],
    )

    print("[residual-balance-summary] balance rows: %d" % len(balance_rows))
    print("[residual-balance-summary] run rows: %d" % len(run_rows))
    print("[residual-balance-summary] refit error rows: %d" % len(error_rows))
    print("[residual-balance-summary] wrote %s" % SUMMARY_DIR)


if __name__ == "__main__":
    main()
