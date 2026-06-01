#!/usr/bin/env python3
"""Summarize C2 score/moment audit outputs.

This script is diagnosis-only. It reads the CSV files produced by
diagnosis/c2/c2_score_moment_audit.R and writes compact summary tables that
make it easier to audit where C2 source-assisted bias enters.
"""

import csv
import glob
import math
import os
import re
from collections import defaultdict


OUTPUT_ROOT = os.environ.get(
    "C2_SCORE_AUDIT_OUTPUT_ROOT",
    os.path.join("diagnosis", "c2", "score_moment_audit"),
)
RESULTS_DIR = os.path.join(OUTPUT_ROOT, "results")
SUMMARY_DIR = os.path.join(OUTPUT_ROOT, "summary")


def read_csv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


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


def max_abs(values):
    values = finite(values)
    return max([abs(value) for value in values]) if values else math.nan


def ratio(numerator, denominator, eps=1e-10):
    numerator = float(numerator)
    denominator = float(denominator)
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        return math.nan
    if abs(denominator) <= eps:
        return math.nan
    return numerator / denominator


def parse_metadata(path, row=None):
    text = os.path.basename(path)
    if row is not None and row.get("tag"):
        text = row.get("tag") + "_" + text

    seed_match = re.search(r"_seed(\d+)", text)
    meta = {
        "K": "",
        "p": "",
        "seed": int(seed_match.group(1)) if seed_match else "",
        "kf": "",
    }

    file_match = re.search(r"_K(\d+)_p(\d+)_seed(\d+)_kf(\d+)", text)
    if file_match:
        meta["K"] = int(file_match.group(1))
        meta["p"] = int(file_match.group(2))
        meta["seed"] = int(file_match.group(3))
        meta["kf"] = int(file_match.group(4))
        return meta

    tag_match = re.search(r"_k(\d+)p(\d+)_seed(\d+)", text)
    if tag_match:
        meta["K"] = int(tag_match.group(1))
        meta["p"] = int(tag_match.group(2))
        meta["seed"] = int(tag_match.group(3))
    return meta


def read_with_metadata(pattern):
    rows = []
    for path in sorted(glob.glob(os.path.join(RESULTS_DIR, pattern))):
        for row in read_csv(path):
            row = dict(row)
            row.update(parse_metadata(path, row))
            row["source_file"] = os.path.basename(path)
            rows.append(row)
    return rows


def keep_keys(row, fields):
    return {field: row.get(field, "") for field in fields}


def group_rows(rows, fields):
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row.get(field, "") for field in fields)].append(row)
    return grouped


def summarize_merged(merged_rows):
    grouped = group_rows(merged_rows, ["K", "mode"])
    out = []
    for key, rows in sorted(grouped.items()):
        K, mode = key
        observed_fill = [
            ratio(f(row, "mean_observed_delta"), f(row, "mean_required_delta"))
            for row in rows
        ]
        model_fill = [
            ratio(f(row, "mean_model_delta_true_pi"), f(row, "mean_required_delta"))
            for row in rows
        ]
        covered = [
            1.0 if abs(f(row, "bias_vs_superpopulation")) <= 1.96 * f(row, "se") else 0.0
            for row in rows if math.isfinite(f(row, "bias_vs_superpopulation")) and math.isfinite(f(row, "se"))
        ]
        out.append({
            "K": K,
            "mode": mode,
            "n_runs": len(rows),
            "mean_bias": mean([f(row, "bias_vs_superpopulation") for row in rows]),
            "mean_abs_bias": mean([abs(f(row, "bias_vs_superpopulation")) for row in rows]),
            "coverage_rate": mean(covered),
            "mean_se": mean([f(row, "se") for row in rows]),
            "mean_target_bias": mean([f(row, "target_bias_vs_superpopulation") for row in rows]),
            "mean_weight_sum": mean([f(row, "weight_sum") for row in rows]),
            "mean_required_delta": mean([f(row, "mean_required_delta") for row in rows]),
            "mean_observed_delta": mean([f(row, "mean_observed_delta") for row in rows]),
            "mean_model_delta_true_pi": mean([f(row, "mean_model_delta_true_pi") for row in rows]),
            "mean_observed_fill_ratio": mean(observed_fill),
            "mean_model_fill_ratio": mean(model_fill),
            "mean_source_bias_observed": mean([f(row, "mean_source_bias_observed") for row in rows]),
            "mean_source_bias_model_true_pi": mean([f(row, "mean_source_bias_model_true_pi") for row in rows]),
            "mean_gamma_linf": mean([f(row, "mean_gamma_linf") for row in rows]),
            "mean_alpha_init_linf": mean([f(row, "mean_alpha_init_linf") for row in rows]),
            "mean_alpha_final_gamma_linf": mean([f(row, "mean_alpha_final_gamma_linf") for row in rows]),
            "mean_intercept_gap_true_pi": mean([f(row, "mean_intercept_gap_true_pi") for row in rows]),
        })
    return out


def summarize_fold_sources(fold_rows):
    grouped = group_rows(fold_rows, ["K", "mode", "source"])
    out = []
    for key, rows in sorted(grouped.items()):
        K, mode, source = key
        out.append({
            "K": K,
            "mode": mode,
            "source": source,
            "n_fold_rows": len(rows),
            "mean_required_delta": mean([f(row, "required_delta") for row in rows]),
            "mean_observed_delta": mean([f(row, "observed_delta") for row in rows]),
            "mean_model_delta_true_pi": mean([f(row, "model_delta_true_pi") for row in rows]),
            "mean_source_bias_observed": mean([f(row, "source_bias_observed") for row in rows]),
            "mean_source_bias_model_true_pi": mean([f(row, "source_bias_model_true_pi") for row in rows]),
            "mean_outcome_noise_delta": mean([f(row, "outcome_noise_delta") for row in rows]),
            "mean_intercept_gap_true_pi": mean([f(row, "intercept_gap_true_pi") for row in rows]),
            "mean_weight_mean_treated": mean([f(row, "weight_mean_treated") for row in rows]),
            "max_weight_max": max_abs([f(row, "weight_max") for row in rows]),
            "mean_gamma_linf": mean([f(row, "gamma_score_linf") for row in rows]),
            "mean_alpha_init_linf": mean([f(row, "alpha_init_score_linf") for row in rows]),
            "mean_alpha_final_gamma_linf": mean([f(row, "alpha_final_gamma_score_linf") for row in rows]),
        })
    return out


def root_cause_flags(merged_rows):
    out = []
    for row in sorted(merged_rows, key=lambda r: (r.get("K", ""), r.get("seed", ""), r.get("mode", ""))):
        observed_bias = f(row, "mean_source_bias_observed")
        model_bias = f(row, "mean_source_bias_model_true_pi")
        noise_bias = observed_bias - model_bias
        if math.isfinite(observed_bias) and abs(observed_bias) > 1e-10:
            model_share = model_bias / observed_bias
            noise_share = noise_bias / observed_bias
        else:
            model_share = math.nan
            noise_share = math.nan

        if math.isfinite(model_share) and abs(model_share) >= 0.5:
            primary_gap = "model_residual_underfill"
        elif math.isfinite(noise_share) and abs(noise_share) >= 0.5:
            primary_gap = "observed_outcome_noise"
        else:
            primary_gap = "mixed_or_small"

        out.append({
            "K": row.get("K", ""),
            "p": row.get("p", ""),
            "seed": row.get("seed", ""),
            "mode": row.get("mode", ""),
            "estimate_bias": f(row, "bias_vs_superpopulation"),
            "se": f(row, "se"),
            "covered_95": abs(f(row, "bias_vs_superpopulation")) <= 1.96 * f(row, "se"),
            "target_bias": f(row, "target_bias_vs_superpopulation"),
            "required_delta": f(row, "mean_required_delta"),
            "observed_delta": f(row, "mean_observed_delta"),
            "model_delta_true_pi": f(row, "mean_model_delta_true_pi"),
            "observed_fill_ratio": ratio(f(row, "mean_observed_delta"), f(row, "mean_required_delta")),
            "model_fill_ratio": ratio(f(row, "mean_model_delta_true_pi"), f(row, "mean_required_delta")),
            "source_bias_observed": observed_bias,
            "source_bias_model_true_pi": model_bias,
            "outcome_noise_share": noise_share,
            "model_gap_share": model_share,
            "mean_gamma_linf": f(row, "mean_gamma_linf"),
            "mean_alpha_init_linf": f(row, "mean_alpha_init_linf"),
            "mean_alpha_final_gamma_linf": f(row, "mean_alpha_final_gamma_linf"),
            "mean_intercept_gap_true_pi": f(row, "mean_intercept_gap_true_pi"),
            "primary_gap": primary_gap,
        })
    return out


def fold_extremes(fold_rows, n=50):
    rows = sorted(
        fold_rows,
        key=lambda row: abs(f(row, "source_bias_observed")) if math.isfinite(f(row, "source_bias_observed")) else -1,
        reverse=True,
    )
    fields = [
        "K", "p", "seed", "kf", "tag", "mode", "fold", "source",
        "required_delta", "observed_delta", "model_delta_true_pi",
        "source_bias_observed", "source_bias_model_true_pi",
        "outcome_noise_delta", "intercept_gap_true_pi", "weight_mean_treated",
        "weight_max", "gamma_score_linf", "alpha_init_score_linf",
        "alpha_final_gamma_score_linf",
    ]
    return [keep_keys(row, fields) for row in rows[:n]]


def main():
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    run_rows = read_with_metadata("*_run_summary.csv")
    score_rows = read_with_metadata("*_score_summary.csv")
    fold_rows = read_with_metadata("*_fold_source_scores.csv")

    if not run_rows or not score_rows or not fold_rows:
        raise SystemExit(
            "missing score audit inputs under {0}: runs={1}, scores={2}, folds={3}".format(
                RESULTS_DIR, len(run_rows), len(score_rows), len(fold_rows)
            )
        )

    score_by_key = {
        (row.get("tag", ""), row.get("mode", "")): row for row in score_rows
    }
    merged_rows = []
    for row in run_rows:
        key = (row.get("tag", ""), row.get("mode", ""))
        merged = dict(row)
        for score_key, score_value in score_by_key.get(key, {}).items():
            if score_key not in merged:
                merged[score_key] = score_value
        merged_rows.append(merged)

    run_fields = [
        "K", "p", "seed", "kf", "tag", "mode", "estimate",
        "truth_superpopulation", "bias_vs_superpopulation", "se",
        "target_estimate", "target_bias_vs_superpopulation", "weight_sum",
        "source_file",
    ]
    score_fields = [
        "K", "p", "seed", "kf", "tag", "mode", "n_fold_source",
        "mean_gamma_intercept_abs", "mean_gamma_linf",
        "mean_alpha_init_intercept_abs", "mean_alpha_init_linf",
        "mean_alpha_final_gamma_intercept_abs", "mean_alpha_final_gamma_linf",
        "mean_required_delta", "mean_observed_delta", "mean_model_delta_true_pi",
        "mean_source_bias_observed", "mean_source_bias_model_true_pi",
        "mean_intercept_gap_true_pi", "source_file",
    ]
    merged_fields = run_fields[:-1] + score_fields[6:-1] + ["source_file"]
    by_mode_fields = [
        "K", "mode", "n_runs", "mean_bias", "mean_abs_bias", "coverage_rate",
        "mean_se", "mean_target_bias", "mean_weight_sum",
        "mean_required_delta", "mean_observed_delta", "mean_model_delta_true_pi",
        "mean_observed_fill_ratio", "mean_model_fill_ratio",
        "mean_source_bias_observed", "mean_source_bias_model_true_pi",
        "mean_gamma_linf", "mean_alpha_init_linf",
        "mean_alpha_final_gamma_linf", "mean_intercept_gap_true_pi",
    ]
    fold_source_fields = [
        "K", "mode", "source", "n_fold_rows", "mean_required_delta",
        "mean_observed_delta", "mean_model_delta_true_pi",
        "mean_source_bias_observed", "mean_source_bias_model_true_pi",
        "mean_outcome_noise_delta", "mean_intercept_gap_true_pi",
        "mean_weight_mean_treated", "max_weight_max",
        "mean_gamma_linf", "mean_alpha_init_linf",
        "mean_alpha_final_gamma_linf",
    ]
    flags_fields = [
        "K", "p", "seed", "mode", "estimate_bias", "se", "covered_95",
        "target_bias", "required_delta", "observed_delta",
        "model_delta_true_pi", "observed_fill_ratio", "model_fill_ratio",
        "source_bias_observed", "source_bias_model_true_pi",
        "outcome_noise_share", "model_gap_share", "mean_gamma_linf",
        "mean_alpha_init_linf", "mean_alpha_final_gamma_linf",
        "mean_intercept_gap_true_pi", "primary_gap",
    ]
    extremes_fields = [
        "K", "p", "seed", "kf", "tag", "mode", "fold", "source",
        "required_delta", "observed_delta", "model_delta_true_pi",
        "source_bias_observed", "source_bias_model_true_pi",
        "outcome_noise_delta", "intercept_gap_true_pi", "weight_mean_treated",
        "weight_max", "gamma_score_linf", "alpha_init_score_linf",
        "alpha_final_gamma_score_linf",
    ]

    write_csv(os.path.join(SUMMARY_DIR, "all_run_summary.csv"), run_rows, run_fields)
    write_csv(os.path.join(SUMMARY_DIR, "all_score_summary.csv"), score_rows, score_fields)
    write_csv(os.path.join(SUMMARY_DIR, "merged_run_score_summary.csv"), merged_rows, merged_fields)
    write_csv(os.path.join(SUMMARY_DIR, "bias_and_moment_by_mode.csv"), summarize_merged(merged_rows), by_mode_fields)
    write_csv(os.path.join(SUMMARY_DIR, "fold_source_by_source.csv"), summarize_fold_sources(fold_rows), fold_source_fields)
    write_csv(os.path.join(SUMMARY_DIR, "root_cause_flags.csv"), root_cause_flags(merged_rows), flags_fields)
    write_csv(os.path.join(SUMMARY_DIR, "fold_source_extremes.csv"), fold_extremes(fold_rows), extremes_fields)

    print("[score-summary] input run rows: {0}".format(len(run_rows)))
    print("[score-summary] input score rows: {0}".format(len(score_rows)))
    print("[score-summary] input fold/source rows: {0}".format(len(fold_rows)))
    print("[score-summary] wrote summary tables to {0}".format(SUMMARY_DIR))


if __name__ == "__main__":
    main()
