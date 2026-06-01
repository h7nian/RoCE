#!/usr/bin/env python3
"""Compare two C2 score/moment audit summary roots.

The intended use is cached-vs-no-cache nuisance diagnostics. It reads
summary/root_cause_flags.csv from each root and writes paired deltas.
"""

import csv
import math
import os


BASE_ROOT = os.environ.get(
    "C2_SCORE_COMPARE_BASE_ROOT",
    os.path.join("diagnosis", "c2", "score_moment_audit"),
)
CANDIDATE_ROOT = os.environ.get(
    "C2_SCORE_COMPARE_CANDIDATE_ROOT",
    os.path.join("diagnosis", "c2", "score_moment_audit_nocache"),
)
OUTPUT_ROOT = os.environ.get(
    "C2_SCORE_COMPARE_OUTPUT_ROOT",
    os.path.join("diagnosis", "c2", "score_moment_compare"),
)


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


def key(row):
    return (row.get("K", ""), row.get("p", ""), row.get("seed", ""), row.get("mode", ""))


def delta_row(base, candidate):
    out = {
        "K": candidate.get("K", ""),
        "p": candidate.get("p", ""),
        "seed": candidate.get("seed", ""),
        "mode": candidate.get("mode", ""),
        "base_primary_gap": base.get("primary_gap", ""),
        "candidate_primary_gap": candidate.get("primary_gap", ""),
    }
    for field in (
        "estimate_bias", "target_bias", "required_delta", "observed_delta",
        "model_delta_true_pi", "observed_fill_ratio", "model_fill_ratio",
        "source_bias_observed", "source_bias_model_true_pi",
        "mean_gamma_linf", "mean_alpha_init_linf",
        "mean_alpha_final_gamma_linf", "mean_intercept_gap_true_pi",
    ):
        out["base_" + field] = f(base, field)
        out["candidate_" + field] = f(candidate, field)
        out["delta_" + field] = f(candidate, field) - f(base, field)
    return out


def summarize_deltas(rows):
    grouped = {}
    for row in rows:
        group_key = (row.get("K", ""), row.get("mode", ""))
        grouped.setdefault(group_key, []).append(row)

    out = []
    for group_key, sub in sorted(grouped.items()):
        K, mode = group_key
        out.append({
            "K": K,
            "mode": mode,
            "n_pairs": len(sub),
            "mean_delta_estimate_bias": mean([f(row, "delta_estimate_bias") for row in sub]),
            "mean_delta_source_bias_observed": mean([f(row, "delta_source_bias_observed") for row in sub]),
            "mean_delta_source_bias_model_true_pi": mean([f(row, "delta_source_bias_model_true_pi") for row in sub]),
            "mean_delta_observed_delta": mean([f(row, "delta_observed_delta") for row in sub]),
            "mean_delta_model_delta_true_pi": mean([f(row, "delta_model_delta_true_pi") for row in sub]),
            "mean_delta_gamma_linf": mean([f(row, "delta_mean_gamma_linf") for row in sub]),
            "mean_delta_alpha_init_linf": mean([f(row, "delta_mean_alpha_init_linf") for row in sub]),
        })
    return out


def main():
    base_path = os.path.join(BASE_ROOT, "summary", "root_cause_flags.csv")
    candidate_path = os.path.join(CANDIDATE_ROOT, "summary", "root_cause_flags.csv")
    if not os.path.exists(base_path):
        raise SystemExit("missing base summary: {0}".format(base_path))
    if not os.path.exists(candidate_path):
        raise SystemExit("missing candidate summary: {0}".format(candidate_path))

    base_rows = {key(row): row for row in read_csv(base_path)}
    candidate_rows = {key(row): row for row in read_csv(candidate_path)}
    paired_keys = sorted(set(base_rows).intersection(candidate_rows))
    if not paired_keys:
        raise SystemExit("no overlapping rows between {0} and {1}".format(base_path, candidate_path))

    rows = [delta_row(base_rows[item], candidate_rows[item]) for item in paired_keys]
    fields = [
        "K", "p", "seed", "mode", "base_primary_gap", "candidate_primary_gap",
        "base_estimate_bias", "candidate_estimate_bias", "delta_estimate_bias",
        "base_target_bias", "candidate_target_bias", "delta_target_bias",
        "base_required_delta", "candidate_required_delta", "delta_required_delta",
        "base_observed_delta", "candidate_observed_delta", "delta_observed_delta",
        "base_model_delta_true_pi", "candidate_model_delta_true_pi", "delta_model_delta_true_pi",
        "base_observed_fill_ratio", "candidate_observed_fill_ratio", "delta_observed_fill_ratio",
        "base_model_fill_ratio", "candidate_model_fill_ratio", "delta_model_fill_ratio",
        "base_source_bias_observed", "candidate_source_bias_observed", "delta_source_bias_observed",
        "base_source_bias_model_true_pi", "candidate_source_bias_model_true_pi", "delta_source_bias_model_true_pi",
        "base_mean_gamma_linf", "candidate_mean_gamma_linf", "delta_mean_gamma_linf",
        "base_mean_alpha_init_linf", "candidate_mean_alpha_init_linf", "delta_mean_alpha_init_linf",
        "base_mean_alpha_final_gamma_linf", "candidate_mean_alpha_final_gamma_linf", "delta_mean_alpha_final_gamma_linf",
        "base_mean_intercept_gap_true_pi", "candidate_mean_intercept_gap_true_pi", "delta_mean_intercept_gap_true_pi",
    ]
    summary_fields = [
        "K", "mode", "n_pairs", "mean_delta_estimate_bias",
        "mean_delta_source_bias_observed",
        "mean_delta_source_bias_model_true_pi",
        "mean_delta_observed_delta", "mean_delta_model_delta_true_pi",
        "mean_delta_gamma_linf", "mean_delta_alpha_init_linf",
    ]
    write_csv(os.path.join(OUTPUT_ROOT, "summary", "root_cause_flag_delta.csv"), rows, fields)
    write_csv(os.path.join(OUTPUT_ROOT, "summary", "delta_by_mode.csv"), summarize_deltas(rows), summary_fields)
    print("[score-compare] paired rows: {0}".format(len(rows)))
    print("[score-compare] wrote comparison under {0}".format(os.path.join(OUTPUT_ROOT, "summary")))


if __name__ == "__main__":
    main()
