#!/usr/bin/env python3
"""Summarize completed C2 diagnosis CSVs.

This script is intentionally read-only. It scans diagnosis/c2/results by
default and reports the quantities needed for the C2 coverage diagnosis:
coverage, center bias, empirical SD, mean SE, bias/SE, SE/SD, and the
normal-approximation coverage after setting bias to zero.
"""

import argparse
import csv
import math
import re
from pathlib import Path


FOCUS_METHODS = {
    "target_only",
    "federated_dr",
    "one_round_crossfit",
    "two_round_crossfit",
    "pooled_dr",
    "tilted_aipw",
    "oracle_dr",
}


def as_float(row, field):
    value = row.get(field, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def centered_normal_coverage(se_mean, bias_sd):
    if not math.isfinite(se_mean) or not math.isfinite(bias_sd) or bias_sd <= 0:
        return float("nan")
    z = 1.96 * se_mean / bias_sd
    return math.erf(z / math.sqrt(2.0))


def iter_summary_rows(results_dir: Path):
    for path in sorted(results_dir.glob("*_summary.csv")):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row.get("method") not in FOCUS_METHODS:
                    continue
                yield path, row


def parse_setting_from_name(path):
    match = re.search(r"_n([0-9]+)_K([0-9]+)_p([0-9]+)_", path.name)
    if not match:
        return {}
    return {
        "n_total": int(match.group(1)),
        "K": int(match.group(2)),
        "p": int(match.group(3)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        default="diagnosis/c2/results",
        help="Directory containing *_summary.csv files.",
    )
    parser.add_argument(
        "--include-kf3",
        action="store_true",
        help="Include non-default K_f=3 sensitivity outputs.",
    )
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    if not results_dir.exists():
        raise SystemExit(f"results directory does not exist: {results_dir}")

    rows = []
    for path, row in iter_summary_rows(results_dir):
        meta = parse_setting_from_name(path)
        n_folds = int(as_float(row, "n_folds")) if row.get("n_folds") else None
        if not args.include_kf3 and n_folds != 10:
            continue
        bias = as_float(row, "bias_mean")
        bias_sd = as_float(row, "bias_sd")
        se_mean = as_float(row, "se_mean")
        coverage = as_float(row, "coverage")
        rows.append(
            {
                "file": path.name,
                "tag": row.get("tag", ""),
                "estimand": row.get("estimand_type", ""),
                "K": int(as_float(row, "K")) if row.get("K") else meta.get("K", -1),
                "p": int(as_float(row, "p")) if row.get("p") else meta.get("p", -1),
                "method": row["method"],
                "bias": bias,
                "bias_sd": bias_sd,
                "se": se_mean,
                "coverage": coverage,
                "bias_over_se": abs(bias) / se_mean if se_mean > 0 else float("nan"),
                "se_over_sd": se_mean / bias_sd if bias_sd > 0 else float("nan"),
                "centered": centered_normal_coverage(se_mean, bias_sd),
                "n": int(as_float(row, "n_success")),
            }
        )

    rows.sort(key=lambda r: (r["tag"], r["estimand"], r["K"], r["p"], r["method"]))
    if not rows:
        print(f"No completed summary rows found in {results_dir}")
        return 0

    current = None
    for row in rows:
        group = (row["tag"], row["estimand"], row["K"], row["p"])
        if group != current:
            current = group
            print()
            print(
                f"### {row['tag']} | estimand={row['estimand']} "
                f"K={row['K']} p={row['p']} n={row['n']}"
            )
            print(
                "method                 bias      sd      se     cov  "
                "|bias|/se  se/sd  centered"
            )
        print(
            f"{row['method']:<20s} "
            f"{row['bias']:+.5f} {row['bias_sd']:.5f} {row['se']:.5f} "
            f"{row['coverage']:.3f} {row['bias_over_se']:.3f} "
            f"{row['se_over_sd']:.3f} {row['centered']:.3f}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
