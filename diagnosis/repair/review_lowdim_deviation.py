#!/usr/bin/env python3
"""Build the low-dimensional deviation report after campaign controllers exit."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    base, output = pilot.scratch_path(args.run_root), pilot.scratch_path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    campaigns = dict(deviation="lowdim_deviation_mc50_cutoff2_v1",
                     control="lowdim_cutoff2_control_mc50_v1",
                     baselines="lowdim_deviation_baselines_mc50_v1")
    counts = []
    for label, name in campaigns.items():
        root = base / name
        script = "summarize_repeat_pilot.py" if label == "baselines" else "diagnose_deviation_results.py"
        result = subprocess.run([sys.executable, "-u", "-B", str(Path(__file__).with_name(script)),
                                 str(root), str(output / label)])
        if result.returncode != 0 and not (label == "baselines" and result.returncode == 1):
            raise RuntimeError("Report failed for " + label)
        with (output / label / "task_status.csv").open() as stream:
            tasks = list(csv.DictReader(stream))
        counts.append(dict(campaign=label, planned=len(tasks),
                           complete=sum(row["state"] == "COMPLETE" for row in tasks)))
    complete = all(row["planned"] == row["complete"] for row in counts)
    pilot.replace_json(output / "report_status.json", dict(
        state="COMPLETE" if complete else "PARTIAL", campaigns=counts))
    lines = ["# Low-dimensional source-deviation validation", "",
             "n=1000/site, p=10/20/50, C1/C2/C3, 50 prescribed seeds per cell, cutoff=2.",
             "Both control and deviation runs use the same three-level calibrated algorithm.",
             "Population truth is fixed within each cell. No seeds are excluded by statistical performance.", ""]
    lines += ["{}: {}/{} committed repeats.".format(row["campaign"], row["complete"], row["planned"])
              for row in counts]
    if not complete:
        lines += ["", "PARTIAL: these results do not represent all planned repeats."]
    lines += ["", "| p | Scenario | rho | N | Bias | RMSE | Coverage | Mean SE / SD |",
              "|---:|---|---:|---:|---:|---:|---:|---:|"]
    for label in ("control", "deviation"):
        with (output / label / "metrics.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        for row in rows:
            if not row["method"].endswith("_joint_tate") or int(row["n_success"]) == 0:
                continue
            ratio = "{:.3f}".format(float(row["se_to_sd_ratio"])) if row["se_to_sd_ratio"] else "NA"
            lines.append("| {} | {} | {} | {} | {:.5f} | {:.5f} | {:.1%} | {} |".format(
                row["p"], row["config"], row["rho"], row["n_success"],
                float(row["bias"]), float(row["rmse"]), float(row["coverage"]), ratio))
    lines += ["", "Interpretation checks:", "",
              "- Coverage intervals and Monte Carlo uncertainty are in each metrics.csv.",
              "- Source-weight diagnostics distinguish treated/control weights, penalty activation, residual bias and target-anchor error.",
              "- Cutoff 2 activates a soft penalty; exceeding 2 does not force zero source weight.",
              "- Baselines use matching data/seed identities; their full results are in baselines/metrics.csv.",
              "- This report does not infer a coverage guarantee from a small number of repeats."]
    (output / "report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
