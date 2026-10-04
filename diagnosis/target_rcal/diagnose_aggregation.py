#!/usr/bin/env python3
"""Compare recorded screening and evaluation discrepancies; no estimator changes."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("production_root", type=Path)
    parser.add_argument("output_root", type=Path)
    args = parser.parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    records = []
    for family, root in (("negative_transfer", args.production_root),
                         ("shared_shift", args.production_root / "shared_shift")):
        with (root / "manifest_main.csv").open() as handle:
            tasks = [row for row in csv.DictReader(handle)
                     if row["config"] == "C1" and row["K"] == "4" and
                     float(row["rho"]) in (0, 2.5)]
        for task in tasks:
            path = root / "raw" / f"task_{int(task['task_id']):06d}.csv"
            if not path.exists():
                continue
            with path.open() as handle:
                rows = {row["method"]: row for row in csv.DictReader(handle)}
            direct = rows["one_round_crossfit_ate"]
            armwise = rows["one_round_crossfit_ate_armwise"]
            for source in range(1, 5):
                prefix = f"source_s{source}_"
                record = dict(family=family, rho=float(task["rho"]), seed=int(task["sim_id"]),
                              source=source)
                for field in ("weight", "wald_mean", "penalty_activation_fraction",
                              "screening_discrepancy_mean", "evaluation_discrepancy_mean",
                              "screening_discrepancy_se_mean"):
                    record[field] = float(direct[prefix + field])
                for arm in (0, 1):
                    record[f"mu{arm}_weight"] = float(armwise[f"mu{arm}_" + prefix + "weight"])
                    record[f"mu{arm}_screening"] = float(
                        armwise[f"mu{arm}_" + prefix + "screening_discrepancy_mean"])
                    record[f"mu{arm}_evaluation"] = float(
                        armwise[f"mu{arm}_" + prefix + "evaluation_discrepancy_mean"])
                records.append(record)
    with (args.output_root / "source_records.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)
    grouped = defaultdict(list)
    for record in records:
        grouped[(record["family"], record["rho"], record["source"])].append(record)
    summary = []
    fields = [field for field in records[0] if field not in ("family", "rho", "seed", "source")]
    for (family, rho, source), group in sorted(grouped.items()):
        row = dict(family=family, rho=rho, source=source, n=len(group))
        row.update({field: statistics.mean(record[field] for record in group) for field in fields})
        summary.append(row)
    with (args.output_root / "source_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys())
        writer.writeheader()
        writer.writerows(summary)
    for row in summary:
        print("{family:18s} rho={rho:g} s={source} n={n} weight={weight:.3f} "
              "mu1_weight={mu1_weight:.3f} mu0_weight={mu0_weight:.3f} "
              "screen={screening_discrepancy_mean:+.4f} "
              "eval={evaluation_discrepancy_mean:+.4f}".format(**row))


if __name__ == "__main__":
    main()
