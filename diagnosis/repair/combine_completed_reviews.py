#!/usr/bin/env python3
"""Combine disjoint seed panels from completed, compatible paired reviews."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import review_completed_highdim as review
from summarize_repeat_pilot import write_csv

PROVENANCE_FIELDS = {"first_seed", "exclude_nodes", "baseline_reference_roots"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("reports", type=Path, nargs="+")
    args = parser.parse_args()
    output = pilot.scratch_path(args.output)
    if output.exists(): raise ValueError("Use a new output directory")
    records, provenance, configurations = [], [], None
    identities, repeats = set(), set()
    for report in args.reports:
        report = pilot.scratch_path(report)
        checks = json.loads((report / "pairing_checks.json").read_text())
        if checks["repeats"] != checks["data_hash_matches"] or checks["repeats"] != checks["target_reference_matches"]:
            raise ValueError("Incomplete pairing checks")
        signature = []
        for field in ("method_root", "baseline_root"):
            config = pilot.verify(pilot.scratch_path(checks[field]))
            # Baseline reference paths identify disjoint seed banks. Every
            # reference value was checked against its own paired method run.
            signature.append({key:value for key,value in config.items() if key not in PROVENANCE_FIELDS})
        if configurations is not None and configurations != signature:
            raise ValueError("Recorded configurations differ beyond seed/reference provenance and node exclusions")
        configurations = signature
        path = report / "paired_estimates.csv"
        with path.open() as stream:
            panel = list(csv.DictReader(stream))
        panel_repeats = set()
        for row in panel:
            for field in ("estimate","se","truth","target_error"):
                row[field] = float(row[field])
                if not math.isfinite(row[field]): raise ValueError("Nonfinite paired estimate")
            if row["se"] <= 0: raise ValueError("Invalid standard error")
            key = review.identity(row)+(int(row["n_deviated_sites"]),)
            identity = key+(row["method"],)
            if identity in identities: raise ValueError("Overlapping seed/method records")
            identities.add(identity);panel_repeats.add(key)
        if len(panel_repeats) != checks["repeats"] or panel_repeats & repeats:
            raise ValueError("Repeat count mismatch or overlapping panels")
        repeats.update(panel_repeats);records.extend(panel)
        provenance.append(dict(report=str(report),repeats=len(panel_repeats),
            estimates_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),checks=checks))
    metrics = review.summarize_paired_estimates(records)
    output.mkdir(parents=True,exist_ok=False)
    write_csv(output / "paired_estimates.csv", records)
    write_csv(output / "metrics.csv", metrics)
    (output / "combined_provenance.json").write_text(json.dumps(dict(repeats=len(repeats),panels=provenance),indent=2)+"\n")
    (output / "README.md").write_text(
        "# Combined completed seed panels\n\n"
        "{} disjoint repeats from compatible method/baseline configurations.\n"
        "Cell-specific repeat counts and Monte Carlo uncertainty are in metrics.csv.\n"
        "No seed is selected by its estimate, coverage or error. Correlated scenarios remain separate.\n".format(len(repeats)))
    print("Combined {} verified, nonoverlapping repeats".format(len(repeats)))


if __name__ == "__main__":
    main()
