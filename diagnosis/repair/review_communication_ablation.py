#!/usr/bin/env python3
"""Compare completed one/two-round protocols on the same generated data."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
import review_completed_highdim as completed
import submit_repeat_pilot as pilot
from summarize_repeat_pilot import write_csv


DESIGN_FIELDS = (
    "dgp_type", "dgp_version", "dgp_configuration", "working_dimension",
    "dgp_max_signal_slopes", "truth_method", "outcome_family",
    "aggregation_lambda", "aggregation_cutoff", "nlambda_init",
    "nuisance_lambda_rule", "nuisance_training_policy", "target_nuisance_method",
    "source_validation_method", "calibration_layout", "nuisance_tol",
    "calibration_recipe", "target_propensity_initialization",
    "target_training_radius", "target_inference_radius", "roce_crossfit_levels",
    "M_tau", "M_tau_inference", "deviation_mechanism", "compiled_nuisance_solver",
)
MODES = ("joint_tate", "separate_arms", "common_tate")


def index_tasks(roots, protocol):
    tasks = {}
    provenance = []
    for root in roots:
        root = pilot.scratch_path(root)
        configuration = pilot.verify(root)
        for task in completed.manifest(root):
            if task["protocol"] != protocol:
                raise ValueError("Unexpected communication protocol: " + str(root))
            key = completed.identity(task)
            if key in tasks:
                raise ValueError("Overlapping seed panels: " + str(key))
            tasks[key] = (root, task, configuration)
        provenance.append(dict(root=str(root), manifest_sha256=hashlib.sha256(
            (root / "manifest.csv").read_bytes()).hexdigest()))
    return tasks, provenance


def pair_protocols(candidate, reference):
    candidate_root, task, configuration = candidate
    reference_root, reference_task, reference_configuration = reference
    for field in ("library", "n_per_site", "n_folds", "n_deviated_sites"):
        default = 1 if field == "n_deviated_sites" else None
        if configuration.get(field, default) != reference_configuration.get(field, default):
            raise ValueError("Protocol configurations differ: " + field)
    two, two_hash = completed.completed_result(candidate_root, task)
    one, one_hash = completed.completed_result(reference_root, reference_task)
    if two_hash != one_hash:
        raise ValueError("Generated data differ: " + str(completed.identity(task)))
    for method in ("target_only_ate", "target_anchor_ate"):
        for field in ("estimate", "se", "truth"):
            if abs(float(two[method][field]) - float(one[method][field])) > 1e-12:
                raise ValueError("Target reference differs: " + method + ":" + field)
    records = []
    for mode in MODES:
        suffix = "_crossfit_ate" + ("" if mode == "common_tate" else "_" + mode)
        two_row = two["two_round" + suffix]
        one_row = one["one_round" + suffix]
        for field in DESIGN_FIELDS:
            if field not in two_row or field not in one_row or two_row[field] != one_row[field]:
                raise ValueError("Scientific setting differs: " + field)
        if abs(float(two_row["truth"]) - float(one_row["truth"])) > 1e-12:
            raise ValueError("Population truths differ")
        for protocol, row in (("one_round", one_row), ("two_round", two_row)):
            records.append(dict(task, protocol=protocol, method=protocol + "_" + mode,
                mode=mode, estimate=float(row["estimate"]), se=float(row["se"]),
                truth=float(row["truth"]), data_sha256=two_hash,
                n_deviated_sites=int(configuration.get("n_deviated_sites", 1)),
                target_error=float(two["target_only_ate"]["estimate"]) - float(row["truth"])))
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--candidate-roots", type=Path, nargs="+", required=True)
    parser.add_argument("--reference-roots", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    output = pilot.scratch_path(args.output)
    if output.exists():
        raise ValueError("Use a new report directory")
    candidates, candidate_provenance = index_tasks(args.candidate_roots, "two_round")
    references, reference_provenance = index_tasks(args.reference_roots, "one_round")
    if not candidates or not set(candidates).issubset(references):
        raise ValueError("Empty comparison or missing one-round reference seeds")
    with ThreadPoolExecutor(max_workers=8) as pool:
        pairs = pool.map(lambda key: pair_protocols(candidates[key], references[key]), candidates)
        records = [row for pair in pairs for row in pair]
    metrics = completed.summarize_paired_estimates(records)
    groups = defaultdict(dict)
    for row in records:
        cell = (row["config"], int(row["p"]), int(row["K"]), float(row["rho"]), row["mode"])
        groups[cell].setdefault(int(row["sim_id"]), {})[row["protocol"]] = row
    contrasts = []
    for cell, seeds in sorted(groups.items()):
        differences = []
        for pair in seeds.values():
            errors = {name: row["estimate"] - row["truth"] for name, row in pair.items()}
            differences.append(errors["two_round"]**2 - errors["one_round"]**2)
        contrasts.append(dict(config=cell[0], p=cell[1], K=cell[2], rho=cell[3], mode=cell[4],
            repeats=len(seeds), mse_difference=statistics.mean(differences),
            mse_difference_mcse=statistics.stdev(differences)/math.sqrt(len(differences))))
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "paired_estimates.csv", records)
    write_csv(output / "metrics.csv", metrics)
    write_csv(output / "paired_mse.csv", contrasts)
    (output / "pairing_checks.json").write_text(json.dumps(dict(
        pairs=len(candidates), data_hash_matches=len(candidates),
        ordinary_target_matches=len(candidates), calibrated_target_matches=len(candidates),
        source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        candidate_roots=candidate_provenance, reference_roots=reference_provenance), indent=2) + "\n")
    (output / "README.md").write_text(
        "# Paired communication ablation\n\n"
        "{} complete pairs; all candidate seeds retained. One/two communication rounds\n"
        "both use three-layer cross-fitting. Data hashes, both target anchors and scientific\n"
        "settings match. Initial source outcome models differ by protocol.\n\n"
        "metrics.csv includes all three aggregation modes and coverage uncertainty.\n"
        "paired_mse.csv reports two-round minus one-round MSE with paired Monte Carlo SE.\n"
        "The protocol comparison does not measure the number of cross-fitting layers.\n".format(len(candidates)))
    print("Verified {} paired protocol repeats; saved {}".format(len(candidates), output))


if __name__ == "__main__":
    main()
