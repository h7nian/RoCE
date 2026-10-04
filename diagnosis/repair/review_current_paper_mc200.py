#!/usr/bin/env python3
"""Review a complete MC200 panel, including its reused baseline results."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import data_identity, read_campaign
from review_completed_highdim import completed_result, summarize_paired_estimates
from summarize_repeat_pilot import write_csv

DATA_FIELDS = ("config", "K", "p", "rho", "sim_id", "deviation_mechanism", "n_deviated_sites")


def index_pairs(rows, panel, dimension, sources=None):
    pairs = defaultdict(dict)
    for row in rows:
        if row["panel"] != panel or int(row["p"]) != dimension or row["role"] not in ("method", "baseline") or (sources is not None and int(row["K"]) not in sources):
            continue
        key = tuple(row[field] for field in DATA_FIELDS)
        if row["role"] in pairs[key]:
            raise ValueError("Duplicate panel role: " + str(key))
        pairs[key][row["role"]] = row
    if not pairs or any(set(pair) != {"method", "baseline"} for pair in pairs.values()):
        raise ValueError("Every prescribed data setting needs one method and baseline")
    counts = Counter(tuple(value for field, value in zip(DATA_FIELDS, key) if field != "sim_id") for key in pairs)
    if any(count != 200 for count in counts.values()):
        raise ValueError("Every panel cell must contain exactly200 prescribed seeds")
    for pair in pairs.values():
        if not 1 <= int(pair["method"]["sim_id"]) <= 200:
            raise ValueError("The primary validation requires seeds1–200")
    return pairs


def primary_records(task, results, data_hash):
    expected_method = {"one_round_crossfit_ate", "one_round_crossfit_ate_separate_arms",
                       "one_round_crossfit_ate_joint_tate", "target_only_ate", "target_anchor_ate"}
    if not expected_method.issubset(results):
        raise ValueError("A requested primary method is missing")
    truth = float(results["target_only_ate"]["truth"])
    target_error = float(results["target_only_ate"]["estimate"])-truth
    rows = []
    for name, result in results.items():
        if abs(float(result["truth"])-truth) > 1e-12:
            raise ValueError("Population truths differ between methods")
        rows.append(dict(task, method=name, estimate=float(result["estimate"]),
            se=float(result["se"]), truth=truth, target_error=target_error, data_sha256=data_hash))
    return rows


def combine_pair(task, method, baseline, method_hash, baseline_hash):
    if method_hash != baseline_hash:
        raise ValueError("Paired generated-data hashes differ")
    expected_baseline = {name + "_ate" for name in pilot.BASELINE_METHODS} | {"target_only_ate"}
    if set(baseline) != expected_baseline:
        raise ValueError("A requested method or baseline is missing")
    for field in ("estimate", "se", "truth"):
        if abs(float(method["target_only_ate"][field])-float(baseline["target_only_ate"][field])) > 1e-12:
            raise ValueError("Paired ordinary-target reference differs: " + field)
    combined = dict(method); combined.update(baseline)
    return primary_records(task, combined, method_hash)


def review(study, output, panel, dimension, sources=None, methods_only=False):
    study, output = map(pilot.scratch_path, (study, output))
    with (study / "panel_index.csv").open() as stream:
        pairs = index_pairs(list(csv.DictReader(stream)), panel, dimension, sources)
    roles = ("method",) if methods_only else ("method", "baseline")
    roots = sorted({pair[role]["root"] for pair in pairs.values() for role in roles})
    configurations, manifests = {}, {}
    for name in roots:
        root = pilot.scratch_path(name)
        configurations[name] = pilot.verify(root)
        _, tasks = read_campaign(root)
        manifests[name] = {str(row["task_id"]): row for row in tasks}

    def read_setting(pair):
        locations = [pair[role] for role in roles]
        tasks, configurations_for_pair, results, hashes = [], [], [], []
        for row in locations:
            configuration = configurations[row["root"]]
            task = manifests[row["root"]][str(row["task_id"])]
            for field in DATA_FIELDS:
                actual = configuration.get("n_deviated_sites", 1) if field == "n_deviated_sites" else task[field]
                expected = row[field]
                matches = float(actual) == float(expected) if field in ("K", "p", "rho", "sim_id", "n_deviated_sites") else actual == expected
                if not matches:
                    raise ValueError("Panel index and frozen task differ: " + field)
            if row["role"] == "method" and (configuration.get("crossfit_layers") != 2 or
                    configuration.get("source_validation_method") != "outer_fit" or
                    configuration.get("source_nuisance_method", "calibrated") != "calibrated" or
                    configuration.get("target_nuisance_method", "hou_calibrated") != "hou_calibrated"):
                raise ValueError("The primary panel requires the full two-layer method")
            result, data_hash = completed_result(Path(row["root"]), task)
            tasks.append(task); configurations_for_pair.append(configuration)
            results.append(result); hashes.append(data_hash)
        record = {field: locations[0][field] for field in DATA_FIELDS}
        record.update(protocol=tasks[0]["protocol"], panel=panel)
        if methods_only:
            return primary_records(record, results[0], hashes[0])
        if data_identity(tasks[0], configurations_for_pair[0]) != data_identity(tasks[1], configurations_for_pair[1]):
            raise ValueError("Paired scientific data specifications differ")
        for field in ("n_per_site", "n_folds", "nlambda", "nuisance_tol"):
            if configurations_for_pair[0][field] != configurations_for_pair[1][field]:
                raise ValueError("Paired fitting controls differ: " + field)
        return combine_pair(record, results[0], results[1], hashes[0], hashes[1])

    with ThreadPoolExecutor(max_workers=8) as pool:
        records = [row for result in pool.map(read_setting, pairs.values()) for row in result]
    metrics = summarize_paired_estimates(records)
    if any(row["repeats"] != 200 for row in metrics):
        raise ValueError("An output cell is not a complete200-repeat comparison")
    output.mkdir(parents=True, exist_ok=False)
    write_csv(output / "paired_estimates.csv", records)
    write_csv(output / "metrics.csv", metrics)
    pilot.write_json(output / "pairing_checks.json", dict(panel=panel, dimension=dimension,
        primary_repeats=len(pairs), methods_only=methods_only,
        source_counts=sorted({int(pair["method"]["K"]) for pair in pairs.values()}),
        paired_repeats=0 if methods_only else len(pairs),
        data_hash_matches=0 if methods_only else len(pairs), target_reference_matches=0 if methods_only else len(pairs),
        crossfit_layers=2, prescribed_repeats_per_cell=200, study=str(study),
        panel_index_sha256=pilot.digest(study / "panel_index.csv"),
        roots=roots, scientific_configuration_hashes={name: pilot.digest(Path(name) / "configuration.json") for name in roots}))
    lines = ["# Complete two-layer MC200 comparison", "",
        ("Panel: {}. p={}. {} completed primary repeats; only the ordinary/calibrated target references from the same worker are included. External baseline pairing is not asserted.".format(panel, dimension, len(pairs)) if methods_only else
         "Panel: {}. p={}. {} method/baseline pairs passed frozen-configuration, data-hash and ordinary-target checks.".format(panel, dimension, len(pairs))),
        "Every displayed cell contains all200 prescribed seeds. Comparisons retain the identical-data target reference.", "",
        "| Scenario | K | rho | Method | Bias | RMSE | Coverage | SE/SD |",
        "|---|---:|---:|---|---:|---:|---:|---:|"]
    for row in metrics:
        lines.append("|{config}|{K}|{rho}|{method}|{bias:.5f}|{rmse:.5f}|{coverage:.1%}|{se_to_sd_ratio:.3f}|".format(**row))
    lines += ["", "Wilson coverage intervals and paired MSE contrasts with Monte Carlo SEs are in metrics.csv.",
        "Do not pool correlated scenarios as independent repeats. Mean SE/empirical SD is descriptive at200 repeats.",
        "Sample-size/IVW estimands and each DR baseline's working-model assumptions must be distinguished from source calibration.", ""]
    (output / "report.md").write_text("\n".join(lines))
    print("Verified {} complete {} repeats; {} metric rows.".format(
        len(pairs), "primary" if methods_only else "paired", len(metrics)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--panel", choices=("main", "c4", "half_invalid"), default="main")
    parser.add_argument("--dimension", choices=(100, 200), type=int, required=True)
    parser.add_argument("--sources", choices=(2,4,6,8), type=int, nargs="+")
    parser.add_argument("--methods-only", action="store_true", help="Review complete primary repeats without asserting completed external baselines")
    args = parser.parse_args()
    review(args.study, args.output, args.panel, args.dimension, args.sources, args.methods_only)
