#!/usr/bin/env python3
"""Prepare the fixed-K, two-layer MC200 panel and reuse matching baselines.

No submissions are performed here. Existing checked launchers and controllers
execute the frozen campaigns; every new repeat remains a separate Slurm job.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
from prepare_highdim_validation import (EXCLUDED_NODES, data_identity,
                                       experiment_identity, read_campaign)

DIMENSIONS = (100, 200)
SOURCE_COUNTS = (2, 4, 6, 8)
SCENARIOS = ("C1", "C2", "C3")
DEVIATIONS = (0, .5, 1, 1.5, 2)
SEEDS = tuple(range(1, 201))


def scientific_configuration(role="method", n_deviated_sites=1):
    configuration = dict(n_per_site=1000, dgp_type="bounded",
        dgp_control=dict(covariate_shift=1, source_treatment_scale=1),
        n_folds=10, nlambda=100, nuisance_tol=1e-10, source_radius=12,
        target_score_radius=12, aggregation_lambda=.5, n_deviated_sites=n_deviated_sites)
    if role == "baseline":
        configuration.update(baseline_methods=list(pilot.BASELINE_METHODS),
                             baseline_variance_method="bootstrap", n_bootstrap=5000)
    else:
        configuration.update(crossfit_layers=2, source_validation_method="outer_fit",
            source_nuisance_method="standard" if role == "standard_source" else "calibrated",
            target_nuisance_method="lasso" if role == "ordinary_target" else "hou_calibrated",
            recipes=["score_derivative"])
    return configuration


def existing_roots(implementation):
    base = implementation / "r7"
    roots = [base / "validation_highdim_mc200_cutoff2_v1", base / "baseline_mc200_v1"]
    for name in ("outcome_shift", "outcome_shift_baselines"):
        roots.append(base / "extended_validation_v1" / name)
    for name in ("highdim_validation_v2", "source_count_validation_v1", "source_count_weak_fraction_v1"):
        registry = json.loads((base / name / "campaign_registry.json").read_text())
        roots.extend(Path(entry["root"]) for entry in registry if entry["origin"] == "new")
    # These seeds are401–450. Keep them in the inventory, but never relabel
    # them as the1–200 primary panel or confuse three layers with two.
    roots.extend(implementation / "r9/layer_validation_v2" / ("layers2_p" + str(p))
                 for p in DIMENSIONS)
    return list(dict.fromkeys(roots))


def catalogue(roots):
    records = {}
    metadata = []
    for root in roots:
        configuration = pilot.verify(root)
        _, tasks = read_campaign(root)
        metadata.append(dict(root=str(root), planned=len(tasks),
            manifest_sha256=pilot.digest(root / "manifest.csv"),
            configuration_sha256=pilot.digest(root / "configuration.json")))
        for task in tasks:
            key = experiment_identity(task, configuration)
            directory = root / "tasks" / str(task["task_id"])
            record = dict(root=str(root), task_id=int(task["task_id"]),
                          complete=(directory / "COMPLETE").is_file())
            if key in records:
                raise ValueError("Duplicate canonical experiment: " + str((record, records[key])))
            records[key] = record
    return records, metadata


def rectangular_groups(tasks):
    """Group missing cells only when their exact prescribed seed sets agree."""
    cells = defaultdict(list)
    for task in tasks:
        cells[(task["panel"], int(task["p"]), int(task["K"]), task["n_deviated_sites"],
               task["config"], float(task["rho"]))].append(int(task["sim_id"]))
    grouped = defaultdict(list)
    for (panel, dimension, sources, deviated, scenario, rho), seeds in cells.items():
        if len(seeds) != len(set(seeds)):
            raise ValueError("Duplicate prescribed seed")
        grouped[(panel, dimension, sources, deviated, tuple(sorted(seeds)))].append((scenario, rho))
    result = []
    for (panel, dimension, sources, deviated, seeds), pairs in sorted(grouped.items()):
        configs = sorted({pair[0] for pair in pairs})
        rhos = sorted({pair[1] for pair in pairs})
        if set(pairs) != {(config, rho) for config in configs for rho in rhos}:
            raise ValueError("Missing cells are not a rectangular campaign; split them explicitly")
        result.append(dict(panel=panel, dimension=dimension, sources=sources, n_deviated_sites=deviated,
                           configs=configs, rhos=rhos, seeds=list(seeds)))
    return result


def resolve_requests(existing):
    requests, reused, missing = [], [], defaultdict(list)
    definitions = []
    for role in ("method", "baseline"):
        definitions.extend([(role, "main", list(SCENARIOS), list(SOURCE_COUNTS), list(DEVIATIONS)),
            (role, "c4", ["C4"], list(SOURCE_COUNTS), list(DEVIATIONS)),
            (role, "half_invalid", list(SCENARIOS), [4, 6, 8], [.5, 1, 1.5, 2])])
    definitions += [(role, "calibration", list(SCENARIOS), [2], [0])
                    for role in ("standard_source", "ordinary_target")]
    for role, panel, configs, sources, rhos in definitions:
        tasks = pilot.task_rows(configs, sources, rhos, SEEDS,
                               "one_round", "treated_arm", list(DIMENSIONS))
        for task in tasks:
            task.update(panel=panel, n_deviated_sites=int(task["K"])//2 if panel == "half_invalid" else 1)
            configuration = scientific_configuration(role, task["n_deviated_sites"])
            identity = experiment_identity(task, configuration)
            record = dict(role=role, task=task)
            if identity in existing:
                previous = existing[identity]
                if not previous["complete"]:
                    raise ValueError("Existing unfinished experiment needs recovery, not duplicate submission: " + str(previous))
                record.update(previous)
                reused.append(record)
            else:
                missing[role].append(task)
            requests.append(record)
    return requests, reused, missing


def prepare(implementation, output, check_root):
    implementation, output, check_root = map(pilot.scratch_path, (implementation, output, check_root))
    for marker in ("TESTS_PASSED", "CROSSFIT_LAYER_CHECKS_PASSED", "CV_CERTIFICATE_CHECKS_PASSED",
                   "SOURCE_STANDARD_CHECKS_PASSED"):
        if not (check_root / marker).is_file():
            raise ValueError("Missing checked-library capability: " + marker)
    existing, metadata = catalogue(existing_roots(implementation))
    requests, reused, missing = resolve_requests(existing)
    counts = Counter(record["role"] for record in requests)
    if counts != Counter(method=46400, baseline=46400, standard_source=1200, ordinary_target=1200):
        raise ValueError("The prescribed MC200 grid is incomplete")
    output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    workflow = output / "coordinator_workflow"
    workflow.mkdir()
    for name in (Path(__file__).name, "submit_repeat_pilot.py", "advance_validation_plan.py",
                 "advance_validation_plan.sh", "prepare_highdim_validation.py"):
        shutil.copy2(source / name, workflow / name)
    pilot.write_json(workflow / "manifest.json", {str(path): pilot.digest(path) for path in workflow.iterdir()})
    pilot.write_json(output / "existing_campaigns.json", metadata)
    pilot.write_json(output / "reused_results.json", reused)
    campaigns, assignments = [], {}
    # First-seed execution checks precede every remaining seed. Supplementary
    # campaigns also wait for the corresponding main pipeline's first wave.
    caps = dict(method=200, baseline=50, standard_source=16, ordinary_target=16)
    for role in ("method", "baseline", "standard_source", "ordinary_target"):
        groups = rectangular_groups(missing[role])
        groups.sort(key=lambda item: (item["panel"] != "main", item["panel"], item["dimension"], item["sources"]))
        for group in groups:
            p, sources = group["dimension"], group["sources"]
            name = "{}_{}_p{}_k{}_s{}_{}".format(group["panel"], role, p, sources,
                                                min(group["seeds"]), max(group["seeds"]))
            root = output / name
            command = [sys.executable, "-B", str(source / "submit_repeat_pilot.py"), "prepare", str(root),
                "--check-root", str(check_root), "--dimensions", str(p), "--sources", str(sources),
                "--configs"] + group["configs"] + ["--rhos"] + [str(rho) for rho in group["rhos"]]
            command += ["--repeats"] + [str(seed) for seed in group["seeds"]]
            command += ["--aggregation-cutoff", "2", "--cpus", "2", "--checkpoint",
                "--n-deviated-sites", str(group["n_deviated_sites"]),
                "--nuisance-cv-certificate", "--exclude-nodes", EXCLUDED_NODES,
                "--partitions", pilot.DEFAULT_PARTITIONS, "--save-fits", "first",
                "--memory", "8G" if role == "baseline" else "16G",
                "--time-limit", "02:00:00" if role == "baseline" else "04:00:00"]
            if role == "baseline":
                command += ["--baseline-methods"] + list(pilot.BASELINE_METHODS)
                references = [campaign["root"] for campaign in campaigns if campaign["role"] == "method"
                              and campaign["dimension"] == p and campaign["sources"] == sources
                              and campaign["panel"] == group["panel"]]
                command += ["--baseline-reference-roots"] + references
            else:
                command += ["--crossfit-layers", "2", "--recipes", "score_derivative"]
                if role == "standard_source":
                    command += ["--source-nuisance-method", "standard"]
                if role == "ordinary_target":
                    command += ["--target-nuisance-method", "lasso"]
            subprocess.run(command, check=True, env=pilot.submission_environment())
            configuration, tasks = read_campaign(root)
            for task in tasks:
                identity = experiment_identity(task, configuration)
                if identity in existing or identity in assignments:
                    raise ValueError("Attempted duplicate scientific experiment")
                assignments[identity] = dict(root=str(root), task_id=int(task["task_id"]))
            gates = []
            if group["panel"] != "main":
                parent = next(item for item in campaigns if item["role"] == "method"
                    and item["panel"] == "main" and item["dimension"] == p and item["sources"] == sources)
                gates.append(dict(root=parent["root"], tasks=list(range(1, len(SCENARIOS)*len(DEVIATIONS)+1))))
            maximum = (8 if role == "baseline" else 16) if group["panel"] in ("c4", "half_invalid") else caps[role]
            campaigns.append(dict(name=name, root=str(root), role=role, panel=group["panel"],
                dimension=p, sources=sources, n_deviated_sites=group["n_deviated_sites"],
                planned=len(tasks), max_in_flight=maximum, gates=gates))
    reused_lookup = {(record["role"], data_identity(record["task"],
                     scientific_configuration(record["role"], record["task"]["n_deviated_sites"]))): record
                     for record in reused}
    panel = []
    for request in requests:
        role, task = request["role"], request["task"]
        configuration = scientific_configuration(role, task["n_deviated_sites"])
        old = reused_lookup.get((role, data_identity(task, configuration)))
        location = old or assignments[experiment_identity(task, configuration)]
        record = {key: task[key] for key in ("config", "K", "p", "rho", "sim_id", "protocol", "deviation_mechanism")}
        record.update(role=role, panel=task["panel"], n_deviated_sites=task["n_deviated_sites"],
                      root=location["root"], task_id=location["task_id"], reused=old is not None)
        panel.append(record)
    with (output / "panel_index.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(panel[0]))
        writer.writeheader()
        writer.writerows(panel)
    paired = defaultdict(set)
    for record in panel:
        key = tuple(record[field] for field in ("config", "K", "p", "rho", "sim_id", "n_deviated_sites"))
        if record["role"] in paired[key]:
            raise ValueError("Duplicated panel role")
        paired[key].add(record["role"])
    if len(paired) != 46400 or any(not {"method", "baseline"}.issubset(roles) for roles in paired.values()):
        raise ValueError("A prescribed method/baseline pair is missing")
    plan = dict(purpose="Current-paper fixed-K two-layer MC200; weak departures are reported boundaries",
        scientific_scope=dict(dimensions=list(DIMENSIONS), source_counts=list(SOURCE_COUNTS),
            configs=list(SCENARIOS), rhos=list(DEVIATIONS), repeats=200, n_per_site=1000,
            crossfit_layers=2, aggregation_cutoff=2), campaigns=campaigns,
        supplementary="C4 and half-invalid panels have lower concurrency and execution gates; K2 half-invalid and rho0 controls reuse the main panel",
        deferred="Growing-K and weak-bias solutions belong to the next paper")
    pilot.write_json(output / "expansion_plan.json", plan)
    pilot.write_json(output / "design_checks.json", dict(requested=dict(counts),
        reused=dict(Counter(row["role"] for row in reused)),
        new=dict(Counter({role: len(tasks) for role, tasks in missing.items()})),
        paired_data_settings=len(paired), main_data_settings=24000, source_calibration_pairs=1200,
        maximum_in_flight=sum(campaign["max_in_flight"] for campaign in campaigns)))
    print(json.dumps(json.loads((output / "design_checks.json").read_text()), sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("implementation", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-root", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.implementation, args.output, args.check_root)


if __name__ == "__main__":
    main()
