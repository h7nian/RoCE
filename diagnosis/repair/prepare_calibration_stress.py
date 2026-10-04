#!/usr/bin/env python3
"""Prepare a separate covariate-shift stress study using the checked workflows."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import prepare_current_paper_mc200 as current
from prepare_highdim_validation import EXCLUDED_NODES, experiment_identity, read_campaign


def configuration(role):
    result = current.scientific_configuration(role)
    result["dgp_control"]["covariate_shift"] = 2
    return result


def prepare(implementation, main_root, output, check_root):
    implementation, main_root, output, check_root = map(pilot.scratch_path,
        (implementation, main_root, output, check_root))
    roots = current.existing_roots(implementation) + [implementation / "r7/extended_validation_v1" / name
        for name in ("covariate_shift_double", "covariate_shift_double_baselines")]
    existing, metadata = current.catalogue(roots)
    tasks = pilot.task_rows(list(current.SCENARIOS), [2], [0], current.SEEDS,
                           "one_round", "treated_arm", list(current.DIMENSIONS))
    missing, reused = {}, []
    for role in ("method", "baseline", "standard_source"):
        missing[role] = []
        for task in tasks:
            key = experiment_identity(task, configuration(role))
            if key in existing:
                if not existing[key]["complete"]:
                    raise ValueError("Existing unfinished stress experiment requires recovery")
                reused.append(dict(role=role, task=task, **existing[key]))
            else:
                missing[role].append(dict(task, panel="calibration_stress", n_deviated_sites=1))
    output.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    workflow = output / "coordinator_workflow"
    workflow.mkdir()
    for name in (Path(__file__).name, "prepare_current_paper_mc200.py", "prepare_highdim_validation.py",
                 "submit_repeat_pilot.py", "advance_validation_plan.py", "advance_validation_plan.sh"):
        shutil.copy2(source / name, workflow / name)
    pilot.write_json(workflow / "manifest.json", {str(path): pilot.digest(path) for path in workflow.iterdir()})
    pilot.write_json(output / "existing_campaigns.json", metadata)
    pilot.write_json(output / "reused_results.json", reused)
    campaigns, assignments = [], {}
    for role in ("method", "baseline", "standard_source"):
        for group in current.rectangular_groups(missing[role]):
            dimension = group["dimension"]
            root = output / "{}_p{}_k2_s{}_{}".format(role, dimension, min(group["seeds"]), max(group["seeds"]))
            command = [sys.executable, "-B", str(source / "submit_repeat_pilot.py"), "prepare", str(root),
                "--check-root", str(check_root), "--dimensions", str(dimension), "--sources", "2",
                "--configs"] + group["configs"] + ["--rhos", "0", "--repeats"] + [str(seed) for seed in group["seeds"]]
            command += ["--covariate-shift", "2", "--aggregation-cutoff", "2", "--cpus", "2",
                "--checkpoint", "--nuisance-cv-certificate", "--partitions", pilot.DEFAULT_PARTITIONS,
                "--exclude-nodes", EXCLUDED_NODES, "--save-fits", "first",
                "--memory", "8G" if role == "baseline" else "16G", "--time-limit", "04:00:00"]
            if role == "baseline":
                reference = next(row["root"] for row in campaigns if row["role"] == "method" and row["dimension"] == dimension)
                command += ["--baseline-methods"] + list(pilot.BASELINE_METHODS) + ["--baseline-reference-roots", reference]
            else:
                command += ["--crossfit-layers", "2", "--recipes", "score_derivative"]
                if role == "standard_source":
                    command += ["--source-nuisance-method", "standard"]
            subprocess.run(command, check=True, env=pilot.submission_environment())
            actual, rows = read_campaign(root)
            for row in rows:
                key = experiment_identity(row, actual)
                if key != experiment_identity(row, configuration(role)) or key in assignments or key in existing:
                    raise ValueError("Stress configuration mismatch or duplicate experiment")
                assignments[key] = dict(root=str(root), task_id=int(row["task_id"]))
            parent = main_root / "main_method_p{}_k2_s1_200".format(dimension)
            gates = [dict(root=str(parent), tasks=list(range(1,16)))]
            if role == "standard_source":
                reference = next(row["root"] for row in campaigns if row["role"] == "method" and row["dimension"] == dimension)
                gates.append(dict(root=reference, tasks=[1,2,3]))
            campaigns.append(dict(name=root.name, root=str(root), role=role, dimension=dimension,
                planned=len(rows), max_in_flight=8 if role == "baseline" else 16, gates=gates))
    index = []
    for role in ("method", "baseline", "standard_source"):
        for task in tasks:
            key = experiment_identity(task, configuration(role))
            location = existing.get(key, assignments.get(key))
            if location is None:
                raise ValueError("Unassigned stress repeat")
            index.append(dict(role=role, **{name: task[name] for name in
                ("config","K","p","rho","sim_id","protocol","deviation_mechanism")},
                root=location["root"], task_id=location["task_id"], reused=key in existing))
    with (output / "panel_index.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(index[0]))
        writer.writeheader(); writer.writerows(index)
    pilot.write_json(output / "expansion_plan.json", dict(
        purpose="Prespecified source-calibration stress comparison: covariate_shift=2, outcome rho=0",
        campaigns=campaigns, scope="Two layers, K2,1000/site,p100/200,C1-C3,200 paired repeats; old shift1 results remain intact"))
    checks = dict(requested=3600, reused=dict(Counter(row["role"] for row in reused)),
        new={role: len(rows) for role, rows in missing.items()},
        maximum_in_flight=sum(row["max_in_flight"] for row in campaigns))
    pilot.write_json(output / "design_checks.json", checks)
    print(json.dumps(checks, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("implementation", type=Path)
    parser.add_argument("main_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--check-root", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.implementation, args.main_root, args.output, args.check_root)
