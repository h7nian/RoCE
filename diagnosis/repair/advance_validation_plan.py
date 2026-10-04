#!/usr/bin/env python3
"""Release prepared validation campaigns only after their execution checks pass."""
import argparse
import csv
import json
import math
import re
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot


def validate_layer_pair_configurations(first, second):
    if sorted([first.get("crossfit_layers", 0), second.get("crossfit_layers", 0)]) != [2, 3]:
        raise ValueError("A layer-pair gate requires one two-layer and one three-layer campaign")
    exclusions = {"scope", "check_root", "library", "cpus", "memory", "time_limit", "partitions",
                  "exclude_nodes", "account", "shared_cache", "checkpoint", "max_restarts", "save_fits",
                  "crossfit_layers", "source_validation_method"}
    scientific = lambda config: {key: value for key, value in config.items() if key not in exclusions}
    if scientific(first) != scientific(second):
        raise ValueError("Layer-pair scientific configurations differ beyond the layer choice")


def check_gate(gate):
    root = pilot.scratch_path(gate["root"])
    configuration = pilot.verify(root)
    paired_root = pilot.scratch_path(gate["paired_root"]) if gate.get("paired_root") else None
    if paired_root is not None:
        paired_configuration = pilot.verify(paired_root)
        validate_layer_pair_configurations(configuration, paired_configuration)
        with (root / "manifest.csv").open() as stream:
            first_manifest = {row["task_id"]: row for row in csv.DictReader(stream)}
        with (paired_root / "manifest.csv").open() as stream:
            second_manifest = {row["task_id"]: row for row in csv.DictReader(stream)}
    if gate.get("require_active_controller"):
        status = root / "supervisor_status.json"
        if not status.is_file():
            return False, "controller has not reported its state"
        try:
            state = json.loads(status.read_text())["state"]
        except (ValueError, KeyError):
            return False, "controller status update is not yet readable"
        if state not in {"RUNNING", "COMPLETE", "CHECKPOINTED_FOR_REQUEUE"}:
            return False, "controller state is " + state
    for task in gate["tasks"]:
        directory = root / "tasks" / str(task)
        if not (directory / "COMPLETE").is_file():
            return False, "task {} has not committed results".format(task)
        if not (directory / "data_sha256.txt").is_file():
            raise ValueError("Completed task is missing its data hash: " + str(directory))
        if ((directory / "status.txt").read_text().strip() != "COMPLETE" or
                not re.fullmatch(r"[0-9a-f]{64}", (directory / "data_sha256.txt").read_text().strip())):
            raise ValueError("Inconsistent completion metadata: " + str(directory))
        expected_hash = gate.get("expected_data_hashes", {}).get(str(task))
        if expected_hash is not None and (directory / "data_sha256.txt").read_text().strip() != expected_hash:
            raise ValueError("Generated data do not match the expected design: " + str(directory))
        if paired_root is not None:
            partner = paired_root / "tasks" / str(task)
            if not (partner / "COMPLETE").is_file():
                return False, "paired task {} has not committed results".format(task)
            if first_manifest.get(str(task)) != second_manifest.get(str(task)):
                raise ValueError("Layer-pair manifest tasks differ")
            for filename in ["data_sha256.txt"] + [recipe + "_outer_fit_sha256.txt" for recipe in configuration["recipes"]]:
                if not (directory / filename).is_file() or not (partner / filename).is_file():
                    raise ValueError("Layer-pair fingerprint is missing: " + filename)
                first = (directory / filename).read_text().strip()
                second = (partner / filename).read_text().strip()
                if not re.fullmatch(r"[0-9a-f]{64}", first) or first != second:
                    raise ValueError("Layer-pair data or outer nuisance fits differ: " + filename)
        with (directory / "results.csv").open() as stream:
            rows = [row for row in csv.DictReader(stream) if row["estimand_scope"] == "tate"]
        if not rows or any(not math.isfinite(float(row[key])) for row in rows
                           for key in ("estimate", "se", "truth")) or any(float(row["se"]) <= 0 for row in rows):
            raise ValueError("Invalid completed validation result: " + str(directory))
        if len({(row["pilot_recipe"], row["method"]) for row in rows}) != len(rows):
            raise ValueError("Duplicated method results: " + str(directory))
        recipes = configuration["recipes"]
        if configuration.get("layer_diagnostics"):
            with (directory / "layer_estimates.csv").open() as stream:
                layer_rows = list(csv.DictReader(stream))
            expected_rows = {(recipe, mode, quantity) for recipe in recipes
                             for mode in ("common_tate", "separate_arms", "joint_tate")
                             for quantity in ("mu1", "mu0", "tate")}
            actual_rows = {(row["recipe"], row["aggregation_mode"], row["quantity"]) for row in layer_rows}
            if len(layer_rows) != len(expected_rows) or actual_rows != expected_rows:
                raise ValueError("Incomplete or duplicated layer variance diagnostics")
            if any(int(row["crossfit_layers"]) != configuration["crossfit_layers"] or
                   any(not math.isfinite(float(row[key])) for key in
                       ("estimate", "truth", "variance", "variance_fixed_weights", "covariance_mu1_mu0",
                        "covariance_mu1_mu0_fixed_weights")) or
                   min(float(row["variance"]), float(row["variance_fixed_weights"])) <= 0 for row in layer_rows):
                raise ValueError("Invalid layer variance diagnostics")
        for recipe in recipes:
            methods = {row["method"] for row in rows if row["pilot_recipe"] == recipe}
            if configuration.get("baseline_methods"):
                expected = {method + "_ate" for method in configuration["baseline_methods"]} | {"target_only_ate"}
            else:
                with (root / "manifest.csv").open() as stream:
                    task_row = next(row for row in csv.DictReader(stream) if row["task_id"] == str(task))
                prefix = task_row["protocol"] + "_crossfit_ate"
                expected = {prefix, prefix + "_separate_arms", prefix + "_joint_tate", "target_only_ate"}
            if not expected.issubset(methods):
                raise ValueError("Missing requested methods: " + str(directory))
    if paired_root is not None:
        ready, reason = check_gate(dict(root=str(paired_root), tasks=gate["tasks"]))
        if not ready:
            return False, "paired campaign: " + reason
    # Coverage and RMSE are outcomes to report, never filters for accepting seeds.
    return True, "all designated execution checks passed"


def advance(plan, status_path):
    states = []
    for campaign in plan["campaigns"]:
        root = pilot.scratch_path(campaign["root"])
        receipt = root / "controller_submission.json"
        if receipt.exists():
            submission = json.loads(receipt.read_text())
            if submission.get("state") != "submitted":
                raise ValueError("Review the existing ambiguous controller receipt: " + str(receipt))
            states.append(dict(name=campaign["name"], state="SUBMITTED",
                               controller_job=submission["stdout"].strip()))
            continue
        failures = []
        for gate in campaign["gates"]:
            ready, reason = check_gate(gate)
            if not ready:
                failures.append(reason)
        if failures:
            states.append(dict(name=campaign["name"], state="WAITING", reasons=failures))
            continue
        workflow = pilot.scratch_path(campaign.get("controller_workflow", root / "workflow"))
        command = [sys.executable, "-B", str(workflow / "submit_repeat_pilot.py"),
                   "supervise", str(root), "--max-in-flight", str(campaign["max_in_flight"])]
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, timeout=120)
        result.check_returncode()
        submission = json.loads(receipt.read_text())
        states.append(dict(name=campaign["name"], state="SUBMITTED",
                           controller_job=submission["stdout"].strip()))
    complete = all(state["state"] == "SUBMITTED" for state in states)
    pilot.replace_json(status_path, dict(state="ALL_CONTROLLERS_SUBMITTED" if complete else "WAITING_FOR_CHECKS",
                                        campaigns=states))
    return complete


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("--hours", type=float, default=8)
    args = parser.parse_args()
    path = pilot.scratch_path(args.plan)
    if not 0 < args.hours <= 24:
        parser.error("hours must be in (0, 24]")
    plan = json.loads(path.read_text())
    status_path = path.with_name("expansion_status.json")
    deadline = time.monotonic() + args.hours * 3600
    try:
        while time.monotonic() < deadline:
            if advance(plan, status_path):
                return 0
            time.sleep(30)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        pilot.replace_json(status_path, dict(state="FAILED", error=str(error)))
        raise
    return 75


if __name__ == "__main__":
    raise SystemExit(main())
