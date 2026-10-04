#!/usr/bin/env python3
"""Prepare and submit a development pilot: one scenario/repeat per Slurm job.

All calibration recipes and aggregation modes within a task use the same seed.
This launcher does not authorize a formal coverage study or select a method.
"""

import argparse
import csv
import fcntl
import hashlib
import itertools
import os
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from uuid import uuid4


SCRATCH = Path("/scratch.global/zhan9381/FACE-HD")
DEFAULT_PARTITIONS = "preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb"
DEFAULT_CONTROLLER_PARTITIONS = DEFAULT_PARTITIONS
DEFAULT_MAX_IN_FLIGHT = 512
MAX_IN_FLIGHT = 2000
MAX_SUBMISSION_BATCH = 256
BASELINE_METHODS = ("sample_size", "inverse_variance", "federated_dr", "pooled_dr")


def scratch_path(value):
    path = Path(value).resolve()
    if SCRATCH not in path.parents:
        raise ValueError("Output must be below the FACE-HD scratch directory")
    return path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def replace_json(path, value):
    temporary = path.with_name(path.name + "." + uuid4().hex + ".tmp")
    try:
        write_json(temporary, value)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def submission_comment(root, task):
    campaign = hashlib.sha256(str(root).encode()).hexdigest()[:20]
    return "roce:" + campaign + ":" + str(task)


def submission_environment():
    """Let each new job obtain its own restart counter from Slurm.

    Slurm only supplies this variable for restarted jobs. Exporting a
    controller's counter can therefore mislabel a fresh worker's first attempt.
    Preserve the user's module/library environment.
    """
    environment = os.environ.copy()
    environment.pop("SLURM_RESTART_COUNT", None)
    return environment


def task_rows(configs, sources, rhos, repeats, protocol, mechanism, dimensions=None):
    if mechanism == "both_arms" and configs != ["C1"]:
        raise ValueError("The prescribed shared-shift design contains C1 only")
    rows = []
    # Interleave cells within each seed so the first wave checks every dimension.
    for repeat, dimension, config, count, rho in itertools.product(
            repeats, dimensions or [100], configs, sources, rhos):
        rows.append(dict(task_id=len(rows) + 1, config=config, K=count, p=dimension, rho=rho,
                         sim_id=repeat, protocol=protocol, deviation_mechanism=mechanism))
    identities = [tuple(row.values())[1:] for row in rows]
    if len(identities) != len(set(identities)):
        raise ValueError("Duplicate scenario/repeat keys are not allowed")
    return rows


def prepare(args):
    root = scratch_path(args.output)
    check = scratch_path(args.check_root)
    if not (check / "TESTS_PASSED").is_file():
        raise ValueError("A completed regression check is required")
    library = check / "Rlib"
    if not (library / "RoCE" / "DESCRIPTION").is_file():
        raise ValueError("The checked installed library is missing")
    dgp_type = getattr(args, "dgp_type", "bounded")
    dimensions = getattr(args, "dimensions", None) or ([10] if dgp_type == "bounded" else [100])
    repeats = args.repeats
    if repeats is None:
        repeats = list(range(args.seed_start, args.seed_start + (args.n_repeats or 1)))
    if any(seed < 1 for seed in repeats):
        raise ValueError("Repeat seeds must be positive")
    baseline_methods = getattr(args, "baseline_methods", None)
    certificate = getattr(args, "nuisance_cv_certificate", False)
    if not isinstance(certificate, bool):
        raise ValueError("nuisance_cv_certificate must be boolean")
    if certificate and not (check / "CV_CERTIFICATE_CHECKS_PASSED").is_file():
        raise ValueError("The checked library must pass the full CV certificate checks")
    source_program = getattr(args, "source_nuisance_method", "calibrated")
    if source_program not in ("calibrated", "standard"):
        raise ValueError("source_nuisance_method must be calibrated or standard")
    layers = getattr(args, "crossfit_layers", None)
    if layers is None and getattr(args, "source_validation_method", "calibrated") == "outer_fit":
        layers = 2
    if layers is not None:
        if layers not in (2, 3) or isinstance(layers, bool):
            raise ValueError("crossfit_layers must be 2 or 3")
        layers = int(layers)
        if baseline_methods:
            raise ValueError("Cross-fitting layer options do not apply to the baseline worker")
        if source_program == "standard" and getattr(args, "target_nuisance_method", "hou_calibrated") == "lasso":
            raise ValueError("Layer selection requires a calibrated source or target")
        if getattr(args, "source_validation_method", "calibrated") == "initial":
            raise ValueError("Layer selection cannot use legacy initial validation")
        if layers == 3 and getattr(args, "source_validation_method", "calibrated") == "outer_fit":
            raise ValueError("outer_fit uses two layers")
        if not (check / "CROSSFIT_LAYER_CHECKS_PASSED").is_file():
            raise ValueError("The checked library must pass the two/three-layer regression checks")
        args.source_validation_method = "outer_fit" if layers == 2 else "calibrated"
    n_deviated_sites = getattr(args, "n_deviated_sites", 1)
    if (not isinstance(n_deviated_sites, int) or isinstance(n_deviated_sites, bool) or
            n_deviated_sites < 0 or n_deviated_sites > min(args.sources)):
        raise ValueError("n_deviated_sites must be an integer between zero and every requested source count")
    if baseline_methods and args.recipes:
        raise ValueError("--recipes cannot be combined with --baseline-methods")
    if getattr(args, "baseline_reference_roots", None) and not baseline_methods:
        raise ValueError("--baseline-reference-roots requires --baseline-methods")
    if baseline_methods and (
            getattr(args, "target_nuisance_method", "hou_calibrated") != "hou_calibrated" or
            getattr(args, "source_validation_method", "calibrated") != "calibrated" or
            source_program != "calibrated"):
        raise ValueError("RoCE nuisance ablation options do not apply to the baseline worker")
    if baseline_methods and (len(set(baseline_methods)) != len(baseline_methods) or
                             not set(baseline_methods).issubset(BASELINE_METHODS)):
        raise ValueError("baseline_methods must contain distinct supported methods")
    if source_program == "standard":
        if args.protocol != "one_round" or getattr(args, "source_validation_method", "calibrated") == "initial":
            raise ValueError("Standard source fitting requires one_round and complete inner validation")
        if not (check / "SOURCE_STANDARD_CHECKS_PASSED").is_file():
            raise ValueError("The checked library must pass the standard-source regression checks")
        if getattr(args, "source_validation_method", "calibrated") != "outer_fit":
            args.source_validation_method = "complete"
    recipes = (["baselines"] if baseline_methods else args.recipes or
               (["score_derivative"] if dgp_type == "bounded" else ["legacy", "score_derivative"]))
    rows = task_rows(args.configs, args.sources, args.rhos, repeats,
                     args.protocol, args.deviation_mechanism, dimensions)
    if len(recipes) != len(set(recipes)):
        raise ValueError("Each calibration recipe must occur only once")
    checkpoint = getattr(args, "checkpoint", True)
    if checkpoint and args.shared_cache:
        raise ValueError("--checkpoint already provides persistent shared caching; omit --shared-cache")
    if "preempt" in args.partitions.split(",") and not checkpoint:
        raise ValueError("preempt requires --checkpoint (which also enables requeue)")
    aggregation_cutoff = getattr(args, "aggregation_cutoff", 1.0)
    if not math.isfinite(aggregation_cutoff) or aggregation_cutoff <= 0:
        raise ValueError("aggregation_cutoff must be finite and positive")
    excluded_nodes = getattr(args, "exclude_nodes", "")
    if not isinstance(excluded_nodes, str) or (excluded_nodes and not re.fullmatch(
            r"[A-Za-z0-9_-]+(?:,[A-Za-z0-9_-]+)*", excluded_nodes)):
        raise ValueError("exclude_nodes must be a comma-separated list of node names")
    root.mkdir(parents=True, exist_ok=False)
    (root / "workflow").mkdir()
    (root / "logs").mkdir()
    (root / "submissions").mkdir()
    (root / "tasks").mkdir()
    for name in (Path(__file__).name, "run_repeat_pilot.R", "run_repeat_pilot.sh",
                 "supervise_repeat_pilot.py", "supervise_repeat_pilot.sh", "summarize_repeat_pilot.py"):
        shutil.copy2(Path(__file__).parent / name, root / "workflow" / name)
    if baseline_methods:
        shutil.copy2(Path(__file__).with_name("run_baseline_repeat.R"), root / "workflow" / "run_baseline_repeat.R")
    if layers is not None:
        shutil.copy2(Path(__file__).with_name("layer_comparison_summary.R"), root / "workflow" / "layer_comparison_summary.R")
    with (root / "manifest.csv").open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    configuration = dict(
        scope="development pilot, not formal MC500", check_root=str(check),
        library=str(library), recipes=recipes, n_per_site=1000,
        p=dimensions[0] if len(dimensions) == 1 else None, dimensions=dimensions,
        dgp_type=dgp_type, dgp_control=dict(
            covariate_shift=getattr(args, "covariate_shift", 1),
            source_treatment_scale=getattr(args, "source_treatment_scale", 1)) if dgp_type == "bounded" else None,
        save_fits=getattr(args, "save_fits", None) or ("first" if dgp_type == "bounded" else "all"),
        first_seed=min(repeats),
        n_folds=10, nlambda=100, nuisance_tol=1e-10, source_radius=12 if dgp_type == "bounded" else 5,
        target_score_radius=12 if dgp_type == "bounded" else 2.1972245773362196,
        aggregation_lambda=1 / aggregation_cutoff,
        checkpoint=checkpoint, max_restarts=getattr(args, "max_restarts", 20),
        cpus=args.cpus, memory=args.memory, time_limit=args.time_limit,
        partitions=args.partitions, exclude_nodes=excluded_nodes,
        account=args.account, shared_cache=args.shared_cache)
    if baseline_methods:
        configuration.update(baseline_methods=baseline_methods,
            baseline_reference_roots=[str(scratch_path(path)) for path in
                                      (getattr(args, "baseline_reference_roots", None) or [])],
            baseline_variance_method="bootstrap", n_bootstrap=5000,
            worker_script="run_baseline_repeat.R")
    if layers is not None:
        configuration.update(crossfit_layers=layers, layer_diagnostics=True)
    if certificate:
        configuration["nuisance_cv_certificate"] = True
    for name, default in (("target_nuisance_method", "hou_calibrated"),
                          ("source_nuisance_method", "calibrated"),
                          ("source_validation_method", "calibrated")):
        value = getattr(args, name, default)
        if value != default:
            configuration[name] = value
    if n_deviated_sites != 1:
        configuration["n_deviated_sites"] = n_deviated_sites
    write_json(root / "configuration.json", configuration)
    imported_paths = import_completed(root, getattr(args, "completed_from", None), rows)
    paths = [root / "configuration.json", root / "manifest.csv",
             check / "source_manifest.json", check / "TESTS_PASSED"]
    if source_program == "standard":
        paths.append(check / "SOURCE_STANDARD_CHECKS_PASSED")
    if layers is not None:
        paths.append(check / "CROSSFIT_LAYER_CHECKS_PASSED")
    if certificate:
        paths.append(check / "CV_CERTIFICATE_CHECKS_PASSED")
    paths += imported_paths
    paths += sorted((root / "workflow").iterdir())
    paths += sorted(path for path in (library / "RoCE").rglob("*") if path.is_file())
    write_json(root / "provenance.json", {str(path): digest(path) for path in paths})
    print(f"Prepared {len(rows)} scenario/repeat jobs under {root}")


def import_completed(root, previous_root, rows):
    if previous_root is None:
        return []
    previous_root = scratch_path(previous_root)
    previous = verify(previous_root)
    current = json.loads((root / "configuration.json").read_text())
    operational = {"scope", "check_root", "library", "cpus", "memory", "time_limit",
                   "partitions", "exclude_nodes", "account", "shared_cache", "checkpoint", "max_restarts", "save_fits"}
    scientific = lambda value: {key: item for key, item in value.items() if key not in operational}
    if scientific(previous) != scientific(current):
        raise ValueError("Cannot import results from a different scientific configuration")
    with (previous_root / "manifest.csv").open() as stream:
        previous_tasks = {row["task_id"]: row for row in csv.DictReader(stream)}
    with (root / "manifest.csv").open() as stream:
        current_tasks = {item["task_id"]: item for item in csv.DictReader(stream)}
    imported = []
    paths = []
    for row in rows:
        task = str(row["task_id"])
        directory = previous_root / "tasks" / task
        if not (directory / "COMPLETE").is_file():
            continue
        if previous_tasks.get(task) != current_tasks[task]:
            raise ValueError("Completed task identity mismatch: " + task)
        if (directory / "status.txt").read_text().strip() != "COMPLETE":
            raise ValueError("Completed task has inconsistent status: " + task)
        required = [directory / name for name in ("COMPLETE", "results.csv", "data_sha256.txt")]
        if not all(path.is_file() for path in required):
            raise ValueError("Incomplete committed result: " + task)
        receipt = previous_root / "submissions" / (task + ".json")
        backend = "slurm"
        if receipt.is_file():
            shutil.copy2(receipt, root / "submissions" / receipt.name)
        else:
            receipt = previous_root / "local_validation.json"
            if not receipt.is_file():
                raise ValueError("Completed import lacks a Slurm receipt or local validation record: " + task)
            local = json.loads(receipt.read_text())
            if (local.get("backend") != "local" or local.get("exit_code") != 0 or
                    task not in [str(value) for value in local.get("completed_tasks", [])]):
                raise ValueError("Invalid local validation record for imported task: " + task)
            backend = "local"
        (root / "tasks" / task).symlink_to(directory, target_is_directory=True)
        imported.append(dict(task_id=task, directory=str(directory), receipt=str(receipt), backend=backend))
        paths += required + [receipt]
    report = root / "imported_results.json"
    write_json(report, dict(source=str(previous_root), tasks=imported,
        note="Completed-result reuse; original libraries and Slurm or local execution records retained."))
    return paths + [report]


def verify(root):
    for name, expected in json.loads((root / "provenance.json").read_text()).items():
        if digest(Path(name)) != expected:
            raise ValueError(f"Changed input/library/workflow: {name}")
    return json.loads((root / "configuration.json").read_text())


def submission_partitions(root, configuration):
    """Read an optional operational override without changing frozen inputs."""
    override = root / "routing.json"
    if override.is_file():
        routing = json.loads(override.read_text())
        if (not isinstance(routing, dict) or "partitions" not in routing or
                not set(routing).issubset({"partitions", "exclude_nodes"})):
            raise ValueError("routing.json permits only the partitions and exclude_nodes fields")
        partitions = routing["partitions"]
    else:
        partitions = configuration["partitions"]
    if not isinstance(partitions, str) or not re.fullmatch(r"[A-Za-z0-9_-]+(?:,[A-Za-z0-9_-]+)*", partitions):
        raise ValueError("partitions must be a comma-separated list of partition names")
    names = partitions.split(",")
    if len(names) != len(set(names)):
        raise ValueError("Duplicate partitions are not allowed")
    if "preempt" in names and not configuration.get("checkpoint"):
        raise ValueError("preempt requires checkpoint/requeue support")
    return partitions


def submission_excluded_nodes(root, configuration):
    override = root / "routing.json"
    routing = json.loads(override.read_text()) if override.is_file() else {}
    excluded = routing.get("exclude_nodes", configuration.get("exclude_nodes", ""))
    if not isinstance(excluded, str) or (excluded and not re.fullmatch(
            r"[A-Za-z0-9_-]+(?:,[A-Za-z0-9_-]+)*", excluded)):
        raise ValueError("exclude_nodes must be a comma-separated list of node names")
    return excluded


def submit(args):
    if not 1 <= args.max_jobs <= MAX_SUBMISSION_BATCH:
        raise ValueError(f"max_jobs must be in [1, {MAX_SUBMISSION_BATCH}]")
    root = scratch_path(args.output)
    configuration = verify(root)
    partitions = submission_partitions(root, configuration)
    excluded_nodes = submission_excluded_nodes(root, configuration)
    with (root / ".submit.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with (root / "manifest.csv").open() as stream:
            rows = list(csv.DictReader(stream))
        submitted = 0
        for row in rows:
            stop_requested = getattr(args, "stop_requested", None)
            if stop_requested is not None and stop_requested():
                break
            task = row["task_id"]
            receipt = root / "submissions" / f"{task}.json"
            if receipt.exists() or (root / "tasks" / task).exists():
                continue
            command = ["sbatch", "--parsable",
                       "--requeue" if configuration.get("checkpoint") else "--no-requeue",
                       "--open-mode=append", "--nodes=1", "--ntasks=1",
                       f"--cpus-per-task={configuration['cpus']}",
                       f"--mem={configuration['memory']}",
                       f"--time={configuration['time_limit']}",
                       f"--partition={partitions}",
                       f"--account={configuration['account']}",
                       "--comment=" + submission_comment(root, task),
                       f"--job-name=roce{configuration.get('crossfit_layers', 3)}_{row['config']}_k{row['K']}_r{row['sim_id']}_t{task}",
                       f"--chdir={root}", f"--output={root}/logs/{task}_%j.out",
                       f"--error={root}/logs/{task}_%j.err",
                       str(root / "workflow" / "run_repeat_pilot.sh"), str(root), task,
                       configuration.get("worker_script", "run_repeat_pilot.R")]
            if configuration.get("checkpoint"):
                command.insert(2, "--signal=B:USR1@120")
            if excluded_nodes:
                command.insert(2, "--exclude=" + excluded_nodes)
            if args.test_only:
                result = subprocess.run(command[:1] + ["--test-only"] + command[1:],
                                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, timeout=30,
                                        env=submission_environment())
                print(result.stdout + result.stderr, end="")
                result.check_returncode()
            else:
                # Persist intent before sbatch. An interrupted/ambiguous submission
                # is never retried automatically and cannot silently duplicate a seed.
                replace_json(receipt, dict(state="submitting", command=command, task=row,
                                        time=datetime.now(timezone.utc).isoformat()))
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, timeout=30,
                                        env=submission_environment())
                replace_json(receipt, dict(state="submitted" if result.returncode == 0 else "rejected",
                    command=command, task=row, stdout=result.stdout,
                    stderr=result.stderr, returncode=result.returncode))
                result.check_returncode()
                if not re.fullmatch(r"\d+(?:;[^\s]+)?", result.stdout.strip()):
                    raise ValueError("Ambiguous sbatch response; inspect receipt before any retry")
                print(f"Task {task}: job {result.stdout.strip()}")
                # Pace scheduler requests independently of worker concurrency.
                time.sleep(0.2)
            submitted += 1
            if submitted == args.max_jobs:
                break
        print(f"{'Validated' if args.test_only else 'Submitted'} {submitted} jobs")


def restore_nuisance_checkpoints(root, task_id, directory):
    """Restore explicitly staged immutable fits; leave restart progress intact.

    Only completed RDS entries are linked. The R cache still verifies their
    exact fitting inputs, installed-code namespace and model checksums.
    """
    source = root / "checkpoint_seeds" / str(task_id)
    destination = directory / "checkpoints"
    if not source.is_dir() or destination.exists():
        return
    staging = directory / (".checkpoint_restore." + uuid4().hex)
    staging.mkdir()
    try:
        for path in source.rglob("*.rds"):
            target = staging / path.relative_to(source)
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(str(path), str(target))
        staging.replace(destination)
    finally:
        if staging.exists():
            shutil.rmtree(str(staging))


def validate_run(args):
    root = scratch_path(args.output)
    configuration = verify(root)
    with (root / "manifest.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    if str(args.task_id) not in {row["task_id"] for row in rows}:
        raise ValueError("Unknown task ID")
    directory = root / "tasks" / str(args.task_id)
    job_id = os.environ.get("SLURM_JOB_ID", "local")
    restart = int(os.environ.get("SLURM_RESTART_COUNT", "0"))
    if restart < 0 or restart > configuration.get("max_restarts", 20):
        raise ValueError("Maximum restart count exceeded; inspect the saved checkpoints")
    fingerprint = digest(root / "provenance.json")
    previous = None
    if directory.exists():
        if not configuration.get("checkpoint") or restart < 1:
            raise FileExistsError("Task already claimed; duplicate execution refused")
        owner = directory / "owner.json"
        if not owner.is_file():
            receipt = json.loads((root / "submissions" / (str(args.task_id) + ".json")).read_text())
            if receipt.get("state") != "submitted" or receipt["stdout"].strip().split(";")[0] != job_id:
                raise ValueError("Cannot prove the interrupted task belongs to this job")
            previous = dict(job_id=job_id, fingerprint=fingerprint, restart=-1)
        else:
            previous = json.loads(owner.read_text())
        if previous["job_id"] != job_id or previous["fingerprint"] != fingerprint:
            raise ValueError("Checkpoint job/configuration identity mismatch")
        if restart <= previous["restart"]:
            raise FileExistsError("This job attempt already claimed the task")
        if (directory / "COMPLETE").is_file():
            print("COMPLETE")
            return
    else:
        directory.mkdir()
    if previous is not None:
        attempt_archive = directory / "attempts" / str(previous["restart"])
        previous_files = list(directory.glob("*_FAILED.txt")) + [directory / "warnings.csv", directory / "status.txt"]
        for path in previous_files:
            if path.is_file():
                attempt_archive.mkdir(parents=True, exist_ok=True)
                shutil.move(str(path), str(attempt_archive / path.name))
    temporary = directory / ("owner." + uuid4().hex + ".tmp")
    write_json(temporary, dict(job_id=job_id, fingerprint=fingerprint, restart=restart))
    temporary.replace(directory / "owner.json")
    restore_nuisance_checkpoints(root, args.task_id, directory)
    if (directory / "INTERRUPTED").exists():
        (directory / "INTERRUPTED").unlink()
    print("RUN")


def submit_supervisor(args):
    root = scratch_path(args.output)
    configuration = verify(root)
    receipt = root / "controller_submission.json"
    partitions = getattr(args, "controller_partitions", DEFAULT_CONTROLLER_PARTITIONS)
    hours = getattr(args, "controller_hours", 1)
    command = ["sbatch", "--parsable", "--requeue", "--signal=B:USR1@120",
               "--open-mode=append", "--nodes=1", "--ntasks=1",
               "--cpus-per-task=1", "--mem=1G", "--time={}:00:00".format(hours),
               "--account=" + configuration["account"], "--partition=" + partitions,
               "--job-name=roce_validation_controller", "--chdir=" + str(root),
               "--output=" + str(root / "controller_%j.out"),
               "--error=" + str(root / "controller_%j.err"),
               str(Path(__file__).resolve().with_name("supervise_repeat_pilot.sh")),
               str(root), str(args.max_in_flight), str(hours - 0.1),
               str(Path(__file__).resolve().parent)]
    excluded_nodes = submission_excluded_nodes(root, configuration)
    if excluded_nodes:
        command.insert(2, "--exclude=" + excluded_nodes)
    write_json(receipt, dict(state="submitting", command=command))
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            universal_newlines=True, timeout=30,
                            env=submission_environment())
    with receipt.open("w") as stream:
        json.dump(dict(state="submitted" if result.returncode == 0 else "rejected", command=command,
                       stdout=result.stdout, stderr=result.stderr, returncode=result.returncode), stream, indent=2)
        stream.write("\n")
    result.check_returncode()
    if not re.fullmatch(r"\d+(?:;[^\s]+)?", result.stdout.strip()):
        raise ValueError("Ambiguous controller submission; inspect receipt before retrying")
    print("Validation controller: " + result.stdout.strip())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action")
    commands.required = True  # MSI's system Python is 3.6.
    prep = commands.add_parser("prepare")
    prep.add_argument("output")
    prep.add_argument("--check-root", required=True)
    prep.add_argument("--completed-from", help="Reuse committed results after an execution-only migration")
    prep.add_argument("--dgp-type", choices=["bounded", "face"], default="bounded")
    prep.add_argument("--dimensions", nargs="+", type=int, choices=[4, 10, 20, 50, 100, 200])
    prep.add_argument("--configs", nargs="+", choices=["C1", "C2", "C3", "C4"], default=["C1", "C2", "C3"],
                      help="C4 is an explicit both-models-misspecified boundary experiment")
    prep.add_argument("--sources", nargs="+", type=int, choices=[2, 4, 6, 8], default=[2])
    prep.add_argument("--n-deviated-sites", type=int, default=1,
                      help="Number of leading sources with outcome deviation; zero rho leaves all outcomes compatible")
    prep.add_argument("--rhos", nargs="+", type=float, choices=[0, .5, 1, 1.5, 2, 2.5], default=[0])
    seeds = prep.add_mutually_exclusive_group()
    seeds.add_argument("--repeats", nargs="+", type=int)
    seeds.add_argument("--n-repeats", type=int, choices=range(1, 2001))
    prep.add_argument("--seed-start", type=int, default=1)
    prep.add_argument("--protocol", choices=["one_round", "two_round"], default="one_round")
    prep.add_argument("--deviation-mechanism", choices=["treated_arm", "both_arms"], default="treated_arm")
    prep.add_argument("--recipes", nargs="+", choices=["legacy", "score_derivative"],
                      default=None)
    prep.add_argument("--covariate-shift", type=float, choices=[0, .5, 1, 2], default=1)
    prep.add_argument("--source-treatment-scale", type=float, choices=[0, 1], default=1)
    prep.add_argument("--save-fits", choices=["none", "first", "all"])
    prep.add_argument("--aggregation-cutoff", type=float, default=1.0,
                      help="Wald soft-penalty cutoff; aggregation lambda is its reciprocal")
    prep.add_argument("--baseline-methods", nargs="+", choices=BASELINE_METHODS,
                      help="Run only these baselines and the paired target-only reference")
    prep.add_argument("--baseline-reference-roots", nargs="+",
                      help="Existing method campaigns whose data hashes and target references must agree")
    prep.add_argument("--target-nuisance-method", choices=["hou_calibrated", "lasso"], default="hou_calibrated")
    prep.add_argument("--source-nuisance-method", choices=["calibrated", "standard"], default="calibrated",
                      help="Source final score calibration or standard outcome/merged-weight fitting")
    prep.add_argument("--source-validation-method", choices=["calibrated", "complete", "initial", "outer_fit"], default="calibrated",
                      help="Complete chosen program, calibrated compatibility alias, or legacy initial-model validation")
    prep.add_argument("--crossfit-layers", type=int, choices=[2, 3],
                      help="Two: eta uses outer-fit training scores; three: independent complete inner validation")
    prep.add_argument("--nuisance-cv-certificate", action="store_true",
                      help="Use the checked density-CV KKT certificate without changing the lambda grid")
    prep.add_argument("--cpus", type=int, choices=range(1, 17), default=2)
    prep.add_argument("--shared-cache", action="store_true",
                      help="Use temporary shared nuisance caches within each repeat")
    prep.add_argument("--memory", default="4G")
    prep.add_argument("--time-limit", default="00:30:00")
    checkpoint_options = prep.add_mutually_exclusive_group()
    checkpoint_options.add_argument("--checkpoint", action="store_true", default=True,
        help="Persist exact nuisance fits and enable Slurm requeue (default)")
    checkpoint_options.add_argument("--no-checkpoint", action="store_false", dest="checkpoint",
        help="Disable persistence for a non-preempting partition")
    prep.add_argument("--max-restarts", type=int, choices=range(0, 101), default=20)
    prep.add_argument("--partitions", default=DEFAULT_PARTITIONS)
    prep.add_argument("--exclude-nodes", default="",
                      help="Comma-separated node names to exclude from workers and controller")
    prep.add_argument("--account", default="hou00123")
    sub = commands.add_parser("submit")
    sub.add_argument("output")
    sub.add_argument("--max-jobs", type=int, choices=range(1, MAX_SUBMISSION_BATCH + 1), default=1)
    sub.add_argument("--test-only", action="store_true")
    run = commands.add_parser("validate-run")
    run.add_argument("output")
    run.add_argument("task_id", type=int)
    supervisor = commands.add_parser("supervise")
    supervisor.add_argument("output")
    supervisor.add_argument("--max-in-flight", type=int, choices=range(1, MAX_IN_FLIGHT + 1),
                            default=DEFAULT_MAX_IN_FLIGHT)
    supervisor.add_argument("--controller-partitions", default=DEFAULT_CONTROLLER_PARTITIONS)
    supervisor.add_argument("--controller-hours", type=int, choices=range(1, 24), default=1)
    args = parser.parse_args()
    {"prepare": prepare, "submit": submit, "validate-run": validate_run,
     "supervise": submit_supervisor}[args.action](args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        sys.exit(str(error))
