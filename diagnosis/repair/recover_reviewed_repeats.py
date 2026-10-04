#!/usr/bin/env python3
"""Retry explicitly reviewed infrastructure failures, retaining every attempt.

This is an operator-invoked recovery tool, not an automatic failure policy.
The audit supplies exact task/job identities; numerical failures are refused.
"""
import argparse
from contextlib import ExitStack
import fcntl
import json
import os
import re
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot

ALLOWED_FAILURES = {"inherited_restart_counter", "node_startup_or_io",
                    "scheduler_startup", "pre_R_startup_reviewed", "reviewed_R_file_io"}
TERMINAL_FAILURES = {"FAILED", "NODE_FAIL", "BOOT_FAIL", "PREEMPTED"}


def run(command):
    return subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          universal_newlines=True, timeout=60, check=True).stdout


def restore_seed_inventory(root, task, archive):
    cache = archive / "checkpoints"
    if not cache.is_dir():
        return
    seed = root / "checkpoint_seeds" / task
    seed.parent.mkdir(exist_ok=True)
    if seed.exists():
        # A restarted worker's cache usually contains the prior seed inventory.
        # Retain any earlier entries and refuse conflicting finalized payloads.
        for path in seed.rglob("*.rds"):
            target = cache / path.relative_to(seed)
            if target.exists():
                if not os.path.samefile(str(path), str(target)) and pilot.digest(path) != pilot.digest(target):
                    raise ValueError("Conflicting immutable checkpoint: " + str(target))
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                os.link(str(path), str(target))
        seed.rename(archive / "checkpoint_seed_before_recovery")
    seed.symlink_to(cache, target_is_directory=True)



def validate_r_failure_review(directory, item):
    """Admit only individually reviewed, hash-pinned file-connection errors.

    A generic R failure stays ineligible. The audit must identify every error
    marker and corroborating I/O logs; node/time correlation is reviewed by the
    operator before invoking this tool, not inferred automatically here.
    """
    failures = sorted(directory.glob("*_FAILED.txt"))
    if not failures:
        if item["failure_kind"] == "reviewed_R_file_io":
            raise ValueError("Reviewed R error markers are missing")
        return
    reviewed = item.get("reviewed_error_files", {})
    if (item["failure_kind"] != "reviewed_R_file_io" or
            set(reviewed) != {path.name for path in failures}):
        raise ValueError("R failure requires separate review: " + str(failures))
    prefix = r"^parallel_lapply: worker\(s\) for item\(s\) [a-zA-Z0-9_, ]+ failed: "
    pattern = (r"(?:Error in (?:gzfile|file)\([^\n]*\) : cannot open the connection|"
               r"Error : Could not atomically store the shared nuisance fit\.)")
    for path in failures:
        message = re.sub(prefix, "", path.read_text().strip())
        errors = message.split(" | ")
        if (pilot.digest(path) != reviewed[path.name] or
                not all(re.fullmatch(pattern, error) for error in errors)):
            raise ValueError("R failure changed or is not a reviewed file-connection error")
    evidence = item.get("io_evidence", [])
    if not evidence:
        raise ValueError("Reviewed R file errors require corroborating I/O evidence")
    for record in evidence:
        path = Path(record["path"])
        if (path.suffix != ".err" or pilot.digest(path) != record["sha256"] or
                not re.search(r"Input/output error|I/O error", path.read_text())):
            raise ValueError("Corroborating I/O evidence changed or does not show an I/O failure")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    audit_path = pilot.scratch_path(args.audit)
    tasks = json.loads(audit_path.read_text())["tasks"]
    tasks = [item for item in tasks if not item["live"]]
    roots = sorted({pilot.scratch_path(item["root"]) for item in tasks})
    with ExitStack() as locks:
        for root in roots:
            pilot.verify(root)
            for name in (".supervisor.lock", ".submit.lock"):
                lock = locks.enter_context((root / name).open("a"))
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        live = {line.split("|")[0] for line in run(["squeue", "-h", "-u", str(os.getuid()), "-o", "%i|%T"]).splitlines()}
        accounting = {}
        if tasks:
            output = run(["sacct", "-n", "-X", "-P", "-j", ",".join(item["job_id"] for item in tasks),
                          "--format=JobIDRaw,State,ExitCode"])
            accounting = {parts[0]: parts[1:] for parts in (line.split("|") for line in output.splitlines())}
        for item in tasks:
            root, task, job = Path(item["root"]), item["task_id"], item["job_id"]
            directory = root / "tasks" / task
            receipt_path = root / "submissions" / (task + ".json")
            receipt = json.loads(receipt_path.read_text())
            if receipt["stdout"].strip().split(";")[0] != job or job in live:
                raise ValueError("Job ownership changed or still active: " + str(receipt_path))
            if (directory / "COMPLETE").exists() or item["failure_kind"] not in ALLOWED_FAILURES:
                raise ValueError("Task is complete or failure needs individual review: " + str(directory))
            if accounting.get(job, [""])[0] not in TERMINAL_FAILURES:
                raise ValueError("Missing terminal infrastructure-failure accounting: " + job)
            validate_r_failure_review(directory, item)
        print("Verified {} reviewed failures; execute={}".format(len(tasks), args.execute), flush=True)
        if not args.execute:
            return
        for item in tasks:
            root, task, job = Path(item["root"]), item["task_id"], item["job_id"]
            directory = root / "tasks" / task
            receipt_path = root / "submissions" / (task + ".json")
            receipt = json.loads(receipt_path.read_text())
            archive = root / "task_attempts" / task / job
            archive.parent.mkdir(parents=True, exist_ok=True)
            if directory.exists():
                with (directory / "worker.lock").open("a") as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    directory.rename(archive)
            else:
                archive.mkdir()
            pilot.write_json(archive / "recovery_review.json", item)
            restore_seed_inventory(root, task, archive)
            history = root / "submission_attempts" / task / (job + ".json")
            history.parent.mkdir(parents=True, exist_ok=True)
            receipt_path.rename(history)
            configuration = json.loads((root / "configuration.json").read_text())
            command = [part for part in receipt["command"]
                       if not part.startswith(("--exclude=", "--comment=", "--partition="))]
            command[2:2] = ["--exclude=" + pilot.submission_excluded_nodes(root, configuration),
                             "--partition=" + pilot.submission_partitions(root, configuration),
                             "--comment=" + pilot.submission_comment(root, task) + ":retry_" + job]
            intent = dict(state="submitting", command=command, task=receipt["task"],
                          previous_job_id=job, reviewed_failure=item["failure_kind"], archive=str(archive))
            pilot.replace_json(receipt_path, intent)
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    universal_newlines=True, timeout=60, env=pilot.submission_environment())
            intent.update(state="submitted" if result.returncode == 0 else "rejected",
                          stdout=result.stdout, stderr=result.stderr, returncode=result.returncode)
            pilot.replace_json(receipt_path, intent)
            result.check_returncode()
            if not result.stdout.strip().split(";")[0].isdigit():
                raise ValueError("Ambiguous recovery submission: " + str(receipt_path))
            print("{} task {}: {} -> {}".format(root.name, task, job, result.stdout.strip()), flush=True)
            time.sleep(.2)


if __name__ == "__main__":
    main()
