#!/usr/bin/env python3
"""Submit one repeat per job in bounded waves and retain scheduler failures.

The first seed of every cell must finish before the remaining repeats are
submitted. A failed worker stops further submissions; existing workers finish.
No failed or ambiguous submission is retried automatically.
"""
import argparse
import csv
import fcntl
import json
import os
from pathlib import Path
import subprocess
import signal
import sys
import time
from types import SimpleNamespace
from datetime import datetime, timedelta

sys.dont_write_bytecode = True
import submit_repeat_pilot as pilot
import summarize_repeat_pilot as review

RESULT_COMMIT_GRACE_SECONDS = 120
MAX_SCHEDULER_POLL_FAILURES = 3
STOP_SIGNAL = None


def request_stop(signum, frame):
    global STOP_SIGNAL
    STOP_SIGNAL = signum


TERMINAL = {"COMPLETED", "FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY",
            "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE", "REVOKED"}


def recover_submission(root, path, receipt):
    if receipt.get("state") != "submitting":
        raise ValueError("Rejected submission requires review: " + str(path))
    comments = [part.split("=", 1)[1] for part in receipt.get("command", [])
                if part.startswith("--comment=")]
    token = comments[0] if comments else pilot.submission_comment(root, path.stem)
    start = (datetime.fromtimestamp((root / "configuration.json").stat().st_mtime)
             - timedelta(days=1)).strftime("%Y-%m-%d")
    commands = [
        ["squeue", "-h", "-u", str(os.getuid()), "-o", "%i|%k"],
        ["sacct", "-n", "-X", "-P", "-u", str(os.getuid()),
         "--starttime=" + start, "--format=JobIDRaw,Comment"]]
    matches = set()
    for command in commands:
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                universal_newlines=True, timeout=30)
        result.check_returncode()
        for line in result.stdout.splitlines():
            fields = line.split("|")
            if len(fields) >= 2 and fields[0].isdigit() and fields[1] == token:
                matches.add(fields[0])
    if len(matches) != 1:
        raise ValueError("Ambiguous interrupted submission; refusing to duplicate task " + path.stem)
    receipt.update(state="submitted", stdout=matches.pop() + "\n", returncode=0,
                   stderr="Recovered the original job ID after controller interruption")
    pilot.replace_json(path, receipt)
    return receipt


def submitted_jobs(root):
    jobs = {}
    for path in sorted((root / "submissions").glob("*.json")):
        receipt = json.loads(path.read_text())
        if receipt.get("state") != "submitted":
            receipt = recover_submission(root, path, receipt)
        job = receipt["stdout"].strip().split(";")[0]
        if not job.isdigit():
            raise ValueError("Invalid job receipt: " + str(path))
        jobs[path.stem] = job
    return jobs


def scheduler_states(job_ids, start_date):
    if not job_ids:
        return {}
    result = subprocess.run(["sacct", "-n", "-X", "-P", "-j", ",".join(job_ids),
        "--starttime=" + start_date, "--format=JobIDRaw,State,Elapsed,ExitCode"], stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, universal_newlines=True, timeout=60)
    result.check_returncode()
    states = {}
    for line in result.stdout.splitlines():
        fields = line.split("|")
        if len(fields) >= 4 and fields[0] in job_ids:
            state = fields[1].split()[0].rstrip("+")
            states[fields[0]] = dict(job_id=fields[0], scheduler_state=state,
                                    elapsed=fields[2], exit_code=fields[3])
    # The live queue is authoritative during a requeue transition; accounting
    # can briefly retain PREEMPTED/NODE_FAIL while the same job is pending again.
    live = subprocess.run(["squeue", "-h", "-u", str(os.getuid()), "-o", "%i|%T"],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, timeout=30)
    live.check_returncode()
    for line in live.stdout.splitlines():
        fields = line.split("|")
        if len(fields) == 2 and fields[0] in job_ids:
            previous = states.get(fields[0], dict(job_id=fields[0], elapsed="", exit_code=""))
            states[fields[0]] = dict(previous, scheduler_state=fields[1])
    return states


def submission_capacity(rows, jobs, terminal, maximum, first_seed, completed_tasks):
    first_wave = {row["task_id"] for row in rows if int(row["sim_id"]) == first_seed}
    in_flight = sum(job not in terminal for job in jobs.values())
    remaining = sum(row["task_id"] not in jobs and row["task_id"] not in completed_tasks for row in rows)
    if not first_wave.issubset(completed_tasks):
        remaining = len(first_wave.difference(jobs).difference(completed_tasks))
    return max(0, min(maximum - in_flight, remaining))


def supervise(root, maximum, hours, interval):
    configuration = pilot.verify(root)
    with (root / "manifest.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    started = time.monotonic()
    start_date = (datetime.fromtimestamp((root / "configuration.json").stat().st_mtime)
                  - timedelta(days=1)).strftime("%Y-%m-%d")
    terminal, observed = {}, {}
    terminal_seen_at = {}
    terminal_states = TERMINAL - ({"PREEMPTED", "NODE_FAIL"} if configuration.get("checkpoint") else set())
    failed = False
    scheduler_poll_failures = 0
    with (root / ".supervisor.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        while True:
            jobs = submitted_jobs(root)
            unknown = [job for job in jobs.values() if job not in terminal]
            try:
                current_states = scheduler_states(unknown, start_date)
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError) as error:
                # A failed read is not evidence that a worker terminated. Do
                # not submit anything until BOTH accounting and live queue
                # queries succeed; existing receipts remain authoritative.
                scheduler_poll_failures += 1
                stopping = (STOP_SIGNAL is not None or
                            time.monotonic() - started >= hours * 3600 or
                            scheduler_poll_failures >= MAX_SCHEDULER_POLL_FAILURES)
                state = ("STOPPED" if STOP_SIGNAL == signal.SIGTERM else
                         "CHECKPOINTED_FOR_REQUEUE" if stopping else "WAITING_FOR_SCHEDULER")
                pilot.replace_json(root / "supervisor_status.json", dict(
                    state=state, planned=len(rows), submitted=len(jobs),
                    completed=sum((root / "tasks" / row["task_id"] / "COMPLETE").is_file()
                                  for row in rows), next_batch=0,
                    scheduler_poll_failures=scheduler_poll_failures,
                    scheduler_error=str(error)))
                print("{}: {}".format(state, error), flush=True)
                if stopping:
                    return None
                time.sleep(interval)
                continue
            scheduler_poll_failures = 0
            observed.update(current_states)
            terminal.update({job: value for job, value in observed.items()
                             if value["scheduler_state"] in terminal_states})
            now = time.monotonic()
            for job in terminal:
                terminal_seen_at.setdefault(job, now)
            ledger = []
            for task, job in jobs.items():
                item = observed.get(job, dict(job_id=job, scheduler_state="ACCOUNTING_PENDING",
                                             elapsed="", exit_code=""))
                ledger.append(dict(task_id=task, **item))
            review.write_csv(root / "scheduler_status.csv", ledger)
            completed = {row["task_id"] for row in rows
                         if (root / "tasks" / row["task_id"] / "COMPLETE").is_file()}
            # Worker status may briefly reflect a killed child during preemption.
            # Only terminal scheduler failures stop the batch.
            failed = failed or any(value["scheduler_state"] != "COMPLETED"
                                   for value in terminal.values())
            # Shared-filesystem metadata can lag accounting. Require the final
            # result marker, but allow a bounded visibility delay after Slurm exit.
            awaiting_results = {task for task, job in jobs.items()
                if job in terminal and terminal[job]["scheduler_state"] == "COMPLETED"
                and task not in completed}
            overdue_results = {task for task in awaiting_results
                if now - terminal_seen_at[jobs[task]] >= RESULT_COMMIT_GRACE_SECONDS}
            failed = failed or bool(overdue_results)
            if failed:
                capacity = 0
                state = "DRAINING_AFTER_FAILURE"
            else:
                capacity = submission_capacity(rows, jobs, terminal, maximum,
                                               configuration["first_seed"], completed)
                state = "RUNNING"
            pilot.replace_json(root / "supervisor_status.json", dict(
                state=state, planned=len(rows), submitted=len(jobs), completed=len(completed),
                terminal=len(terminal), next_batch=capacity, awaiting_results=len(awaiting_results),
                overdue_results=sorted(overdue_results)))
            print("planned={} submitted={} complete={} terminal={} state={}".format(
                len(rows), len(jobs), len(completed), len(terminal), state), flush=True)
            if len(terminal) == len(jobs) and (failed or (len(completed) == len(rows) and not awaiting_results)):
                passed = review.summarize(root, root / "review")
                pilot.replace_json(root / "supervisor_status.json", dict(
                    state="COMPLETE" if passed else "NEEDS_REVIEW", planned=len(rows),
                    submitted=len(jobs), completed=len(completed)))
                return passed
            if STOP_SIGNAL is not None or time.monotonic() - started >= hours * 3600:
                review.summarize(root, root / "review")
                pilot.replace_json(root / "supervisor_status.json", dict(
                    state="STOPPED" if STOP_SIGNAL == signal.SIGTERM else "CHECKPOINTED_FOR_REQUEUE", planned=len(rows), submitted=len(jobs),
                    completed=len(completed)))
                return None
            if capacity:
                pilot.submit(SimpleNamespace(output=str(root),
                    max_jobs=min(capacity, pilot.MAX_SUBMISSION_BATCH), test_only=False,
                    stop_requested=lambda: STOP_SIGNAL is not None))
            time.sleep(interval)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--max-in-flight", type=int, choices=range(1, pilot.MAX_IN_FLIGHT + 1),
                        default=pilot.DEFAULT_MAX_IN_FLIGHT)
    parser.add_argument("--hours", type=float, default=0.9)
    parser.add_argument("--poll-seconds", type=int, choices=range(10, 61), default=30)
    args = parser.parse_args()
    if not 0 < args.hours <= 95:
        parser.error("--hours must be in (0, 95]")
    root = pilot.scratch_path(args.output)
    signal.signal(signal.SIGUSR1, request_stop)
    signal.signal(signal.SIGTERM, request_stop)
    try:
        success = supervise(root, args.max_in_flight, args.hours, args.poll_seconds)
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        (root / "SUPERVISOR_FAILED.txt").write_text(str(error) + "\n")
        raise
    raise SystemExit(143 if STOP_SIGNAL == signal.SIGTERM else
                     (75 if success is None else (0 if success else 1)))
