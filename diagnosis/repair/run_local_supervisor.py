#!/usr/bin/env python3
"""Run the existing monitor within an allocation, with a queued Slurm fallback."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("workflow", type=Path)
    parser.add_argument("--hours", type=float, default=1.25)
    args = parser.parse_args()
    root, workflow = args.root.resolve(), args.workflow.resolve()
    scratch = Path("/scratch.global/zhan9381/FACE-HD")
    if scratch not in root.parents or scratch not in workflow.parents:
        parser.error("Use FACE-HD scratch directories")
    allocation = os.environ.get("SLURM_JOB_ID", "")
    if not allocation.isdigit() or not 0 < args.hours <= 8:
        parser.error("An existing Slurm allocation and a bounded duration are required")
    def run(command):
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, timeout=30)
        result.check_returncode()
        return result.stdout
    description = run(["scontrol", "-o", "show", "job", allocation])
    end = re.search(r"(?:^| )EndTime=(\S+)", description)
    if end is None or "JobState=RUNNING" not in description:
        raise ValueError("The local allocation is not running")
    available = (datetime.strptime(end.group(1), "%Y-%m-%dT%H:%M:%S") - datetime.now()).total_seconds()
    hours = min(args.hours, (available - 300) / 3600)
    if hours <= 0:
        raise ValueError("Insufficient time remains in the current allocation")
    receipt = json.loads((root / "controller_submission.json").read_text())
    backup = receipt["stdout"].strip().split(";")[0]
    if not backup.isdigit() or run(["squeue", "-h", "-j", backup, "-o", "%T"]).strip() != "PENDING":
        raise ValueError("The backup controller must still be pending")
    run(["scontrol", "update", "JobId=" + backup, "Dependency=afterany:" + allocation])
    command = [sys.executable, "-B", str(workflow / "supervise_repeat_pilot.py"),
               str(root), "--max-in-flight", "512", "--hours", str(hours)]
    metadata = dict(pid=os.getpid(), node=socket.gethostname(), allocation=allocation,
                    backup_job=backup, command=command, state="RUNNING")
    (root / "local_controller.json").write_text(json.dumps(metadata, indent=2) + "\n")
    result = subprocess.run(command)
    status_file = root / "supervisor_status.json"
    status = json.loads(status_file.read_text()).get("state") if status_file.exists() else "UNKNOWN"
    if status in {"COMPLETE", "NEEDS_REVIEW"}:
        run(["scancel", backup])
    else:
        run(["scontrol", "update", "JobId=" + backup, "Dependency="])
    metadata.update(state="EXITED", returncode=result.returncode, campaign_state=status)
    (root / "local_controller.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
