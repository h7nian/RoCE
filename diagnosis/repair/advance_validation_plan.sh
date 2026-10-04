#!/bin/bash
# Release gated campaigns, preserving their submission receipts across requeue.
set -euo pipefail
PLAN="${1:?Pass the prepared expansion plan}"
WORKFLOW="${2:?Pass the frozen coordinator workflow}"
python3 - "$WORKFLOW" <<'PY'
import hashlib, json, sys
from pathlib import Path
root = Path(sys.argv[1])
for name, expected in json.loads((root / "manifest.json").read_text()).items():
    if hashlib.sha256(Path(name).read_bytes()).hexdigest() != expected:
        raise SystemExit("Changed coordinator input: " + name)
PY
set +e
python3 -B "$WORKFLOW/advance_validation_plan.py" "$PLAN" --hours 0.8
status=$?
set -e
if [[ "$status" == 75 ]]; then
    scontrol requeue "${SLURM_JOB_ID:?Requeue requires a Slurm coordinator job}"
    exit 0
fi
exit "$status"
