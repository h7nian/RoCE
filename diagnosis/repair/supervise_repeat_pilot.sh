#!/bin/bash
set -euo pipefail
PILOT_ROOT="${1:?Pass the frozen validation root}"
CONCURRENCY_ARGS=()
if [[ -n "${2:-}" ]]; then CONCURRENCY_ARGS=(--max-in-flight "$2"); fi
HOURS="${3:-0.9}"
# Slurm copies this shell script into its spool directory; $0 is not the
# directory containing the frozen Python workflow.
WORKFLOW_ROOT="${4:-${PILOT_ROOT}/workflow}"
export PYTHONDONTWRITEBYTECODE=1
python3 "${WORKFLOW_ROOT}/supervise_repeat_pilot.py" "$PILOT_ROOT" \
  "${CONCURRENCY_ARGS[@]}" --hours "$HOURS" &
controller_pid=$!
request_stop() { kill -USR1 "$controller_pid" 2>/dev/null || true; }
trap request_stop USR1
trap 'kill -TERM "$controller_pid" 2>/dev/null || true; wait "$controller_pid" 2>/dev/null || true; exit 143' TERM
set +e
while true; do
    wait "$controller_pid"
    status=$?
    if ! kill -0 "$controller_pid" 2>/dev/null; then break; fi
done
set -e
if [[ "$status" == 75 ]]; then
    scontrol requeue "${SLURM_JOB_ID:?Requeue requires a Slurm controller job}"
    exit 0
fi
exit "$status"
