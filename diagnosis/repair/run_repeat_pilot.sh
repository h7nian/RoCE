#!/bin/bash
# Exactly one repeat per job. Requeue resumes only the same validated job ID.
set -euo pipefail
PILOT_ROOT="${1:?Pass the pilot root}"
TASK_ID="${2:?Pass one manifest task ID}"
WORKER_SCRIPT="${3:-run_repeat_pilot.R}"
TASK_ROOT="${PILOT_ROOT}/tasks/${TASK_ID}"
claim=$(python3 -B "${PILOT_ROOT}/workflow/submit_repeat_pilot.py" validate-run "$PILOT_ROOT" "$TASK_ID")
if [[ "$claim" == COMPLETE ]]; then exit 0; fi
exec 9>"${TASK_ROOT}/worker.lock"
flock -n 9 || { echo "Another worker holds this repeat" >&2; exit 1; }
if [[ -f "${TASK_ROOT}/COMPLETE" ]]; then exit 0; fi
mark_interrupted() {
    touch "${TASK_ROOT}/INTERRUPTED"
}
request_requeue() {
    mark_interrupted
    scontrol requeue "${SLURM_JOB_ID:?Requeue requires a Slurm job}"
}
trap 'request_requeue' USR1
trap 'mark_interrupted; exit 143' TERM
export TMPDIR="${TASK_ROOT}/tmp"
mkdir -p "$TMPDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_NUISANCE_CV_THREADS=1
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
Rscript --vanilla "${PILOT_ROOT}/workflow/${WORKER_SCRIPT}" "$PILOT_ROOT" "$TASK_ID" &
worker_pid=$!
set +e
while true; do
    wait "$worker_pid"
    status=$?
    if ! kill -0 "$worker_pid" 2>/dev/null; then break; fi
done
exit "$status"
