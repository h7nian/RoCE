#!/bin/bash
# Arguments match run_honest_population_case.R; one repeat per Slurm job.
set -euo pipefail
case_output="${3:?Pass the scratch case directory}"
case_workflow="${2:?Pass the frozen workflow directory}"
[[ "$case_output" == /scratch.global/zhan9381/FACE-HD/* ]] || exit 2
mkdir -p "$case_output"
exec 9>"$case_output/worker.lock"
flock -n 9 || { echo "Another worker holds this research repeat" >&2; exit 1; }
request_requeue() {
    if (( ${SLURM_RESTART_COUNT:-0} >= 20 )); then exit 1; fi
    scontrol requeue "${SLURM_JOB_ID:?Requeue requires a Slurm job}"
}
trap 'request_requeue' USR1
trap 'exit 143' TERM
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export TMPDIR=/scratch.global/zhan9381/tmp
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 ROCE_NUISANCE_CV_THREADS=1
Rscript --vanilla "$case_workflow/run_honest_population_case.R" "$@" &
case_pid=$!
set +e
while true; do
    wait "$case_pid"
    case_status=$?
    if ! kill -0 "$case_pid" 2>/dev/null; then break; fi
done
exit "$case_status"
