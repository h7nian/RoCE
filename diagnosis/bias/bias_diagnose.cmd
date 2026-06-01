#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --job-name=FACE-bias-diag
#SBATCH --output=diagnosis/bias/log/%x_%A_%a.out
#SBATCH --error=diagnosis/bias/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5
#SBATCH --requeue
#SBATCH --signal=B:USR1@120
# Other --time/--mem/--cpus-per-task/--partition/--array come from
# bias_diagnose.sh via the sbatch CLI.

# ============================================================================
# FACE-HD bias-diagnosis array script
# ----------------------------------------------------------------------------
# Each array task executes ONE setting from the table below by calling
# diagnosis/bias/run_setting.R. All tasks are independent; results are
# written to diagnosis/bias/results/<tag>_..._summary.csv.
#
# Checkpoint / Requeue (mirrors main.cmd):
#   1. SLURM sends SIGUSR1 120 s before termination on preempt-class partitions.
#   2. handle_preemption() touches CHECKPOINT_DIR/.preempt_signal_${SLURM_JOB_ID};
#      R's run_simulation_study() polls this file, saves an .RData checkpoint
#      via save_checkpoint(), and creates .checkpoint_saved_${SLURM_JOB_ID}.
#   3. The handler waits up to ~100 s for that flag, then calls
#      `scontrol requeue ${SLURM_JOB_ID}` to relaunch the SAME array task with
#      the SAME SLURM_JOB_ID. SLURM array tasks keep their JOB_ID across
#      requeue, so the checkpoint path is stable.
#   4. On restart, load_checkpoint() resumes the sim loop from the saved index.
#   5. cleanup_checkpoint() in run_setting.R wipes artifacts on success.
#
# THE ALGORITHM IS NOT MODIFIED — this is a measurement-only suite.
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"

mkdir -p diagnosis/bias/log diagnosis/bias/results diagnosis/bias/checkpoints diagnosis/bias/tmp

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NLAMBDA_INIT=100

N_SIMS="${DIAG_N_SIMS:-100}"

# ----------------------------------------------------------------------------
# Checkpoint directory + job id are exported so run_setting.R sees them.
# Using a per-task subdir keeps the array tasks' checkpoints isolated.
# ----------------------------------------------------------------------------
export CHECKPOINT_DIR="diagnosis/bias/checkpoints/task_${SLURM_ARRAY_TASK_ID:-na}"
export CHECKPOINT_JOB_ID="${SLURM_JOB_ID}"
mkdir -p "${CHECKPOINT_DIR}"

# ----------------------------------------------------------------------------
# Preemption handler (mirrors main.cmd).
#   - touch .preempt_signal_<JOBID> to ask the R loop to flush.
#   - wait for .checkpoint_saved_<JOBID> (max ~100 s).
#   - scontrol requeue, same JOB_ID -> resume next attempt.
# ----------------------------------------------------------------------------
handle_preemption() {
    echo "$(date): Received SIGUSR1. Asking R to checkpoint..."
    mkdir -p "${CHECKPOINT_DIR}"
    touch "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}"
    for _ in $(seq 1 100); do
        if [[ -f "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}" ]]; then
            echo "$(date): Checkpoint saved successfully."
            rm -f "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}"
            break
        fi
        sleep 1
    done
    rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}"
    echo "$(date): Requesting requeue (JOB_ID=${SLURM_JOB_ID})..."
    if ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo "$(date): WARNING: scontrol requeue failed; exiting." >&2
    fi
    exit 0
}
trap 'handle_preemption' USR1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

# ----------------------------------------------------------------------------
# Grid definition. Columns: tag  n_total  K  p  shift_strength
# Anchor point is task 1 (n=1000, K=3, p=20, shift=0.5); other tasks vary
# one axis at a time. Keep in sync with the comment block in bias_diagnose.sh.
#
# Revision after the first run (job 9276320):
#   * Dropped K_2: the FACE-HD two-layer crossfit is degenerate at K=2 (see
#     diagnosis/bias/results/K_2_n1000_K2_p20_ss0.5_summary.csv; the
#     one/two-round estimators diverge with bias ~±100).
#   * Replaced shift_0 (shift=0.0, rejected by validate_simulation_params)
#     with shift_005 (shift=0.05) — keeps the bias-vs-shift axis meaningful
#     while passing the positivity check.
# ----------------------------------------------------------------------------
declare -A GRID_TAG GRID_N GRID_K GRID_P GRID_SHIFT
GRID_TAG[1]="anchor";    GRID_N[1]=1000;  GRID_K[1]=3;  GRID_P[1]=20;   GRID_SHIFT[1]=0.5
GRID_TAG[2]="n_small";   GRID_N[2]=500;   GRID_K[2]=3;  GRID_P[2]=20;   GRID_SHIFT[2]=0.5
GRID_TAG[3]="n_large";   GRID_N[3]=2000;  GRID_K[3]=3;  GRID_P[3]=20;   GRID_SHIFT[3]=0.5
GRID_TAG[4]="p_50";      GRID_N[4]=1000;  GRID_K[4]=3;  GRID_P[4]=50;   GRID_SHIFT[4]=0.5
GRID_TAG[5]="p_100";     GRID_N[5]=1000;  GRID_K[5]=3;  GRID_P[5]=100;  GRID_SHIFT[5]=0.5
GRID_TAG[6]="p_200";     GRID_N[6]=1000;  GRID_K[6]=3;  GRID_P[6]=200;  GRID_SHIFT[6]=0.5
GRID_TAG[7]="K_4";       GRID_N[7]=1000;  GRID_K[7]=4;  GRID_P[7]=20;   GRID_SHIFT[7]=0.5
GRID_TAG[8]="K_5";       GRID_N[8]=1000;  GRID_K[8]=5;  GRID_P[8]=20;   GRID_SHIFT[8]=0.5
GRID_TAG[9]="shift_005"; GRID_N[9]=1000;  GRID_K[9]=3;  GRID_P[9]=20;   GRID_SHIFT[9]=0.05
GRID_TAG[10]="shift_1";  GRID_N[10]=1000; GRID_K[10]=3; GRID_P[10]=20;  GRID_SHIFT[10]=1.0

if [[ -z "${GRID_TAG[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

TAG="${GRID_TAG[${TASK_ID}]}"
N_TOTAL="${GRID_N[${TASK_ID}]}"
K="${GRID_K[${TASK_ID}]}"
P="${GRID_P[${TASK_ID}]}"
SHIFT="${GRID_SHIFT[${TASK_ID}]}"

echo "======================================================"
echo " FACE-HD Bias Diagnosis | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:        ${TAG}"
echo " n_total:    ${N_TOTAL}"
echo " K:          ${K}"
echo " p:          ${P}"
echo " shift:      ${SHIFT}"
echo " n_sims:     ${N_SIMS}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  (SLURM_JOB_ID=${SLURM_JOB_ID})"
echo " Checkpoint: ${CHECKPOINT_DIR}"
echo "======================================================"

# Run R in BACKGROUND so the USR1 trap is serviced (mirrors main.cmd).
# main.cmd uses R CMD BATCH; we use Rscript so stdout/stderr stream to the
# SLURM log directly without an extra .Rout file.
Rscript --vanilla diagnosis/bias/run_setting.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SHIFT}" "${N_SIMS}" &
R_PID=$!
wait "${R_PID}"
R_EXIT=$?

# If we got here without preemption, clean up any stray signal files.
rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}" \
      "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}"

echo
echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date) (exit=${R_EXIT})"
echo "======================================================"

exit "${R_EXIT}"
