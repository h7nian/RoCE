#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=diagnosis/bias/log/%x_%j.out
#SBATCH --error=diagnosis/bias/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5
#SBATCH --requeue
#SBATCH --signal=B:USR1@120

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/bias/log diagnosis/bias/verify_out diagnosis/bias/checkpoints

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NLAMBDA_INIT=100

N_FOLDS="${VERIFY_N_FOLDS:-3}"
N_SIMS="${VERIFY_N_SIMS:-200}"
CONFIG="${VERIFY_CONFIG:-C1}"

# Checkpoint dir keyed to (verify, config, K_f). Stable across requeues so the
# in-R load_checkpoint() picks up exactly where the last attempt stopped.
export CHECKPOINT_DIR="diagnosis/bias/checkpoints/verify_${CONFIG}_kf${N_FOLDS}"
export CHECKPOINT_JOB_ID="${SLURM_JOB_ID}"
mkdir -p "${CHECKPOINT_DIR}"

# Preemption handler mirrors main.cmd: touch preempt file, wait up to 100s
# for R to flush the save flag, then scontrol requeue.
handle_preemption() {
    echo "$(date): Received SIGUSR1. Asking R to checkpoint..."
    touch "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}"
    for _ in $(seq 1 100); do
        if [[ -f "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}" ]]; then
            echo "$(date): Checkpoint saved."
            rm -f "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}"
            break
        fi
        sleep 1
    done
    rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}"
    echo "$(date): Requesting requeue (JOB_ID=${SLURM_JOB_ID})..."
    scontrol requeue "${SLURM_JOB_ID}" || true
    exit 0
}
trap 'handle_preemption' USR1

echo "======================================================"
echo " FACE-HD K_f Verification Run"
echo "======================================================"
echo " Config:     ${CONFIG}"
echo " K_f:        ${N_FOLDS}"
echo " n_sims:     ${N_SIMS}"
echo " Job:        ${SLURM_JOB_ID}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo " Checkpoint: ${CHECKPOINT_DIR}"
echo "======================================================"

Rscript --vanilla diagnosis/bias/verify_kf3.R "${N_FOLDS}" "${N_SIMS}" "${CONFIG}" &
R_PID=$!
wait "${R_PID}"
R_EXIT=$?

rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}" \
      "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}"

echo
echo "======================================================"
echo " verify_kf done at $(date) (exit=${R_EXIT})"
echo "======================================================"
exit "${R_EXIT}"
