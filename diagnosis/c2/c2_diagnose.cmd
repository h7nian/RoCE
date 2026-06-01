#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --job-name=FACE-c2-diag
#SBATCH --output=diagnosis/c2/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5
#SBATCH --requeue
#SBATCH --signal=B:USR1@120
# Other --time/--mem/--cpus-per-task/--partition/--array come from
# c2_diagnose.sh via the sbatch CLI.

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/c2/log diagnosis/c2/results diagnosis/c2/checkpoints diagnosis/c2/tmp

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1
export NLAMBDA_INIT=100

N_SIMS="${C2_DIAG_N_SIMS:-100}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

export CHECKPOINT_DIR="diagnosis/c2/checkpoints/task_${TASK_ID}"
export CHECKPOINT_JOB_ID="${SLURM_JOB_ID}"
mkdir -p "${CHECKPOINT_DIR}"

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

# Grid columns:
#   tag n_total K p config estimand allocation transform heterogeneity shift n_folds
declare -A GRID
GRID[1]="c2_k3p10_super_kf10 5000 3 10 C2 superpopulation model mild none 0.5 10"
GRID[2]="c2_k3p10_sample_kf10 5000 3 10 C2 sample model mild none 0.5 10"
GRID[3]="c2_k3p50_super_kf10 5000 3 50 C2 superpopulation model mild none 0.5 10"
GRID[4]="c2_k3p50_sample_kf10 5000 3 50 C2 sample model mild none 0.5 10"

# These match the C2 rows in the result table where coverage was low.  Tasks
# are paired by superpopulation/sample estimand to separate Monte Carlo target
# variation from nuisance/aggregation bias.
GRID[5]="c2_k4p10_super_kf10 5000 4 10 C2 superpopulation model mild none 0.5 10"
GRID[6]="c2_k4p10_sample_kf10 5000 4 10 C2 sample model mild none 0.5 10"
GRID[7]="c2_k4p50_super_kf10 5000 4 50 C2 superpopulation model mild none 0.5 10"
GRID[8]="c2_k4p50_sample_kf10 5000 4 50 C2 sample model mild none 0.5 10"
GRID[9]="c2_k5p10_super_kf10 5000 5 10 C2 superpopulation model mild none 0.5 10"
GRID[10]="c2_k5p10_sample_kf10 5000 5 10 C2 sample model mild none 0.5 10"
GRID[11]="c2_k5p50_super_kf10 5000 5 50 C2 superpopulation model mild none 0.5 10"
GRID[12]="c2_k5p50_sample_kf10 5000 5 50 C2 sample model mild none 0.5 10"

# Larger-n check.  If the C2 residual bias is a second-order nuisance term, the
# bias/se ratio should improve materially when n doubles at fixed p,K.
GRID[13]="c2_n10000_k3p10_super_kf10 10000 3 10 C2 superpopulation model mild none 0.5 10"
GRID[14]="c2_n10000_k3p10_sample_kf10 10000 3 10 C2 sample model mild none 0.5 10"

# Hard-design probes for paper positioning, all at the default K_f=10.  These
# stress source borrowing under valid fitted bases (C1) before reintroducing the
# C2 outcome-basis issue.  Partial heterogeneity leaves some transportable
# sources; strong heterogeneity checks whether adaptive aggregation backs away.
GRID[15]="hard_independent_partial_C1_kf10 5000 5 10 C1 superpopulation independent mild partial 1.0 10"
GRID[16]="hard_sourceheavy_partial_C1_kf10 5000 5 10 C1 superpopulation source_heavy mild partial 0.5 10"
GRID[17]="hard_independent_strong_C1_kf10 5000 5 10 C1 superpopulation independent mild strong 1.0 10"
GRID[18]="hard_sourceheavy_strong_C1_kf10 5000 5 10 C1 superpopulation source_heavy mild strong 0.5 10"
GRID[19]="hard_independent_partial_C2_kf10 5000 5 10 C2 superpopulation independent mild partial 1.0 10"
GRID[20]="hard_sourceheavy_partial_C2_kf10 5000 5 10 C2 superpopulation source_heavy mild partial 0.5 10"

if [[ -z "${GRID[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

read -r TAG N_TOTAL K P CONFIG ESTIMAND ALLOCATION TRANSFORM HET SHIFT N_FOLDS <<< "${GRID[${TASK_ID}]}"

echo "======================================================"
echo " FACE-HD C2 Diagnosis | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:        ${TAG}"
echo " n_total:    ${N_TOTAL}"
echo " K:          ${K}"
echo " p:          ${P}"
echo " config:     ${CONFIG}"
echo " estimand:   ${ESTIMAND}"
echo " allocation: ${ALLOCATION}"
echo " transform:  ${TRANSFORM}"
echo " hetero:     ${HET}"
echo " shift:      ${SHIFT}"
echo " K_f:        ${N_FOLDS}"
echo " n_sims:     ${N_SIMS}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}  (SLURM_JOB_ID=${SLURM_JOB_ID})"
echo " Checkpoint: ${CHECKPOINT_DIR}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/run_setting.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${CONFIG}" "${ESTIMAND}" \
    "${ALLOCATION}" "${TRANSFORM}" "${HET}" "${SHIFT}" "${N_FOLDS}" "${N_SIMS}" &
R_PID=$!
wait "${R_PID}"
R_EXIT=$?

rm -f "${CHECKPOINT_DIR}/.preempt_signal_${SLURM_JOB_ID}" \
      "${CHECKPOINT_DIR}/.checkpoint_saved_${SLURM_JOB_ID}"

echo
echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date) (exit=${R_EXIT})"
echo "======================================================"
exit "${R_EXIT}"
