#!/bin/bash
#SBATCH --job-name=FACE-c2-scale
#SBATCH --output=diagnosis/c2/score_moment_scaling/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/score_moment_scaling/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="${C2_SCORE_AUDIT_OUTPUT_ROOT:-diagnosis/c2/score_moment_scaling}"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
N_FOLDS="${C2_SCORE_AUDIT_FOLDS:-10}"
N_CORES="${C2_SCORE_AUDIT_CORES:-1}"
NLAMBDA_INIT="${C2_SCORE_AUDIT_NLAMBDA_INIT:-100}"
MODES="${C2_SCORE_AUDIT_MODES:-one_round,two_round}"
MODES="${MODES//:/,}"
USE_LAMBDA_CACHE="${C2_SCORE_AUDIT_USE_LAMBDA_CACHE:-true}"

# n-scaling sweep of the (proven) C2 score/moment audit on the failing p=10 settings.
# Goal: test whether out_noise = E_s[w (Y - m_true)] and |bias|/se SHRINK with n
# (finite-sample higher-order bias -> coverage recovers) or PERSIST (structural).
# Grid columns:
#   tag n_total K p seed
declare -A GRID
# --- K=3, p=10 ---
GRID[1]="c2scale_k3p10_seed120_n5000  5000  3 10 120"
GRID[2]="c2scale_k3p10_seed120_n10000 10000 3 10 120"
GRID[3]="c2scale_k3p10_seed120_n20000 20000 3 10 120"
GRID[4]="c2scale_k3p10_seed220_n5000  5000  3 10 220"
GRID[5]="c2scale_k3p10_seed220_n10000 10000 3 10 220"
GRID[6]="c2scale_k3p10_seed220_n20000 20000 3 10 220"
# --- K=4, p=10 ---
GRID[7]="c2scale_k4p10_seed364_n5000  5000  4 10 364"
GRID[8]="c2scale_k4p10_seed364_n10000 10000 4 10 364"
GRID[9]="c2scale_k4p10_seed364_n20000 20000 4 10 364"
GRID[10]="c2scale_k4p10_seed464_n5000  5000  4 10 464"
GRID[11]="c2scale_k4p10_seed464_n10000 10000 4 10 464"
GRID[12]="c2scale_k4p10_seed464_n20000 20000 4 10 464"

if [[ -z "${GRID[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no scaling grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

read -r TAG N_TOTAL K P SEED <<< "${GRID[${TASK_ID}]}"

echo "======================================================"
echo " RoCE C2 Score/Moment n-SCALING | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:         ${TAG}"
echo " n_total:     ${N_TOTAL}"
echo " K:           ${K}"
echo " p:           ${P}"
echo " seed:        ${SEED}"
echo " K_f:         ${N_FOLDS}"
echo " nlambda:     ${NLAMBDA_INIT}"
echo " n_cores:     ${N_CORES}"
echo " modes:       ${MODES}"
echo " lambda cache:${USE_LAMBDA_CACHE}"
echo " output root: ${C2_SCORE_AUDIT_OUTPUT_ROOT}"
echo " Start:       $(date)"
echo " Host:        $(hostname)"
echo " Job:         ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_score_moment_audit.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" \
    "${N_FOLDS}" "${NLAMBDA_INIT}" "${N_CORES}" "${MODES}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
