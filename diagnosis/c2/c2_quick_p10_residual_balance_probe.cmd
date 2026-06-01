#!/bin/bash
#SBATCH --job-name=FACE-c2-p10quick
#SBATCH --output=diagnosis/c2/residual_balance_p10_quick_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/residual_balance_p10_quick_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_RBAL_OUTPUT_ROOT="${C2_RBAL_OUTPUT_ROOT:-diagnosis/c2/residual_balance_p10_quick_probe}"
mkdir -p "${C2_RBAL_OUTPUT_ROOT}/log" "${C2_RBAL_OUTPUT_ROOT}/results"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
N_FOLDS="${C2_RBAL_FOLDS:-10}"
N_CORES="${C2_RBAL_CORES:-1}"
NLAMBDA_INIT="${C2_RBAL_NLAMBDA_INIT:-100}"
MODES="${C2_RBAL_MODES:-one_round,two_round}"
MODES="${MODES//:/,}"

# Fast p=10 C2 checks.  Include known problematic seeds plus nearby seeds so
# we can see whether the center bias is seed-specific or systematic.
# Grid columns:
#   tag n_total K p seed
declare -A GRID
GRID[1]="c2quick_k3p10_seed40 5000 3 10 40"
GRID[2]="c2quick_k3p10_seed120 5000 3 10 120"
GRID[3]="c2quick_k3p10_seed220 5000 3 10 220"
GRID[4]="c2quick_k3p10_seed340 5000 3 10 340"
GRID[5]="c2quick_k3p10_seed460 5000 3 10 460"
GRID[6]="c2quick_k4p10_seed364 5000 4 10 364"
GRID[7]="c2quick_k4p10_seed464 5000 4 10 464"
GRID[8]="c2quick_k4p10_seed564 5000 4 10 564"

if [[ -z "${GRID[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no quick p=10 grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

read -r TAG N_TOTAL K P SEED <<< "${GRID[${TASK_ID}]}"

echo "======================================================"
echo " FACE-HD C2 p=10 Quick Residual Balance | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:        ${TAG}"
echo " n_total:    ${N_TOTAL}"
echo " K:          ${K}"
echo " p:          ${P}"
echo " seed:       ${SEED}"
echo " K_f:        ${N_FOLDS}"
echo " nlambda:    ${NLAMBDA_INIT}"
echo " n_cores:    ${N_CORES}"
echo " modes:      ${MODES}"
echo " output root:${C2_RBAL_OUTPUT_ROOT}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_residual_balance_probe.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" \
    "${N_FOLDS}" "${NLAMBDA_INIT}" "${N_CORES}" "${MODES}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
