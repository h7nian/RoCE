#!/bin/bash
#SBATCH --job-name=FACE-c2-nrefit
#SBATCH --output=diagnosis/c2/nuisance_refit_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/nuisance_refit_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/c2/nuisance_refit_probe/log diagnosis/c2/nuisance_refit_probe/results

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
N_FOLDS="${C2_NREFIT_FOLDS:-10}"
N_CORES="${C2_NREFIT_CORES:-1}"
NLAMBDA_INIT="${C2_NREFIT_NLAMBDA_INIT:-100}"
MODES="${C2_NREFIT_MODES:-one_round,two_round}"
MODES="${MODES//:/,}"

# Grid columns:
#   tag n_total K p seed
declare -A GRID
GRID[1]="c2_k3p10_seed40 5000 3 10 40"
GRID[2]="c2_k3p10_seed120 5000 3 10 120"
GRID[3]="c2_k3p50_seed9 5000 3 50 9"
GRID[4]="c2_k4p10_seed364 5000 4 10 364"
GRID[5]="c2_k4p50_seed484 5000 4 50 484"

if [[ -z "${GRID[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

read -r TAG N_TOTAL K P SEED <<< "${GRID[${TASK_ID}]}"

echo "======================================================"
echo " RoCE C2 Nuisance Refit Probe | Task ${TASK_ID}"
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
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_nuisance_refit_probe.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" \
    "${N_FOLDS}" "${NLAMBDA_INIT}" "${N_CORES}" "${MODES}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
