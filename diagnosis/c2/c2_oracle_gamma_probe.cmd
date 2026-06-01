#!/bin/bash
#SBATCH --job-name=FACE-c2-oracle-gamma
#SBATCH --output=diagnosis/c2/oracle_gamma_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/oracle_gamma_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_ORACLE_GAMMA_OUTPUT_ROOT="${C2_ORACLE_GAMMA_OUTPUT_ROOT:-diagnosis/c2/oracle_gamma_probe}"
mkdir -p "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/log" "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/results"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"

declare -A GRID
GRID[1]="c2oracle_k3p10_seed40 5000 3 10 40"
GRID[2]="c2oracle_k3p10_seed120 5000 3 10 120"
GRID[3]="c2oracle_k3p10_seed220 5000 3 10 220"
GRID[4]="c2oracle_k3p10_seed340 5000 3 10 340"
GRID[5]="c2oracle_k3p10_seed460 5000 3 10 460"
GRID[6]="c2oracle_k4p10_seed364 5000 4 10 364"
GRID[7]="c2oracle_k4p10_seed464 5000 4 10 464"
GRID[8]="c2oracle_k4p10_seed564 5000 4 10 564"

if [[ -z "${GRID[${TASK_ID}]:-}" ]]; then
    echo "ERROR: no oracle gamma grid entry for SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
    exit 64
fi

read -r TAG N_TOTAL K P SEED <<< "${GRID[${TASK_ID}]}"

echo "======================================================"
echo " FACE-HD C2 Oracle Gamma Probe | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:        ${TAG}"
echo " n_total:    ${N_TOTAL}"
echo " K:          ${K}"
echo " p:          ${P}"
echo " seed:       ${SEED}"
echo " output root:${C2_ORACLE_GAMMA_OUTPUT_ROOT}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_oracle_gamma_probe.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
