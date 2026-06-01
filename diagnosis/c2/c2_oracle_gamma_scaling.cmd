#!/bin/bash
#SBATCH --job-name=FACE-c2-orcl-scale
#SBATCH --output=diagnosis/c2/oracle_gamma_scaling/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/oracle_gamma_scaling/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=8

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_ORACLE_GAMMA_OUTPUT_ROOT="${C2_ORACLE_GAMMA_OUTPUT_ROOT:-diagnosis/c2/oracle_gamma_scaling}"
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

# Oracle-gamma coverage vs n: does the ORACLE (true density-ratio) estimator keep
# ~0.95 coverage at large n while the FITTED-gamma estimator degrades (0.94->0.84)?
# If yes, identification/estimand is sound and the bug is gamma-hat ESTIMATION
# (the density-ratio CV validation-loss scale in cv_utils.hpp). Reuses the proven,
# UNMODIFIED diagnosis/c2/c2_oracle_gamma_probe.R.
#
# Grid computed from TASK_ID: cells (K,n) x SEEDS_PER_CELL seeds.
SEEDS_PER_CELL="${C2_ORC_SEEDS_PER_CELL:-60}"
SEED_BASE="${C2_ORC_SEED_BASE:-2000}"
CELLS=("3 5000" "3 20000" "4 5000" "4 20000")
P=10

cell_idx=$(( (TASK_ID - 1) / SEEDS_PER_CELL ))
seed_in_cell=$(( (TASK_ID - 1) % SEEDS_PER_CELL + 1 ))
if (( cell_idx < 0 || cell_idx >= ${#CELLS[@]} )); then
    echo "ERROR: TASK_ID=${TASK_ID} maps to out-of-range cell ${cell_idx}" >&2
    exit 64
fi
read -r K N_TOTAL <<< "${CELLS[${cell_idx}]}"
SEED=$(( SEED_BASE + seed_in_cell ))
TAG="c2orc_k${K}_n${N_TOTAL}_seed${SEED}"

echo "======================================================"
echo " FACE-HD C2 Oracle-Gamma Scaling | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:        ${TAG}  (K=${K}, n=${N_TOTAL}, seed=${SEED})"
echo " output root:${C2_ORACLE_GAMMA_OUTPUT_ROOT}"
echo " Start:      $(date)"
echo " Job:        ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_oracle_gamma_probe.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
