#!/bin/bash
#SBATCH --job-name=FACE-c2lamgrid
#SBATCH --output=diagnosis/c2/lambda_grid_compare/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/lambda_grid_compare/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_LAMBDA_GRID_OUTPUT_ROOT="diagnosis/c2/lambda_grid_compare"
mkdir -p "${C2_LAMBDA_GRID_OUTPUT_ROOT}/log" "${C2_LAMBDA_GRID_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
# Defaults = full C2 p50 cell; smoke overrides via --export.
N_TOTAL="${C2LG_N_TOTAL:-5000}"; N_FOLDS="${C2LG_N_FOLDS:-10}"; P="${C2LG_P:-50}"; K="${C2LG_K:-3}"
NLAMBDA="${C2LG_NLAMBDA:-100}"; SEED_BASE="${C2LG_SEED_BASE:-60000}"; TAGPFX="${C2LG_TAGPFX:-c2lamgrid}"
SEED=$(( SEED_BASE + TASK_ID )); TAG="${TAGPFX}_seed${SEED}"
echo "=== C2 lambda-grid | task ${TASK_ID} | n=${N_TOTAL} K=${K} p=${P} kf=${N_FOLDS} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_lambda_grid_compare.R "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" 1 one_round
echo "=== done $(date) ==="
