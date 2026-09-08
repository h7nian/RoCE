#!/bin/bash
#SBATCH --job-name=FACE-c2-cov
#SBATCH --output=diagnosis/c2/coverage_scaling/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/coverage_scaling/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=8

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="${C2_SCORE_AUDIT_OUTPUT_ROOT:-diagnosis/c2/coverage_scaling}"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

# Optional: prepend an isolated library holding a rebuilt RoCE (fix validation),
# so library(RoCE) loads the patched build while dependencies resolve from ~/Rlibs.
if [[ -n "${C2_FIX_LIB:-}" ]]; then
    export R_LIBS_USER="${C2_FIX_LIB}:${HOME}/Rlibs"
else
    export R_LIBS_USER="${HOME}/Rlibs"
fi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
N_FOLDS="${C2_SCORE_AUDIT_FOLDS:-10}"
N_CORES="${C2_SCORE_AUDIT_CORES:-1}"
NLAMBDA_INIT="${C2_SCORE_AUDIT_NLAMBDA_INIT:-100}"
MODES="${C2_SCORE_AUDIT_MODES:-one_round}"
MODES="${MODES//:/,}"
USE_LAMBDA_CACHE="${C2_SCORE_AUDIT_USE_LAMBDA_CACHE:-true}"

# Coverage study: many seeds per cell to estimate actual CI coverage at small vs
# large n, confirming whether the C2 p=10 coverage dip is finite-sample (recovers
# with n). Reuses the proven, UNMODIFIED diagnosis/c2/c2_score_moment_audit.R.
#
# Grid is computed from TASK_ID (no giant hardcoded table):
#   cells (K, n): {(3,5000),(3,20000),(4,5000),(4,20000)}
#   SEEDS_PER_CELL seeds per cell, starting at SEED_BASE+1.
SEEDS_PER_CELL="${C2_COV_SEEDS_PER_CELL:-50}"
SEED_BASE="${C2_COV_SEED_BASE:-1000}"
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
TAG="c2cov_k${K}_n${N_TOTAL}_seed${SEED}"

echo "======================================================"
echo " RoCE C2 Coverage Scaling | Task ${TASK_ID}"
echo "======================================================"
echo " Tag:         ${TAG}"
echo " cell:        ${cell_idx} (K=${K}, n=${N_TOTAL})"
echo " seed:        ${SEED}"
echo " p:           ${P}"
echo " K_f:         ${N_FOLDS}"
echo " modes:       ${MODES}"
echo " output root: ${C2_SCORE_AUDIT_OUTPUT_ROOT}"
echo " Start:       $(date)"
echo " Job:         ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================"

Rscript --vanilla diagnosis/c2/c2_score_moment_audit.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" \
    "${N_FOLDS}" "${NLAMBDA_INIT}" "${N_CORES}" "${MODES}"

echo "======================================================"
echo " Task ${TASK_ID} (${TAG}) completed at $(date)"
echo "======================================================"
