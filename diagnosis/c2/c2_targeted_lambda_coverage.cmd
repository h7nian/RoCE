#!/bin/bash
#SBATCH --job-name=FACE-c2tgtcov
#SBATCH --output=diagnosis/c2/targeted_lambda_coverage/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/targeted_lambda_coverage/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2TLC_OUTPUT_ROOT="diagnosis/c2/targeted_lambda_coverage"
mkdir -p "${C2TLC_OUTPUT_ROOT}/log" "${C2TLC_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
N_TOTAL="${C2TLC_N_TOTAL:-5000}"; N_FOLDS="${C2TLC_N_FOLDS:-10}"; P="${C2TLC_P:-50}"; K="${C2TLC_K:-3}"
NLAMBDA="${C2TLC_NLAMBDA:-100}"; SEED_BASE="${C2TLC_SEED_BASE:-90000}"; TAGPFX="${C2TLC_TAGPFX:-c2tgtcov}"
SEED=$(( SEED_BASE + TASK_ID )); TAG="${TAGPFX}_seed${SEED}"
echo "=== C2 targeted-lambda COVERAGE | task ${TASK_ID} | n=${N_TOTAL} K=${K} p=${P} kf=${N_FOLDS} seed=${SEED} | stride=${C2TLC_BAL_STRIDE:-4} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_targeted_lambda_coverage.R "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" 1 one_round
echo "=== done $(date) ==="
