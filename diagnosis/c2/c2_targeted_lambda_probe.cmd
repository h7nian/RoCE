#!/bin/bash
#SBATCH --job-name=FACE-c2tgtlam
#SBATCH --output=diagnosis/c2/targeted_lambda_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/targeted_lambda_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_TARGETED_LAMBDA_OUTPUT_ROOT="diagnosis/c2/targeted_lambda_probe"
mkdir -p "${C2_TARGETED_LAMBDA_OUTPUT_ROOT}/log" "${C2_TARGETED_LAMBDA_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
N_TOTAL="${C2TL_N_TOTAL:-5000}"; N_FOLDS="${C2TL_N_FOLDS:-10}"; P="${C2TL_P:-50}"; K="${C2TL_K:-3}"
NLAMBDA="${C2TL_NLAMBDA:-100}"; SEED_BASE="${C2TL_SEED_BASE:-80000}"; TAGPFX="${C2TL_TAGPFX:-c2tgtlam}"
SEED=$(( SEED_BASE + TASK_ID )); TAG="${TAGPFX}_seed${SEED}"
echo "=== C2 targeted-lambda | task ${TASK_ID} | n=${N_TOTAL} K=${K} p=${P} kf=${N_FOLDS} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_targeted_lambda_probe.R "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" 1 one_round
echo "=== done $(date) ==="
