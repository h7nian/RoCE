#!/bin/bash
#SBATCH --job-name=FACE-c2p50lam
#SBATCH --output=diagnosis/c2/residual_balance_p50_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/residual_balance_p50_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_RBAL_OUTPUT_ROOT="diagnosis/c2/residual_balance_p50_probe"
mkdir -p "${C2_RBAL_OUTPUT_ROOT}/log" "${C2_RBAL_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
N_FOLDS=10; NLAMBDA=100; NCORES=1; MODES="one_round"
K=3; P=50; N_TOTAL=5000
SEED_BASE="${C2P50L_SEED_BASE:-30000}"
SEED=$(( SEED_BASE + TASK_ID ))
TAG="c2p50lam_seed${SEED}"

# "Can the p=50 C2 bias be improved by better nuisance regularization?" The residual
# balance probe re-fits gamma/alpha under lambda variants {current, x0.1, x10, 1se, zero, true}
# on the SAME C2 folds and reports the resulting source-correction bias. If x10/1se/true-gamma
# materially shrink the bias, gamma-hat is improvably regularized (SMMAL truncation/calibration lever);
# if only true-gamma helps, the residual is the irreducible high-dim estimation rate.
echo "=== C2 p50 lambda-sensitivity | task ${TASK_ID} | n=${N_TOTAL} K=${K} p=${P} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_residual_balance_probe.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" "${NCORES}" "${MODES}"
echo "=== task ${TASK_ID} done $(date) ==="
