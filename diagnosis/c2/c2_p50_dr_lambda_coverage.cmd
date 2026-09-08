#!/bin/bash
#SBATCH --job-name=FACE-c2p50lamcov
#SBATCH --output=diagnosis/c2/dr_lambda_coverage_p50/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/dr_lambda_coverage_p50/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
LIB="${C2_HOOK_LIB:-${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}/diagnosis/c2/lambda_scale_lib}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="diagnosis/c2/dr_lambda_coverage_p50"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
# load the HOOKED RoCE (with ROCE_DR_LAMBDA_SCALE) from the isolated lib
export R_LIBS_USER="${LIB}:${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
N_FOLDS=10; NLAMBDA=100; NCORES=1; MODES="one_round"
K=3; P=50; N_TOTAL=5000
SEEDS_PER_CELL="${C2DLC_SEEDS_PER_CELL:-20}"
SEED_BASE="${C2DLC_SEED_BASE:-40000}"

# Cells = gamma-lambda scale factor. Tests whether LOOSENING (x0.1) or TIGHTENING (x10) the
# density-ratio lambda relative to the CV-selected value improves FULL-estimator COVERAGE
# (not just moment bias) at C2 p=50. Uses the hooked RoCE via ROCE_DR_LAMBDA_SCALE.
SCALES=(0.1 1.0 10)
cell=$(( (TASK_ID-1)/SEEDS_PER_CELL ))
sic=$(( (TASK_ID-1)%SEEDS_PER_CELL + 1 ))
if (( cell<0 || cell>=${#SCALES[@]} )); then echo "ERR cell $cell" >&2; exit 64; fi
export ROCE_DR_LAMBDA_SCALE="${SCALES[$cell]}"
SEED=$(( SEED_BASE + sic ))
TAG="c2p50dlc_x${ROCE_DR_LAMBDA_SCALE}_seed${SEED}"

echo "=== C2 p50 dr-lambda coverage | task ${TASK_ID} | scale=${ROCE_DR_LAMBDA_SCALE} n=${N_TOTAL} K=${K} p=${P} seed=${SEED} | lib=${LIB} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_score_moment_audit.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" "${NCORES}" "${MODES}"
echo "=== task ${TASK_ID} done $(date) ==="
