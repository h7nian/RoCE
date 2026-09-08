#!/bin/bash
#SBATCH --job-name=FACE-c2p50ns
#SBATCH --output=diagnosis/c2/score_moment_p50_nscaling/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/score_moment_p50_nscaling/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="diagnosis/c2/score_moment_p50_nscaling"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID required}"
N_FOLDS=10; NLAMBDA=100; NCORES=1; MODES="one_round"
K=3; P=50
SEEDS_PER_CELL="${C2P50_SEEDS_PER_CELL:-20}"
SEED_BASE="${C2P50_SEED_BASE:-20000}"

# Evidence that the C2 p=50 under-coverage is NUISANCE ESTIMATION error: if so, the
# RoCE bias should SHRINK as n grows (nuisances converge) at fixed p=50.
# Reuses the proven, UNMODIFIED diagnosis/c2/c2_score_moment_audit.R (config=C2 inside).
NS=(5000 10000 20000)
cell=$(( (TASK_ID-1)/SEEDS_PER_CELL ))
sic=$(( (TASK_ID-1)%SEEDS_PER_CELL + 1 ))
if (( cell<0 || cell>=${#NS[@]} )); then echo "ERR cell $cell" >&2; exit 64; fi
N_TOTAL=${NS[$cell]}
SEED=$(( SEED_BASE + sic ))
TAG="c2p50_n${N_TOTAL}_seed${SEED}"

echo "=== C2 p50 n-scaling | task ${TASK_ID} | n=${N_TOTAL} K=${K} p=${P} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_score_moment_audit.R \
    "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" "${NCORES}" "${MODES}"
echo "=== task ${TASK_ID} done $(date) ==="
