#!/bin/bash
#SBATCH --job-name=FACE-c2p50mtau
#SBATCH --output=diagnosis/c2/mtau_p50/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/mtau_p50/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_SCORE_AUDIT_OUTPUT_ROOT="diagnosis/c2/mtau_p50"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
TASK_ID="${SLURM_ARRAY_TASK_ID:?}"
N_FOLDS=10; NLAMBDA=100; NCORES=1; MODES="one_round"; K=3; P=50; N_TOTAL=5000
SEEDS_PER_CELL="${C2MT_SEEDS_PER_CELL:-20}"; SEED_BASE="${C2MT_SEED_BASE:-50000}"
MTAUS=(3 4.4 10)
cell=$(( (TASK_ID-1)/SEEDS_PER_CELL )); sic=$(( (TASK_ID-1)%SEEDS_PER_CELL + 1 ))
if (( cell<0 || cell>=${#MTAUS[@]} )); then echo "ERR cell $cell">&2; exit 64; fi
export C2_M_TAU="${MTAUS[$cell]}"
SEED=$(( SEED_BASE + sic )); TAG="c2p50mtau_M${C2_M_TAU}_seed${SEED}"
echo "=== C2 p50 M_tau | task ${TASK_ID} | M_tau=${C2_M_TAU} n=${N_TOTAL} K=${K} p=${P} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/c2_score_moment_audit.R "${TAG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}" "${NLAMBDA}" "${NCORES}" "${MODES}"
echo "=== done $(date) ==="
