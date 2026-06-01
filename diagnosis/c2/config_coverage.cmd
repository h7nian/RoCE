#!/bin/bash
#SBATCH --job-name=FACE-cfg-cov
#SBATCH --output=diagnosis/c2/config_coverage/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/config_coverage/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=8

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_CFG_OUTPUT_ROOT="${C2_CFG_OUTPUT_ROOT:-diagnosis/c2/config_coverage}"
mkdir -p "${C2_CFG_OUTPUT_ROOT}/log" "${C2_CFG_OUTPUT_ROOT}/results"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

# Optional isolated lib holding the patched FACEHD (fix validation).
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
N_FOLDS="${C2_CFG_FOLDS:-10}"
N_CORES="${C2_CFG_CORES:-1}"
MODES="${C2_CFG_MODES:-one_round}"
MODES="${MODES//:/,}"
K="${C2_CFG_K:-3}"
P="${C2_CFG_P:-10}"
SEEDS_PER_CELL="${C2_CFG_SEEDS_PER_CELL:-40}"
SEED_BASE="${C2_CFG_SEED_BASE:-3000}"

# Cells = config x n. Regression check: do C1/C3 (currently good) stay ~0.95 and
# C2 recover, under the patched density-ratio CV scale?
CELLS=("C1 5000" "C1 20000" "C2 5000" "C2 20000" "C3 5000" "C3 20000" "C4 5000" "C4 20000")

cell_idx=$(( (TASK_ID - 1) / SEEDS_PER_CELL ))
seed_in_cell=$(( (TASK_ID - 1) % SEEDS_PER_CELL + 1 ))
if (( cell_idx < 0 || cell_idx >= ${#CELLS[@]} )); then
    echo "ERROR: TASK_ID=${TASK_ID} -> out-of-range cell ${cell_idx}" >&2
    exit 64
fi
read -r CONFIG N_TOTAL <<< "${CELLS[${cell_idx}]}"
SEED=$(( SEED_BASE + seed_in_cell ))
TAG="cfgcov"

echo "=== FACE-HD config coverage | task ${TASK_ID} | ${CONFIG} n=${N_TOTAL} K=${K} p=${P} seed=${SEED} ==="
echo " lib: ${R_LIBS_USER}   out: ${C2_CFG_OUTPUT_ROOT}   start: $(date)"

Rscript --vanilla diagnosis/c2/config_coverage.R \
    "${TAG}" "${CONFIG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" \
    "${N_FOLDS}" "${N_CORES}" "${MODES}"

echo "=== task ${TASK_ID} (${CONFIG} n=${N_TOTAL} seed=${SEED}) done $(date) ==="
