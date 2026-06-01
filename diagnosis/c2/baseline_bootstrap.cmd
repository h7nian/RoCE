#!/bin/bash
#SBATCH --job-name=FACE-base-boot
#SBATCH --output=diagnosis/c2/baseline_bootstrap/log/%x_%A_%a.out
#SBATCH --error=diagnosis/c2/baseline_bootstrap/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=8

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export C2_BOOT_OUTPUT_ROOT="${C2_BOOT_OUTPUT_ROOT:-diagnosis/c2/baseline_bootstrap}"
mkdir -p "${C2_BOOT_OUTPUT_ROOT}/log" "${C2_BOOT_OUTPUT_ROOT}/results"
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
N_FOLDS="${C2_BOOT_FOLDS:-10}"
K="${C2_BOOT_K:-3}"; P="${C2_BOOT_P:-10}"
SEEDS_PER_CELL="${C2_BOOT_SEEDS_PER_CELL:-50}"
SEED_BASE="${C2_BOOT_SEED_BASE:-9000}"
export C2_BOOT_B="${C2_BOOT_B:-1000}"

CELLS=("C1 5000" "C2 5000" "C3 5000" "C4 5000" "C1 20000" "C2 20000" "C3 20000" "C4 20000")
cell_idx=$(( (TASK_ID - 1) / SEEDS_PER_CELL ))
seed_in_cell=$(( (TASK_ID - 1) % SEEDS_PER_CELL + 1 ))
if (( cell_idx < 0 || cell_idx >= ${#CELLS[@]} )); then echo "ERR cell $cell_idx" >&2; exit 64; fi
read -r CONFIG N_TOTAL <<< "${CELLS[${cell_idx}]}"
SEED=$(( SEED_BASE + seed_in_cell ))

echo "=== baseline variance (analytic vs exact-FE vs bootstrap) | task ${TASK_ID} | ${CONFIG} n=${N_TOTAL} seed=${SEED} | $(date) ==="
Rscript --vanilla diagnosis/c2/baseline_bootstrap.R "basevar" "${CONFIG}" "${N_TOTAL}" "${K}" "${P}" "${SEED}" "${N_FOLDS}"
echo "=== task ${TASK_ID} done $(date) ==="
