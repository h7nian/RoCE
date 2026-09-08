#!/bin/bash
#SBATCH --job-name=roce_nlambda_compare
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_nlambda_compare.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_nlambda_compare.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${ROCE_RESULT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
GRID_A_RESULT="${1:-${RESULT_ROOT}/smoke_nlambda50_raw/task_000001.csv}"
GRID_B_RESULT="${2:-${RESULT_ROOT}/smoke_single_raw/task_000001.csv}"
OUTPUT_PATH="${3:-${RESULT_ROOT}/nlambda50_vs_100.csv}"

for result_path in "${GRID_A_RESULT}" "${GRID_B_RESULT}"; do
  if [[ ! -f "${result_path}" ]]; then
    echo "smoke result not found: ${result_path}" >&2
    exit 1
  fi
done

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "$(dirname "${OUTPUT_PATH}")" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/compare_nlambda_smokes.R \
  "${GRID_A_RESULT}" "${GRID_B_RESULT}" "${OUTPUT_PATH}"
