#!/bin/bash
#SBATCH --job-name=roce_sensitivity_audit
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_sensitivity_audit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_sensitivity_audit.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
EXPECTED_REPLICATIONS="${2:-500}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}/logs"
Rscript scripts/slurm/aggregate_reused_sensitivities.R \
  "${RESULT_ROOT}" "${EXPECTED_REPLICATIONS}"
