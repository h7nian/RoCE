#!/bin/bash
#SBATCH --job-name=roce_rho0_truth
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rho0_truth.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rho0_truth.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${ROCE_RESULT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
OUTPUT_DIRECTORY="${ROCE_OUTPUT_ROOT:-${RESULT_ROOT}/rho0_site_tate}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"
N_REFERENCE="${ROCE_N_REFERENCE:-500000}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}/logs"
Rscript scripts/slurm/diagnose_rho0_site_tate.R \
  "${OUTPUT_DIRECTORY}" "${N_REFERENCE}"
