#!/bin/bash
#SBATCH --job-name=roce_nlambda_paired
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_nlambda_paired.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_nlambda_paired.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${PROJECT_ROOT}/results/direct_tate_mc500_b5000"
if [[ "$#" -ne 3 ]]; then
  echo "usage: $0 GRID_A_RAW GRID_B_RAW OUTPUT_DIR" >&2
  exit 1
fi

PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_EXPECT_NUISANCE_CV_THREADS="${ROCE_EXPECT_NUISANCE_CV_THREADS:-5}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/compare_nlambda_validation.R "$@"
