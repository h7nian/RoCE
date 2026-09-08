#!/bin/bash
#SBATCH --job-name=roce_checkpoint
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_checkpoint.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_checkpoint.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
MANIFEST="${2:-${ROOT}/manifest_main.csv}"

if [[ "$#" -lt 6 || "$#" -gt 7 ]]; then
  echo "usage: $0 ROOT MANIFEST CONFIG K RHO N_REPLICATIONS [CUTOFF]" >&2
  exit 1
fi

PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${ROOT}/Rlib_current}"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
MANIFEST_FINGERPRINT="$(sha256sum "$2" | awk '{print $1}')"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_MANIFEST_FINGERPRINT="${MANIFEST_FINGERPRINT}"
export ROCE_EXPECT_NUISANCE_CV_THREADS="${ROCE_EXPECT_NUISANCE_CV_THREADS:-5}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${ROOT}/checkpoints" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

Rscript scripts/slurm/audit_direct_tate_checkpoint.R "$@"
