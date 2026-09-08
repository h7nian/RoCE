#!/bin/bash
#SBATCH --job-name=roce_rhc_audit
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rhc_audit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rhc_audit.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/rhc_supported_primary}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
RHC_DATA_PATH="${PROJECT_LIBRARY}/RoCE/extdata/rhc.csv"
if [[ ! -f "${RHC_DATA_PATH}" ]]; then
  echo "installed RHC data file is missing: ${RHC_DATA_PATH}" >&2
  exit 1
fi
RHC_DATA_FINGERPRINT="$(sha256sum "${RHC_DATA_PATH}" | awk '{print $1}')"
RHC_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/run_rhc_direct_tate.sh" \
  "${PROJECT_ROOT}/scripts/slurm/run_rhc_direct_tate.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/resource_topology.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_EXPECT_RHC_DATA_FINGERPRINT="${RHC_DATA_FINGERPRINT}"
export ROCE_EXPECT_RHC_WORKFLOW_FINGERPRINT="${RHC_WORKFLOW_FINGERPRINT}"
export ROCE_EXPECT_RHC_EXCLUDED_SITES="${ROCE_EXPECT_RHC_EXCLUDED_SITES-Medicare & Medicaid,No insurance}"
export ROCE_EXPECT_N_BOOTSTRAP="${ROCE_EXPECT_N_BOOTSTRAP:-5000}"
export ROCE_EXPECT_NLAMBDA_INIT="${ROCE_EXPECT_NLAMBDA_INIT:-100}"
export ROCE_EXPECT_NUISANCE_LAMBDA_RULE="${ROCE_EXPECT_NUISANCE_LAMBDA_RULE:-min}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
Rscript scripts/slurm/audit_rhc_direct_tate.R "${OUTPUT_ROOT}"
