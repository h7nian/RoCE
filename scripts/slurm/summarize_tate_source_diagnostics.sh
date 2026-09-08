#!/bin/bash
#SBATCH --job-name=roce_tate_diag_summary
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_tate_diag_summary.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_tate_diag_summary.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
USAGE="usage: summarize_tate_source_diagnostics.sh RAW_DIR OUTPUT.csv PRIMARY_MANIFEST.csv GROUP_MANIFEST.csv CONFIG K N_REPLICATIONS CUTOFF EXPECT_HARD_THRESHOLD"
RAW_DIRECTORY="${1:?${USAGE}}"
OUTPUT_PATH="${2:?${USAGE}}"
PRIMARY_MANIFEST="${3:?${USAGE}}"
GROUP_MANIFEST="${4:?${USAGE}}"
CONFIG="${5:?${USAGE}}"
SOURCE_COUNT="${6:?${USAGE}}"
N_REPLICATIONS="${7:?${USAGE}}"
CUTOFF="${8:?${USAGE}}"
EXPECT_HARD_THRESHOLD="${9:?${USAGE}}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
export ROCE_WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
export ROCE_MANIFEST_FINGERPRINT="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
export ROCE_GROUP_MANIFEST_FINGERPRINT="$(sha256sum "${GROUP_MANIFEST}" | awk '{print $1}')"

cd "${PROJECT_ROOT}"
Rscript scripts/slurm/summarize_tate_source_diagnostics.R \
  "${RAW_DIRECTORY}" "${OUTPUT_PATH}" \
  "${PRIMARY_MANIFEST}" "${GROUP_MANIFEST}" \
  "${CONFIG}" "${SOURCE_COUNT}" "${N_REPLICATIONS}" "${CUTOFF}" \
  "${EXPECT_HARD_THRESHOLD}"
