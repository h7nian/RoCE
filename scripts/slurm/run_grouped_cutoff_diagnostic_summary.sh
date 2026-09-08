#!/bin/bash
#SBATCH --job-name=roce_group_cutoff_summary
#SBATCH --time=00:15:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_group_cutoff_summary.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_group_cutoff_summary.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
DIAGNOSTIC_ROOT="${1:?DIAGNOSTIC_ROOT is required}"
EXPECTED_REPLICATIONS="${2:?EXPECTED_REPLICATIONS is required}"
IDENTITY_MODE="${3:-require-production-identity}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_candidate_20260824_v32_cv_tail_font20}"
PRIMARY_ROOT="${PROJECT_ROOT}/results/direct_tate_mc500_b5000"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
SUMMARY_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/run_grouped_cutoff_diagnostic_summary.sh" \
  "${PROJECT_ROOT}/scripts/slurm/summarize_grouped_cutoff_diagnostic.R" \
  "${PROJECT_ROOT}/scripts/slurm/grouped_cutoff_diagnostic_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/atomic_output.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_CUTOFF_SUMMARY_WORKFLOW_FINGERPRINT="${SUMMARY_WORKFLOW_FINGERPRINT}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
echo "[summary workflow] fingerprint=${SUMMARY_WORKFLOW_FINGERPRINT}"
Rscript scripts/slurm/summarize_grouped_cutoff_diagnostic.R \
  "${DIAGNOSTIC_ROOT}/raw" \
  "${DIAGNOSTIC_ROOT}/manifest.csv" \
  "${EXPECTED_REPLICATIONS}" \
  "${PRIMARY_ROOT}/raw" \
  "${PRIMARY_ROOT}/manifest_main.csv" \
  "${IDENTITY_MODE}"
