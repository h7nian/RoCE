#!/bin/bash
#SBATCH --job-name=roce_rhc_tate
# Ten outer folds are substantially more expensive than the five-fold
# simulations because both treatment arms and all comparison methods are fit.
# The 24-hour limit avoids losing a nearly complete single RHC analysis.
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=30
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rhc.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rhc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/rhc_supported_primary}"
PACKAGE_CHECK_GATE="${ROCE_PACKAGE_CHECK_GATE:?ROCE_PACKAGE_CHECK_GATE is required}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
PACKAGE_SOURCE_FINGERPRINT="$(roce_package_source_fingerprint "${PROJECT_ROOT}")"
PACKAGE_TEST_GATE="${PROJECT_LIBRARY}/audit_tests_passed.txt"
TEST_SUITE_FINGERPRINT="$({
  find "${PROJECT_ROOT}/tests" -type f -name '*.R' -print0 \
    | sort -z | xargs -0 sha256sum
} | sha256sum | awk '{print $1}')"

# A present installation is not evidence that these exact bytes passed tests.
# Check before loading R or claiming/creating an analysis output directory.
require_gate_value() {
  local gate="$1" key="$2" expected="$3" observed
  if [[ ! -f "${gate}" ]]; then
    echo "required RHC package gate is missing: ${gate}" >&2
    return 1
  fi
  observed="$(awk -F= -v key="${key}" \
    '$1 == key {print substr($0, length($1) + 2)}' "${gate}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "RHC package gate mismatch: ${gate}, field ${key}" >&2
    return 1
  fi
}
require_gate_value "${PACKAGE_TEST_GATE}" package_tests passed
require_gate_value "${PACKAGE_TEST_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value "${PACKAGE_TEST_GATE}" package_source_fingerprint "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value "${PACKAGE_TEST_GATE}" test_suite_fingerprint "${TEST_SUITE_FINGERPRINT}"
require_gate_value "${PACKAGE_TEST_GATE}" test_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
require_gate_value "${PACKAGE_CHECK_GATE}" r_cmd_check passed
require_gate_value "${PACKAGE_CHECK_GATE}" package_source_fingerprint "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value "${PACKAGE_CHECK_GATE}" check_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"

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
export ROCE_RHC_WORKFLOW_FINGERPRINT="${RHC_WORKFLOW_FINGERPRINT}"
export ROCE_RHC_DATA_FINGERPRINT="${RHC_DATA_FINGERPRINT}"
export ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}"
export ROCE_RHC_TOTAL_SITES="${ROCE_RHC_TOTAL_SITES:-4}"
export ROCE_RHC_EXCLUDE_SITES="${ROCE_RHC_EXCLUDE_SITES-Medicare & Medicaid,No insurance}"
export ROCE_N_BOOTSTRAP="${ROCE_N_BOOTSTRAP:-5000}"
export ROCE_NLAMBDA_INIT="${ROCE_NLAMBDA_INIT:-100}"
export ROCE_NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS:-5}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
echo "[resource] cpus=${SLURM_CPUS_PER_TASK} nuisance_cv_threads=${ROCE_NUISANCE_CV_THREADS}"
echo "[package] library=${PROJECT_LIBRARY} fingerprint=${PACKAGE_FINGERPRINT}"
echo "[workflow] fingerprint=${RHC_WORKFLOW_FINGERPRINT}"
echo "[data] rhc_sha256=${RHC_DATA_FINGERPRINT}"
echo "[analysis] nlambda=${ROCE_NLAMBDA_INIT} nuisance_rule=${ROCE_NUISANCE_LAMBDA_RULE:-min} bootstrap=${ROCE_N_BOOTSTRAP}"
Rscript scripts/slurm/run_rhc_direct_tate.R
