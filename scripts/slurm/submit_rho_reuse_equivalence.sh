#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
VALIDATION_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/rho_reuse_equivalence_candidate}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB must name the final tested library}"
PACKAGE_TEST_GATE="${ROCE_PACKAGE_TEST_GATE:-${PROJECT_LIBRARY}/audit_tests_passed.txt}"
PACKAGE_CHECK_GATE="${ROCE_PACKAGE_CHECK_GATE:?ROCE_PACKAGE_CHECK_GATE is required}"
PRIMARY_MANIFEST="${MANIFEST_ROOT}/manifest_main.csv"
GROUP_MANIFEST="${MANIFEST_ROOT}/manifest_main_rho_groups.csv"
MANIFEST_GATE="${ROCE_MANIFEST_GATE:-${MANIFEST_ROOT}/manifest_audit_passed.txt}"
SLURM_PARTITION="${ROCE_SLURM_PARTITION:-msismall}"
SUBMIT_DRY_RUN="${ROCE_SUBMIT_DRY_RUN:-0}"
if [[ ! "${SUBMIT_DRY_RUN}" =~ ^[01]$ ]]; then
  echo "ROCE_SUBMIT_DRY_RUN must be 0 or 1." >&2
  exit 1
fi
if [[ ! "${SLURM_PARTITION}" =~ ^[A-Za-z0-9_-]+$ ]]; then
  echo "ROCE_SLURM_PARTITION contains invalid characters." >&2
  exit 1
fi

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
for required_path in \
  "${PRIMARY_MANIFEST}" "${GROUP_MANIFEST}" "${MANIFEST_GATE}" \
  "${PACKAGE_TEST_GATE}" "${PACKAGE_CHECK_GATE}"; do
  if [[ ! -f "${required_path}" ]]; then
    echo "reuse validation input/gate is missing: ${required_path}" >&2
    exit 1
  fi
done
if [[ -e "${VALIDATION_ROOT}" ]]; then
  echo "refusing to mix reuse-validation artifacts: ${VALIDATION_ROOT}" >&2
  exit 2
fi

roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
PACKAGE_SOURCE_FINGERPRINT="$(roce_package_source_fingerprint "${PROJECT_ROOT}")"
TEST_SUITE_FINGERPRINT="$({
  find "${PROJECT_ROOT}/tests" -type f -name '*.R' -print0 \
    | sort -z | xargs -0 sha256sum
} | sha256sum | awk '{print $1}')"
require_gate_value() {
  local gate_path="$1"
  local key="$2"
  local expected="$3"
  if ! grep -Fqx -- "${key}=${expected}" "${gate_path}"; then
    echo "gate mismatch in ${gate_path}: ${key}=${expected}" >&2
    exit 1
  fi
}
require_gate_value "${MANIFEST_GATE}" manifest_audit passed
require_gate_value "${MANIFEST_GATE}" replications_per_setting 500
require_gate_value "${MANIFEST_GATE}" bootstrap_draws 5000
require_gate_value "${MANIFEST_GATE}" p 100
require_gate_value \
  "${MANIFEST_GATE}" manifest_main_md5 \
  "$(md5sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
require_gate_value \
  "${MANIFEST_GATE}" manifest_main_rho_groups_md5 \
  "$(md5sum "${GROUP_MANIFEST}" | awk '{print $1}')"
require_gate_value \
  "${MANIFEST_GATE}" audit_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_mc500_manifests.R" | awk '{print $1}')"
require_gate_value "${PACKAGE_TEST_GATE}" package_tests passed
require_gate_value \
  "${PACKAGE_TEST_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" package_source_fingerprint \
  "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" test_suite_fingerprint "${TEST_SUITE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" test_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
require_gate_value "${PACKAGE_CHECK_GATE}" r_cmd_check passed
require_gate_value \
  "${PACKAGE_CHECK_GATE}" package_source_fingerprint \
  "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_CHECK_GATE}" check_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"

INDEPENDENT_TASK_ID="$(awk -F, '
  NR > 1 {
    for (i = 1; i <= 7; i++) gsub(/"/, "", $i)
    if ($3 == 1 && $4 == "C3" && $5 == 100 && $6 == 4 && $7 == 2.5) {
      print $1
    }
  }
' "${PRIMARY_MANIFEST}")"
GROUP_TASK_ID="$(awk -F, '
  NR > 1 {
    for (i = 1; i <= 6; i++) gsub(/"/, "", $i)
    if ($3 == 1 && $4 == "C3" && $5 == 100 && $6 == 4) print $1
  }
' "${GROUP_MANIFEST}")"
if [[ ! "${INDEPENDENT_TASK_ID}" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "${GROUP_TASK_ID}" =~ ^[1-9][0-9]*$ ]]; then
  echo "could not resolve unique locked p=100 reuse-audit tasks." >&2
  exit 1
fi
if [[ "${SUBMIT_DRY_RUN}" == "1" ]]; then
  echo "dry run: independent primary task ${INDEPENDENT_TASK_ID}"
  echo "dry run: six-rho grouped task ${GROUP_TASK_ID}"
  echo "dry run: p=100, C3, K=4, rho=2.5, nlambda=100, B=5000"
  echo "dry run: partition=${SLURM_PARTITION}"
  echo "dry run: tested package ${PROJECT_LIBRARY}"
  exit 0
fi

INDEPENDENT_ROOT="${VALIDATION_ROOT}/independent_raw"
INDEPENDENT_SENSITIVITY_ROOT="${VALIDATION_ROOT}/independent_sensitivity"
GROUP_ROOT="${VALIDATION_ROOT}/grouped_raw"
GROUP_SENSITIVITY_ROOT="${VALIDATION_ROOT}/grouped_sensitivity"
AUDIT_ROOT="${VALIDATION_ROOT}/audit"
INDEPENDENT_RESULT="${INDEPENDENT_ROOT}/$(printf 'task_%06d.csv' "${INDEPENDENT_TASK_ID}")"
GROUPED_RESULT="${GROUP_ROOT}/$(printf 'task_%06d.csv' "${INDEPENDENT_TASK_ID}")"
mkdir -p "${INDEPENDENT_ROOT}" "${INDEPENDENT_SENSITIVITY_ROOT}" \
  "${GROUP_ROOT}" "${GROUP_SENSITIVITY_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

roce_resolve_nuisance_cv_threads 5 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"
CPUS_PER_TASK=$((2 * 4 * NUISANCE_CV_THREADS))
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
cd "${PROJECT_ROOT}"

INDEPENDENT_JOB="$(
  ROCE_MANIFEST="${PRIMARY_MANIFEST}" \
  ROCE_OUTPUT_ROOT="${INDEPENDENT_ROOT}" \
  ROCE_SENSITIVITY_OUTPUT_ROOT="${INDEPENDENT_SENSITIVITY_ROOT}" \
  sbatch --parsable --partition="${SLURM_PARTITION}" \
    --array="${INDEPENDENT_TASK_ID}-${INDEPENDENT_TASK_ID}%1" \
    --cpus-per-task="${CPUS_PER_TASK}" --mem=16G --time=18:00:00 \
    scripts/slurm/run_direct_tate_array.sh
)"
GROUP_JOB="$(
  ROCE_GROUP_MANIFEST="${GROUP_MANIFEST}" \
  ROCE_PRIMARY_MANIFEST="${PRIMARY_MANIFEST}" \
  ROCE_OUTPUT_ROOT="${GROUP_ROOT}" \
  ROCE_SENSITIVITY_OUTPUT_ROOT="${GROUP_SENSITIVITY_ROOT}" \
  sbatch --parsable --partition="${SLURM_PARTITION}" \
    --array="${GROUP_TASK_ID}-${GROUP_TASK_ID}%1" \
    --cpus-per-task="${CPUS_PER_TASK}" --mem=16G --time=24:00:00 \
    scripts/slurm/run_direct_tate_rho_group_array.sh
)"
AUDIT_JOB="$(
  ROCE_INDEPENDENT_RESULT="${INDEPENDENT_RESULT}" \
  ROCE_GROUPED_RESULT="${GROUPED_RESULT}" \
  ROCE_REUSE_AUDIT_OUTPUT="${AUDIT_ROOT}" \
  sbatch --parsable --partition="${SLURM_PARTITION}" \
    --dependency="afterok:${INDEPENDENT_JOB}:${GROUP_JOB}" \
    scripts/slurm/run_rho_reuse_equivalence_audit.sh
)"

echo "submitted independent p=100 task ${INDEPENDENT_TASK_ID}: ${INDEPENDENT_JOB}"
echo "submitted six-rho grouped task ${GROUP_TASK_ID}: ${GROUP_JOB}"
echo "submitted dependent exact-equivalence audit: ${AUDIT_JOB}"
