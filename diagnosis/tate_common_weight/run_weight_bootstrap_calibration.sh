#!/bin/bash
#SBATCH --job-name=roce_weight_calibration
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_weight_calibration.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_weight_calibration.err

set -euo pipefail

if [[ "$#" -ne 5 ]]; then
  echo "usage: run_weight_bootstrap_calibration.sh GROUP_MANIFEST PRIMARY_MANIFEST GROUP_TASK_ID OUTPUT_DIR B" >&2
  exit 1
fi
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"
CHECK_GATE="${ROCE_PACKAGE_CHECK_GATE:?ROCE_PACKAGE_CHECK_GATE is required}"
GROUP_MANIFEST="$(readlink -f "$1")"
PRIMARY_MANIFEST="$(readlink -f "$2")"
GROUP_TASK_ID="$3"
OUTPUT_DIRECTORY="$4"
BOOTSTRAP_DRAWS="$5"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
TEST_GATE="${PROJECT_LIBRARY}/audit_tests_passed.txt"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
SOURCE_FINGERPRINT="$(roce_package_source_fingerprint "${PROJECT_ROOT}")"

require_gate_value() {
  local gate="$1" key="$2" expected="$3" observed
  if [[ ! -f "${gate}" ]]; then
    echo "required package gate is missing: ${gate}" >&2
    exit 1
  fi
  observed="$(awk -F= -v key="${key}" '$1 == key {print substr($0, length($1) + 2)}' "${gate}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "package gate mismatch: ${gate}, field ${key}" >&2
    exit 1
  fi
}

require_gate_value "${TEST_GATE}" package_tests passed
require_gate_value "${TEST_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value "${TEST_GATE}" package_source_fingerprint "${SOURCE_FINGERPRINT}"
require_gate_value "${TEST_GATE}" test_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
require_gate_value "${CHECK_GATE}" r_cmd_check passed
require_gate_value "${CHECK_GATE}" package_source_fingerprint "${SOURCE_FINGERPRINT}"
require_gate_value "${CHECK_GATE}" check_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/diagnosis/tate_common_weight/run_weight_bootstrap_calibration.sh" \
  "${PROJECT_ROOT}/diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R" \
  "${PROJECT_ROOT}/scripts/slurm/atomic_output.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/resource_topology.R" \
  "${PROJECT_ROOT}/scripts/slurm/simulation_qc_policy.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"
export ROCE_MANIFEST_FINGERPRINT="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
export ROCE_GROUP_MANIFEST_FINGERPRINT="$(sha256sum "${GROUP_MANIFEST}" | awk '{print $1}')"
export ROCE_NUISANCE_CV_THREADS=5
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
Rscript diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R \
  "${GROUP_MANIFEST}" "${PRIMARY_MANIFEST}" "${GROUP_TASK_ID}" \
  "${OUTPUT_DIRECTORY}" "${BOOTSTRAP_DRAWS}"
