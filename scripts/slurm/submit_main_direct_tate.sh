#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw}"
SOURCE_COUNT="${3:-}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
SMOKE_GATE="${ROCE_SMOKE_GATE:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/smoke_final_direct_tate_raw/direct_tate_smoke_audit_passed.txt}"
MANIFEST_GATE="${ROCE_MANIFEST_GATE:-${MANIFEST_ROOT}/manifest_audit_passed.txt}"
CUTOFF_SELECTION_GATE="${ROCE_CUTOFF_SELECTION_GATE:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot/cutoff_decision_n010/cutoff_selection_passed.txt}"
SMOKE_RESULT="${ROCE_SMOKE_RESULT:-$(dirname "${SMOKE_GATE}")/task_000001.csv}"
SMOKE_SENSITIVITY_RESULT="${ROCE_SMOKE_SENSITIVITY_RESULT:-$(dirname "${SMOKE_GATE}")_reused_sensitivity/task_000001.csv}"

source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"

if [[ ! "${SOURCE_COUNT}" =~ ^(2|4|8)$ ]]; then
  echo "usage: $0 [MANIFEST_ROOT] [OUTPUT_ROOT] K" >&2
  echo "K must be one of 2, 4, or 8; one bounded batch is submitted per call." >&2
  exit 1
fi
if [[ -n "${ROCE_NUISANCE_CV_THREADS:-}" ]]; then
  roce_resolve_nuisance_cv_threads 1 5
elif [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  ROCE_NUISANCE_CV_THREADS=2
  export ROCE_NUISANCE_CV_THREADS
  roce_resolve_nuisance_cv_threads 2 5
else
  roce_resolve_nuisance_cv_threads 5 5
fi
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"

MANIFEST_PATH="${MANIFEST_ROOT}/manifest_main_K${SOURCE_COUNT}.csv"
if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "split manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi
if [[ ! -f "${MANIFEST_GATE}" ]]; then
  echo "audited-manifest gate not found: ${MANIFEST_GATE}" >&2
  echo "regenerate and audit the MC500/B5000 manifests before submission" >&2
  exit 1
fi
if [[ ! -f "${CUTOFF_SELECTION_GATE}" ]]; then
  echo "audited cutoff-selection gate not found: ${CUTOFF_SELECTION_GATE}" >&2
  exit 1
fi
MANIFEST_AUDIT_DRIVER_MD5="$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_mc500_manifests.R" | awk '{print $1}')"
MANIFEST_MD5="$(md5sum "${MANIFEST_PATH}" | awk '{print $1}')"
if ! grep -Fqx -- "manifest_audit=passed" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- "replications_per_setting=500" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- "bootstrap_draws=5000" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- "p=100" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- "audit_driver_md5=${MANIFEST_AUDIT_DRIVER_MD5}" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- "manifest_main_K${SOURCE_COUNT}_md5=${MANIFEST_MD5}" "${MANIFEST_GATE}"; then
  echo "audited-manifest gate does not match ${MANIFEST_PATH}" >&2
  exit 1
fi
SELECTED_CUTOFF="$(awk -F= '$1 == "selected_cutoff" {
  print substr($0, length($1) + 2)
}' "${CUTOFF_SELECTION_GATE}")"
CUTOFF_SELECTION_FINGERPRINT="$(sha256sum \
  "${CUTOFF_SELECTION_GATE}" | awk '{print $1}')"
if [[ "$(awk -F= '$1 == "cutoff_selection" {
       print substr($0, length($1) + 2)
     }' "${CUTOFF_SELECTION_GATE}")" != "passed" ]] ||
   [[ ! "${SELECTED_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! grep -Fqx -- "primary_cutoff=${SELECTED_CUTOFF}" "${MANIFEST_GATE}" ||
   ! grep -Fqx -- \
     "cutoff_selection_fingerprint=${CUTOFF_SELECTION_FINGERPRINT}" \
     "${MANIFEST_GATE}"; then
  echo "manifest and cutoff-selection gates are inconsistent." >&2
  exit 1
fi

roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
AUDIT_DRIVER_MD5="$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_direct_tate_smoke.R" | awk '{print $1}')"
SUBMIT_DRIVER_MD5="$(md5sum "${PROJECT_ROOT}/scripts/slurm/submit_main_direct_tate.sh" | awk '{print $1}')"
EXPECTED_NLAMBDA="$(awk -F, '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      name = $i; gsub(/"/, "", name)
      if (name == "nlambda_init") column = i
    }
    next
  }
  NR == 2 && column > 0 {
    value = $column; gsub(/"/, "", value); print value; exit
  }
' "${MANIFEST_PATH}")"
if [[ ! "${EXPECTED_NLAMBDA}" =~ ^[0-9]+$ ]] ||
   [[ "${EXPECTED_NLAMBDA}" -lt 2 ]]; then
  echo "could not resolve a valid nlambda_init from ${MANIFEST_PATH}" >&2
  exit 1
fi
if [[ ! -f "${SMOKE_GATE}" ]]; then
  echo "production smoke gate not found: ${SMOKE_GATE}" >&2
  echo "run and audit the final p=100 smoke with the selected package/grid first" >&2
  exit 1
fi
if [[ ! -f "${SMOKE_RESULT}" || ! -f "${SMOKE_SENSITIVITY_RESULT}" ]]; then
  echo "smoke gate inputs are missing: ${SMOKE_RESULT} or ${SMOKE_SENSITIVITY_RESULT}" >&2
  exit 1
fi
require_gate_value() {
  local key="$1"
  local expected="$2"
  if ! grep -Fqx -- "${key}=${expected}" "${SMOKE_GATE}"; then
    echo "smoke gate mismatch for ${key}; expected ${expected}" >&2
    exit 1
  fi
}
require_gate_value implementation_audit passed
require_gate_value package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value workflow_fingerprint "${WORKFLOW_FINGERPRINT}"
require_gate_value audit_driver_md5 "${AUDIT_DRIVER_MD5}"
require_gate_value submit_driver_md5 "${SUBMIT_DRIVER_MD5}"
require_gate_value smoke_result_md5 "$(md5sum "${SMOKE_RESULT}" | awk '{print $1}')"
require_gate_value smoke_sensitivity_result_md5 "$(md5sum "${SMOKE_SENSITIVITY_RESULT}" | awk '{print $1}')"
require_gate_value nlambda_init "${EXPECTED_NLAMBDA}"
EXPECTED_PRIMARY_CUTOFF="$(awk -F, '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      name = $i; gsub(/"/, "", name)
      if (name == "cutoff") column = i
    }
    next
  }
  NR == 2 && column > 0 {
    value = $column; gsub(/"/, "", value); print value; exit
  }
' "${MANIFEST_PATH}")"
if [[ ! "${EXPECTED_PRIMARY_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v cutoff="${EXPECTED_PRIMARY_CUTOFF}" \
      'BEGIN { exit !(cutoff > 0) }'; then
  echo "could not resolve a positive primary cutoff from ${MANIFEST_PATH}" >&2
  exit 1
fi
require_gate_value primary_cutoff "${EXPECTED_PRIMARY_CUTOFF}"
if [[ "${EXPECTED_PRIMARY_CUTOFF}" != "${SELECTED_CUTOFF}" ]]; then
  echo "split manifest cutoff does not match the audited selected cutoff." >&2
  exit 1
fi
require_gate_value p 100
require_gate_value config C3
require_gate_value K 4
require_gate_value rho 0

if [[ "${ROCE_PARALLEL_ARMS:-1}" == "1" ]]; then
  CPUS_PER_TASK=$((2 * SOURCE_COUNT * NUISANCE_CV_THREADS))
else
  CPUS_PER_TASK=$((SOURCE_COUNT * NUISANCE_CV_THREADS))
fi
# The K=8 job can fork 16 R workers across the two treatment arms. Give that
# configuration additional headroom while keeping the smaller pilots at 8G.
if [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-16G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-24:00:00}"
elif [[ "${SOURCE_COUNT}" -eq 4 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-8G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-12:00:00}"
else
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-8G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-08:00:00}"
fi

export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
"${PROJECT_ROOT}/scripts/slurm/submit_direct_tate.sh" \
  "${MANIFEST_PATH}" "${OUTPUT_ROOT}" "${CPUS_PER_TASK}" \
  "${MEMORY_PER_TASK}" "${TIME_PER_TASK}"
