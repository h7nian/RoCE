#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_PATH="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/manifest.csv}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw}"
CPUS_PER_TASK="${3:-${ROCE_CPUS_PER_TASK:-4}}"
MEMORY_PER_TASK="${4:-${ROCE_MEMORY_PER_TASK:-8G}}"
TIME_PER_TASK="${5:-${ROCE_TIME_PER_TASK:-12:00:00}}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
roce_resolve_nuisance_cv_threads 1 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"

if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi
MANIFEST_PATH="$(readlink -f -- "${MANIFEST_PATH}")"
OUTPUT_ROOT="$(readlink -m -- "${OUTPUT_ROOT}")"
if ! roce_require_package_library "${PROJECT_LIBRARY}"; then
  echo "set ROCE_PROJECT_LIB to an isolated library that passed the audit suite" >&2
  exit 1
fi
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"

# Keep implementation-smoke sidecars separate from production sensitivities.
# Production uses the canonical raw/ + reused_sensitivity_raw pair; any other
# output root receives a sibling, output-specific sidecar directory unless the
# caller explicitly supplies ROCE_SENSITIVITY_OUTPUT_ROOT.
if [[ -n "${ROCE_SENSITIVITY_OUTPUT_ROOT:-}" ]]; then
  SENSITIVITY_OUTPUT_ROOT="$(readlink -m -- "${ROCE_SENSITIVITY_OUTPUT_ROOT}")"
elif [[ "${OUTPUT_ROOT%/}" == "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw" ]]; then
  SENSITIVITY_OUTPUT_ROOT="${PROJECT_ROOT}/results/direct_tate_mc500_b5000/reused_sensitivity_raw"
else
  SENSITIVITY_OUTPUT_ROOT="${OUTPUT_ROOT%/}_reused_sensitivity"
fi

N_TASKS="$(awk 'END { print NR - 1 }' "${MANIFEST_PATH}")"
if [[ "${N_TASKS}" -lt 1 ]]; then
  echo "manifest has no task rows: ${MANIFEST_PATH}" >&2
  exit 1
fi
roce_require_nuisance_cv_cpu_capacity \
  "${CPUS_PER_TASK}" "${NUISANCE_CV_THREADS}"
if [[ ! "${MEMORY_PER_TASK}" =~ ^[1-9][0-9]*[GM]$ ]]; then
  echo "memory per task must use a positive Slurm value such as 8G: ${MEMORY_PER_TASK}" >&2
  exit 1
fi
if [[ ! "${TIME_PER_TASK}" =~ ^([0-9]+-)?[0-9]{1,2}:[0-9]{2}:[0-9]{2}$ ]]; then
  echo "time per task must use a Slurm value such as 12:00:00: ${TIME_PER_TASK}" >&2
  exit 1
fi

source "${PROJECT_ROOT}/scripts/slurm/array_batch_utils.sh"
# MSI-safe default: submit only a very small pilot-sized batch.  Larger batches
# require an explicit ROCE_BATCH_SIZE / ROCE_MAX_CONCURRENT override after
# reviewing the preceding batch diagnostics.
resolve_roce_array_start_from_outputs "${MANIFEST_PATH}" "${OUTPUT_ROOT}"
if [[ "${ROCE_ALL_OUTPUTS_COMPLETE}" == "1" ]]; then
  echo "all ${N_TASKS} manifest rows already have nonempty outputs in ${OUTPUT_ROOT}"
  exit 0
fi
configure_roce_array_batch "${N_TASKS}" 5 2
enforce_roce_array_safety_cap 5 2
cap_roce_array_batch_to_manifest_setting "${MANIFEST_PATH}"

SUBMIT_DRY_RUN="${ROCE_SUBMIT_DRY_RUN:-0}"
if [[ ! "${SUBMIT_DRY_RUN}" =~ ^[01]$ ]]; then
  echo "ROCE_SUBMIT_DRY_RUN must be 0 or 1: ${SUBMIT_DRY_RUN}" >&2
  exit 1
fi
if [[ "${SUBMIT_DRY_RUN}" == "1" ]]; then
  echo "dry run: manifest rows ${ROCE_BATCH_START}-${ROCE_BATCH_END} of ${N_TASKS}"
  echo "dry run: Slurm array ${ROCE_ARRAY_SPEC}, ${CPUS_PER_TASK} CPUs, ${MEMORY_PER_TASK}, ${TIME_PER_TASK}"
  echo "dry run: nuisance CV uses ${NUISANCE_CV_THREADS} OpenMP threads"
  echo "dry run: tested package library ${PROJECT_LIBRARY}"
  echo "dry run: sensitivity sidecars ${SENSITIVITY_OUTPUT_ROOT}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p results/direct_tate_mc500_b5000/logs "${OUTPUT_ROOT}" \
  "${SENSITIVITY_OUTPUT_ROOT}"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"

ROCE_MANIFEST="${MANIFEST_PATH}" \
ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}" \
ROCE_SENSITIVITY_OUTPUT_ROOT="${SENSITIVITY_OUTPUT_ROOT}" \
sbatch --array="${ROCE_ARRAY_SPEC}" \
  --cpus-per-task="${CPUS_PER_TASK}" \
  --mem="${MEMORY_PER_TASK}" \
  --time="${TIME_PER_TASK}" \
  scripts/slurm/run_direct_tate_array.sh

echo "submitted manifest rows ${ROCE_BATCH_START}-${ROCE_BATCH_END} of ${N_TASKS} (${CPUS_PER_TASK} CPUs, ${MEMORY_PER_TASK}, ${TIME_PER_TASK})"
echo "nuisance CV threads: ${NUISANCE_CV_THREADS}"
echo "tested package library: ${PROJECT_LIBRARY}"
echo "sensitivity sidecars: ${SENSITIVITY_OUTPUT_ROOT}"
