#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_PATH="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/manifest_cutoff_diagnostic.csv}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/cutoff_diagnostic/raw}"
CPUS_PER_TASK="${3:-${ROCE_CPUS_PER_TASK:-}}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
roce_resolve_nuisance_cv_threads 5 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"

if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi

N_TASKS="$(awk 'END { print NR - 1 }' "${MANIFEST_PATH}")"
if [[ "${N_TASKS}" -lt 1 ]]; then
  echo "manifest has no task rows: ${MANIFEST_PATH}" >&2
  exit 1
fi
if [[ -z "${CPUS_PER_TASK}" ]]; then
  CPUS_PER_TASK=$((8 * NUISANCE_CV_THREADS))
fi
roce_require_nuisance_cv_cpu_capacity \
  "${CPUS_PER_TASK}" "${NUISANCE_CV_THREADS}"

source "${PROJECT_ROOT}/scripts/slurm/array_batch_utils.sh"
# Cutoff jobs reuse nuisance fits but remain computationally expensive. Keep
# the same conservative five-task, two-concurrent default as the main batches.
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
  echo "dry run: Slurm array ${ROCE_ARRAY_SPEC}, ${CPUS_PER_TASK} CPUs"
  echo "dry run: nuisance CV uses ${NUISANCE_CV_THREADS} OpenMP threads"
  echo "dry run: tested package library ${PROJECT_LIBRARY}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" results/direct_tate_mc500_b5000/logs
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"

ROCE_MANIFEST="${MANIFEST_PATH}" \
ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}" \
sbatch --array="${ROCE_ARRAY_SPEC}" \
  --cpus-per-task="${CPUS_PER_TASK}" \
  scripts/slurm/run_direct_tate_cutoff_array.sh

echo "submitted manifest rows ${ROCE_BATCH_START}-${ROCE_BATCH_END} of ${N_TASKS}"
echo "resource plan: ${CPUS_PER_TASK} CPUs, ${NUISANCE_CV_THREADS} CV threads"
