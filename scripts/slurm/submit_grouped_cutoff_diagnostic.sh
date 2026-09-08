#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_PATH="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot/manifest.csv}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot/raw}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_candidate_20260824_v32_cv_tail_font20}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
source "${PROJECT_ROOT}/scripts/slurm/array_batch_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi

N_TASKS="$(awk 'END { print NR - 1 }' "${MANIFEST_PATH}")"
if [[ ! "${N_TASKS}" =~ ^[1-9][0-9]*$ ]] || [[ "${N_TASKS}" -gt 25 ]]; then
  echo "grouped cutoff manifest must contain 1 through 25 jobs." >&2
  exit 1
fi
if ! PRIMARY_CUTOFF="$(awk -F, '
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      name = $i; gsub(/"/, "", name)
      if (name == "primary_cutoff") column = i
    }
    next
  }
  column > 0 {
    value = $column; gsub(/"/, "", value)
    if (first == "") first = value
    if (value != first) exit 2
  }
  END { if (column > 0 && first != "") print first; else exit 1 }
' "${MANIFEST_PATH}")"; then
  echo "manifest has no unique primary_cutoff." >&2
  exit 1
fi
if [[ ! "${PRIMARY_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v cutoff="${PRIMARY_CUTOFF}" 'BEGIN { exit !(cutoff > 0) }'; then
  echo "manifest has no unique positive primary_cutoff." >&2
  exit 1
fi
if ! awk -F, -v primary_cutoff="${PRIMARY_CUTOFF}" '
  NR == 1 { next }
  {
    p=$5; gsub(/"/, "", p)
    k=$6; gsub(/"/, "", k)
    r=$7; gsub(/"/, "", r)
    c=$8; gsub(/"/, "", c)
    if (p != 100 || (k != 2 && k != 4 && k != 8) ||
        r !~ /^0;/) exit 1
    n = split(c, values, ";")
    found = 0
    for (i = 1; i <= n; i++) {
      if (values[i] == primary_cutoff) found = 1
    }
    if (!found) exit 1
  }
' "${MANIFEST_PATH}"; then
  echo "manifest failed the p=100, K, rho-grid, or primary-cutoff safety audit." >&2
  exit 1
fi

SOURCE_COUNT="$(awk -F, 'NR == 2 { v=$6; gsub(/"/, "", v); print v }' "${MANIFEST_PATH}")"
if [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS:-2}"
else
  NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS:-5}"
fi
CPUS_PER_TASK=$((2 * SOURCE_COUNT * NUISANCE_CV_THREADS))
POSITIVE_RHO_WORKERS=2
if [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-32G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-12:00:00}"
elif [[ "${SOURCE_COUNT}" -eq 4 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-16G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-08:00:00}"
else
  # Completed p=100, K=2 production jobs peaked at 3.31--3.47 GB.  Eight GB
  # retains more than twofold headroom while improving queue fit on MSI.
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-8G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-04:00:00}"
fi

resolve_roce_array_start_from_outputs "${MANIFEST_PATH}" "${OUTPUT_ROOT}"
if [[ "${ROCE_ALL_OUTPUTS_COMPLETE}" == "1" ]]; then
  echo "all ${N_TASKS} grouped cutoff jobs already have nonempty outputs"
  exit 0
fi
configure_roce_array_batch "${N_TASKS}" 5 2
enforce_roce_array_safety_cap 5 2

if [[ "${ROCE_SUBMIT_DRY_RUN:-0}" == "1" ]]; then
  echo "dry run: rows ${ROCE_BATCH_START}-${ROCE_BATCH_END}/${N_TASKS}"
  echo "dry run: array ${ROCE_ARRAY_SPEC}; p=100 K=${SOURCE_COUNT}"
  echo "dry run: ${CPUS_PER_TASK} CPUs, ${NUISANCE_CV_THREADS} CV threads, ${MEMORY_PER_TASK}, ${TIME_PER_TASK}"
  echo "dry run: tested package ${PROJECT_LIBRARY}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" results/direct_tate_mc500_b5000/logs
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export ROCE_POSITIVE_RHO_WORKERS="${POSITIVE_RHO_WORKERS}"
ROCE_MANIFEST="${MANIFEST_PATH}" \
ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}" \
sbatch --array="${ROCE_ARRAY_SPEC}" \
  --cpus-per-task="${CPUS_PER_TASK}" \
  --mem="${MEMORY_PER_TASK}" \
  --time="${TIME_PER_TASK}" \
  scripts/slurm/run_grouped_cutoff_diagnostic_array.sh

echo "submitted grouped cutoff rows ${ROCE_BATCH_START}-${ROCE_BATCH_END}/${N_TASKS}"
echo "resource plan: ${CPUS_PER_TASK} CPUs, ${MEMORY_PER_TASK}, ${TIME_PER_TASK}; max concurrency ${ROCE_MAX_CONCURRENT_RESOLVED}"
