#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"
MANIFEST_PATH="${RESULT_ROOT}/manifest_nlambda100_validation.csv"
RAW_50="${RESULT_ROOT}/nlambda50_validation_raw"
RAW_100="${RESULT_ROOT}/nlambda100_validation_raw"
PREREQUISITE_SENTINEL="${RESULT_ROOT}/nlambda20_vs_50_seed1/implementation_audit_passed.txt"
COMPARISON_ROOT="${RESULT_ROOT}/nlambda50_vs_100_corrected_seed1"

for required_file in \
  "${MANIFEST_PATH}" \
  "${RAW_50}/task_000001.csv" \
  "${PREREQUISITE_SENTINEL}"; do
  if [[ ! -s "${required_file}" ]]; then
    echo "required audited seed-1 artifact is missing: ${required_file}" >&2
    exit 1
  fi
done

mapfile -t GRID_50_RESULTS < <(
  find "${RAW_50}" -maxdepth 1 -type f -name 'task_*.csv' -print | sort
)
if [[ "${#GRID_50_RESULTS[@]}" -ne 1 || \
      "${GRID_50_RESULTS[0]}" != "${RAW_50}/task_000001.csv" ]]; then
  echo "nlambda=50 directory must contain exactly seed 1 before this paired spot check" >&2
  exit 1
fi
if [[ -e "${RAW_100}/task_000001.csv" || -e "${COMPARISON_ROOT}" ]]; then
  echo "corrected nlambda=100 spot check already has an output; refusing to resubmit" >&2
  exit 1
fi

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"

SUBMIT_DRY_RUN="${ROCE_SUBMIT_DRY_RUN:-0}"
if [[ ! "${SUBMIT_DRY_RUN}" =~ ^[01]$ ]]; then
  echo "ROCE_SUBMIT_DRY_RUN must be 0 or 1: ${SUBMIT_DRY_RUN}" >&2
  exit 1
fi
if [[ "${SUBMIT_DRY_RUN}" == "1" ]]; then
  echo "dry run: corrected-code nlambda=100 seed-1 spot check"
  echo "dry run: exactly one one-element array, 8 CPUs and 8G"
  echo "dry run: dependent 50:100 audit output ${COMPARISON_ROOT}"
  echo "dry run: tested package library ${PROJECT_LIBRARY}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p "${RAW_100}" "${RESULT_ROOT}/logs"
JOB_100="$({
  ROCE_MANIFEST="${MANIFEST_PATH}" \
  ROCE_OUTPUT_ROOT="${RAW_100}" \
  sbatch --parsable --array=1-1%1 --cpus-per-task=8 --mem=8G \
    --time=12:00:00 \
    scripts/slurm/run_direct_tate_array.sh
})"
JOB_100="${JOB_100%%;*}"
if [[ ! "${JOB_100}" =~ ^[0-9]+$ ]]; then
  echo "unexpected nlambda=100 job ID: ${JOB_100}" >&2
  exit 1
fi

COMPARISON_JOB="$({
  sbatch --parsable --dependency="afterok:${JOB_100}" \
    scripts/slurm/compare_nlambda_validation.sh \
    "${RAW_50}" "${RAW_100}" "${COMPARISON_ROOT}"
})"
COMPARISON_JOB="${COMPARISON_JOB%%;*}"
if [[ ! "${COMPARISON_JOB}" =~ ^[0-9]+$ ]]; then
  echo "unexpected comparison job ID: ${COMPARISON_JOB}" >&2
  exit 1
fi

echo "submitted corrected-code nlambda=100 seed-1 spot check: ${JOB_100}"
echo "submitted dependent nlambda=50-vs-100 audit: ${COMPARISON_JOB}"
echo "tested package library: ${PROJECT_LIBRARY}"
