#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
SEED="${2:-}"
GRID_A="${3:-20}"
GRID_B="${4:-50}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
roce_resolve_nuisance_cv_threads 5 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"

if [[ ! "${SEED}" =~ ^[1-5]$ ]]; then
  echo "usage: $0 [RESULT_ROOT] SEED [GRID_A GRID_B]" >&2
  echo "SEED must be the next paired validation seed in 1:5." >&2
  exit 1
fi
PAIR="${GRID_A}:${GRID_B}"
if [[ "${PAIR}" != "20:50" && "${PAIR}" != "50:100" ]]; then
  echo "GRID_A:GRID_B must be 20:50 or 50:100." >&2
  exit 1
fi
CPUS_PER_TASK="${ROCE_CPUS_PER_TASK:-$((2 * 4 * NUISANCE_CV_THREADS))}"
roce_require_nuisance_cv_cpu_capacity \
  "${CPUS_PER_TASK}" "${NUISANCE_CV_THREADS}"

MANIFEST_A="${RESULT_ROOT}/manifest_nlambda${GRID_A}_validation.csv"
MANIFEST_B="${RESULT_ROOT}/manifest_nlambda${GRID_B}_validation.csv"
RAW_A="${RESULT_ROOT}/nlambda${GRID_A}_validation_raw"
RAW_B="${RESULT_ROOT}/nlambda${GRID_B}_validation_raw"
for manifest in "${MANIFEST_A}" "${MANIFEST_B}"; do
  if [[ ! -f "${manifest}" ]]; then
    echo "validation manifest not found: ${manifest}" >&2
    exit 1
  fi
done

if [[ "${SEED}" -gt 1 ]]; then
  for prior_seed in $(seq 1 $((SEED - 1))); do
    task_name="$(printf 'task_%06d.csv' "${prior_seed}")"
    if [[ ! -s "${RAW_A}/${task_name}" || ! -s "${RAW_B}/${task_name}" ]]; then
      echo "paired ${GRID_A}:${GRID_B} seed ${prior_seed} is incomplete; refusing to advance" >&2
      exit 1
    fi
  done
  prior_comparison_root="${RESULT_ROOT}/nlambda${GRID_A}_vs_${GRID_B}_seed1_to_$((SEED - 1))"
  if [[ ! -s "${prior_comparison_root}/implementation_audit_passed.txt" ]]; then
    echo "prior paired implementation audit has not passed: ${prior_comparison_root}" >&2
    exit 1
  fi
fi
current_task="$(printf 'task_%06d.csv' "${SEED}")"
if [[ -e "${RAW_A}/${current_task}" || -e "${RAW_B}/${current_task}" ]]; then
  echo "seed ${SEED} already has at least one output; refusing to resubmit" >&2
  exit 1
fi
comparison_root="${RESULT_ROOT}/nlambda${GRID_A}_vs_${GRID_B}_seed1_to_${SEED}"
if [[ -e "${comparison_root}" ]]; then
  echo "comparison output already exists: ${comparison_root}" >&2
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
  echo "dry run: paired validation seed ${SEED}, nlambda=${GRID_A} and ${GRID_B}"
  echo "dry run: exactly two one-element arrays, ${CPUS_PER_TASK} CPUs and 8G each"
  echo "dry run: nuisance CV uses ${NUISANCE_CV_THREADS} OpenMP threads"
  echo "dry run: dependent audit output ${comparison_root}"
  echo "dry run: tested package library ${PROJECT_LIBRARY}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p "${RAW_A}" "${RAW_B}" "${RESULT_ROOT}/logs"

job_a="$({
  ROCE_MANIFEST="${MANIFEST_A}" \
  ROCE_OUTPUT_ROOT="${RAW_A}" \
  ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}" \
  sbatch --parsable --array="${SEED}-${SEED}%1" \
    --cpus-per-task="${CPUS_PER_TASK}" --mem=8G --time=12:00:00 \
    scripts/slurm/run_direct_tate_array.sh
})"
job_b="$({
  ROCE_MANIFEST="${MANIFEST_B}" \
  ROCE_OUTPUT_ROOT="${RAW_B}" \
  ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}" \
  sbatch --parsable --array="${SEED}-${SEED}%1" \
    --cpus-per-task="${CPUS_PER_TASK}" --mem=8G --time=12:00:00 \
    scripts/slurm/run_direct_tate_array.sh
})"
job_a="${job_a%%;*}"
job_b="${job_b%%;*}"
if [[ ! "${job_a}" =~ ^[0-9]+$ || ! "${job_b}" =~ ^[0-9]+$ ]]; then
  echo "unexpected sbatch job IDs: nlambda${GRID_A}=${job_a}, nlambda${GRID_B}=${job_b}" >&2
  exit 1
fi

comparison_job="$({
  sbatch --parsable --dependency="afterok:${job_a}:${job_b}" \
    scripts/slurm/compare_nlambda_validation.sh \
    "${RAW_A}" "${RAW_B}" "${comparison_root}"
})"
comparison_job="${comparison_job%%;*}"
if [[ ! "${comparison_job}" =~ ^[0-9]+$ ]]; then
  echo "unexpected comparison job ID: ${comparison_job}" >&2
  exit 1
fi

echo "submitted paired seed ${SEED}: nlambda${GRID_A}=${job_a}, nlambda${GRID_B}=${job_b}"
echo "submitted dependent paired audit: ${comparison_job}"
echo "resource plan: ${CPUS_PER_TASK} CPUs, ${NUISANCE_CV_THREADS} CV threads"
echo "tested package library: ${PROJECT_LIBRARY}"
