#!/bin/bash
#SBATCH --job-name=roce_refit_pilot
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=5
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a_nuisance_refit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a_nuisance_refit.err
set -euo pipefail

# Scheduling-only adapter for the grouped-CV single-draw implementation.
# Submit IDs with --array=2-20%N; the caller must respect the four-job total.
if [[ "$#" -ne 5 ]]; then
  echo "usage: run_full_refit_array.sh BUNDLE_DIR RHO OUTPUT_ROOT IDENTITY_AUDIT FIRST_DRAW_AUDIT" >&2
  exit 1
fi
if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^[0-9]+$ ]]; then
  echo "a Slurm array task ID is required" >&2
  exit 1
fi
DRAW_ID=$((10#${SLURM_ARRAY_TASK_ID}))
if (( DRAW_ID < 2 || DRAW_ID > 20 )); then
  echo "this bounded pilot permits only remaining draw IDs 2--20" >&2
  exit 1
fi
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
BUNDLE_DIRECTORY="$(readlink -f "$1")"
RHO_VALUE="$2"
OUTPUT_ROOT="$3"
if [[ "$(basename "${BUNDLE_DIRECTORY}")" != seed_000001 ]] ||
   [[ "${RHO_VALUE}" != 0 && "${RHO_VALUE}" != 1 ]]; then
  echo "the declared pilot is restricted to seed 1 and rho 0/1" >&2
  exit 1
fi
if [[ ! -d "${OUTPUT_ROOT}" ]]; then
  echo "OUTPUT_ROOT must already exist" >&2
  exit 1
fi
OUTPUT_ROOT="$(readlink -f "${OUTPUT_ROOT}")"

require_refit_audit() {
  local audit_directory="$1" expected_draw="$2" gate observed
  audit_directory="$(readlink -f "${audit_directory}")"
  gate="${audit_directory}/audit_passed.txt"
  (cd "${audit_directory}" && sha256sum -c sha256.txt)
  observed="$(awk -F= '$1 == "full_refit_intermediate_audit" {print $2}' "${gate}")"
  [[ "${observed}" == passed ]] || { echo "prerequisite audit did not pass" >&2; exit 1; }
  observed="$(awk -F= '$1 == "draw_bundle" {print substr($0,length($1)+2)}' "${gate}")"
  [[ "${observed}" == "${expected_draw}" ]] || { echo "prerequisite draw mismatch" >&2; exit 1; }
  observed="$(awk -F= '$1 == "source_bundle" {print substr($0,length($1)+2)}' "${gate}")"
  [[ "${observed}" == "${BUNDLE_DIRECTORY}" ]] || { echo "prerequisite source mismatch" >&2; exit 1; }
  observed="$(awk -F= '$1 == "installed_package_fingerprint" {print $2}' "${gate}")"
  [[ "${observed}" == 2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7 ]] || {
    echo "prerequisite package mismatch" >&2; exit 1;
  }
  observed="$(awk -F= '$1 == "nuisance_cv_partition" {print $2}' "${gate}")"
  [[ "${observed}" == origin_grouped ]] || {
    echo "prerequisite nuisance-CV protocol mismatch" >&2; exit 1;
  }
  observed="$(awk -F= '$1 == "nuisance_cv_partition_records" {print $2}' "${gate}")"
  [[ "${observed}" =~ ^[0-9]+$ ]] && (( 10#${observed} > 0 )) || {
    echo "prerequisite nuisance-CV partition records must be positive" >&2; exit 1;
  }
}
DRAW_PREFIX="${OUTPUT_ROOT}/seed_000001_rho_${RHO_VALUE}_draw_"
require_refit_audit "$4" "${DRAW_PREFIX}0000"
require_refit_audit "$5" "${DRAW_PREFIX}0001"
printf -v DRAW_LABEL '%04d' "${DRAW_ID}"
exec bash "${PROJECT_ROOT}/diagnosis/tate_common_weight/run_full_refit_draw.sh" \
  "${BUNDLE_DIRECTORY}" "${RHO_VALUE}" "${DRAW_ID}" "${DRAW_PREFIX}${DRAW_LABEL}"
