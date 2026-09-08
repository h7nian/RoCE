#!/bin/bash
#SBATCH --job-name=roce_inference_seed_audit
#SBATCH --time=00:15:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a_independent_audit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a_independent_audit.err
set -euo pipefail

# Scheduling only: submit with aftercorr to the matching simulation array.
# The existing R auditor validates all scientific payloads and publishes its
# own atomic audit bundle. This wrapper never fits a model or changes a CI.
if [[ "$#" -ne 2 ]]; then
  echo "usage: run_independent_inference_pilot_audit_array.sh PILOT_ROOT PACKAGE_LIBRARY" >&2
  exit 1
fi
if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^([1-9]|[1-9][0-9]|100)$ ]]; then
  echo "SLURM_ARRAY_TASK_ID must be one integer in [1,100]" >&2
  exit 1
fi
TASK_ID="${SLURM_ARRAY_TASK_ID}"
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PILOT_ROOT="$(readlink -f "$1")"
if [[ ! -d "${PILOT_ROOT}" || ! -f "${PILOT_ROOT}/manifest.csv" ]]; then
  echo "pilot root or manifest is missing" >&2
  exit 1
fi
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "$2"
PACKAGE_LIBRARY="$(roce_resolve_package_library "$2")"
SIM_ID=$((10000 + TASK_ID))
printf -v SEED_LABEL 'seed_%06d' "${SIM_ID}"
BUNDLE="${PILOT_ROOT}/${SEED_LABEL}"
OUTPUT="${PILOT_ROOT}/audits/${SEED_LABEL}"
if [[ ! -d "${BUNDLE}" ]]; then
  echo "completed source bundle is missing: ${BUNDLE}" >&2
  exit 1
fi
if [[ -e "${OUTPUT}" || -e "${OUTPUT}.lock" ]]; then
  echo "audit output or lock already exists: ${OUTPUT}" >&2
  exit 1
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "${PROJECT_ROOT}"
Rscript diagnosis/tate_common_weight/audit_independent_inference_pilot.R \
  "${BUNDLE}" "${PACKAGE_LIBRARY}" "${OUTPUT}"
