#!/bin/bash
#SBATCH --job-name=roce_nlambda_manifest
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_nlambda_manifest.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_nlambda_manifest.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
OUTPUT_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
GRID_A="${2:-20}"
GRID_B="${3:-50}"
N_SEEDS="${ROCE_NLAMBDA_VALIDATION_SEEDS:-5}"

if [[ ! "${N_SEEDS}" =~ ^[1-9][0-9]*$ ]] || [[ "${N_SEEDS}" -gt 10 ]]; then
  echo "ROCE_NLAMBDA_VALIDATION_SEEDS must be an integer from 1 through 10" >&2
  exit 1
fi
for grid in "${GRID_A}" "${GRID_B}"; do
  if [[ ! "${grid}" =~ ^[1-9][0-9]*$ ]] || [[ "${grid}" -lt 2 ]]; then
    echo "nuisance-grid sizes must be integers of at least 2: ${grid}" >&2
    exit 1
  fi
done
if [[ "${GRID_A}" -ge "${GRID_B}" ]]; then
  echo "GRID_A must be smaller than GRID_B: ${GRID_A} >= ${GRID_B}" >&2
  exit 1
fi

MANIFEST_A="${OUTPUT_ROOT}/manifest_nlambda${GRID_A}_validation.csv"
MANIFEST_B="${OUTPUT_ROOT}/manifest_nlambda${GRID_B}_validation.csv"
if [[ -e "${MANIFEST_A}" || -e "${MANIFEST_B}" ]]; then
  echo "refusing to overwrite an existing nlambda validation manifest" >&2
  exit 2
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
STAGING_DIRECTORY="$(mktemp -d "${OUTPUT_ROOT}/.nlambda_manifest.XXXXXX")"

Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${STAGING_DIRECTORY}/manifest_nlambda${GRID_A}_validation.csv" \
  nlambda_validation "${N_SEEDS}" "${GRID_A}"
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${STAGING_DIRECTORY}/manifest_nlambda${GRID_B}_validation.csv" \
  nlambda_validation "${N_SEEDS}" "${GRID_B}"

mv "${STAGING_DIRECTORY}/manifest_nlambda${GRID_A}_validation.csv" \
  "${MANIFEST_A}"
mv "${STAGING_DIRECTORY}/manifest_nlambda${GRID_B}_validation.csv" \
  "${MANIFEST_B}"
rmdir "${STAGING_DIRECTORY}"

echo "wrote paired nlambda=${GRID_A}:${GRID_B} validation manifests with ${N_SEEDS} seeds each"
