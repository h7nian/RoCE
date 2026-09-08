#!/bin/bash
#SBATCH --job-name=roce_nlambda100_manifest
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_nlambda100_manifest.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_nlambda100_manifest.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
MANIFEST_PATH="${RESULT_ROOT}/manifest_nlambda100_validation.csv"

if [[ -e "${MANIFEST_PATH}" ]]; then
  echo "refusing to overwrite existing validation manifest: ${MANIFEST_PATH}" >&2
  exit 2
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}" "${RESULT_ROOT}/logs"
STAGING_DIRECTORY="$(mktemp -d "${RESULT_ROOT}/.nlambda100_manifest.XXXXXX")"
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${STAGING_DIRECTORY}/manifest.csv" nlambda_validation 5 100
mv "${STAGING_DIRECTORY}/manifest.csv" "${MANIFEST_PATH}"
rmdir "${STAGING_DIRECTORY}"

echo "wrote five-seed nlambda=100 validation manifest; seed 1 is the spot check"
