#!/bin/bash
#SBATCH --job-name=roce_prepare_group_cutoff
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_prepare_group_cutoff.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_prepare_group_cutoff.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
DIAGNOSTIC_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot}"
CONFIG="${2:-C1}"
SOURCE_COUNT="${3:-2}"
N_SIMULATIONS="${4:-10}"
RHO_VALUES="${5:-0;1;2.5}"
CUTOFF_VALUES="${6:-1;1.5;2;2.5;3}"
SIMULATION_START="${7:-501}"
MANIFEST_PATH="${DIAGNOSTIC_ROOT}/manifest.csv"
RAW_DIRECTORY="${DIAGNOSTIC_ROOT}/raw"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${DIAGNOSTIC_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

Rscript scripts/slurm/build_grouped_cutoff_diagnostic_manifest.R \
  "${MANIFEST_PATH}" "${CONFIG}" "${SOURCE_COUNT}" "${N_SIMULATIONS}" \
  "${RHO_VALUES}" "${CUTOFF_VALUES}" "${SIMULATION_START}"

ROCE_SUBMIT_DRY_RUN=1 \
ROCE_BATCH_SIZE=1 \
ROCE_MAX_CONCURRENT=1 \
  bash scripts/slurm/submit_grouped_cutoff_diagnostic.sh \
    "${MANIFEST_PATH}" "${RAW_DIRECTORY}"

sha256sum "${MANIFEST_PATH}" > "${DIAGNOSTIC_ROOT}/manifest.sha256"
echo "prepared and dry-run audited ${MANIFEST_PATH}"
