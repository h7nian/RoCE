#!/bin/bash
#SBATCH --job-name=roce_cutoff_select
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_cutoff_select.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_cutoff_select.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
SUMMARY_DIRECTORY="${1:?SUMMARY_DIRECTORY is required}"
OUTPUT_DIRECTORY="${2:?OUTPUT_DIRECTORY is required}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
CUTOFF_SELECTION_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/select_grouped_cutoff.sh" \
  "${PROJECT_ROOT}/scripts/slurm/select_grouped_cutoff.R" \
  "${PROJECT_ROOT}/scripts/slurm/atomic_output.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export ROCE_CUTOFF_SELECTION_WORKFLOW_FINGERPRINT="${CUTOFF_SELECTION_WORKFLOW_FINGERPRINT}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
echo "[selection workflow] fingerprint=${CUTOFF_SELECTION_WORKFLOW_FINGERPRINT}"
Rscript scripts/slurm/select_grouped_cutoff.R \
  "${SUMMARY_DIRECTORY}" "${OUTPUT_DIRECTORY}"
