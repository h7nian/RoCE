#!/bin/bash
#SBATCH --job-name=roce_rhc_cutoff
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rhc_cutoff.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rhc_cutoff.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/rhc_supported_primary}"
INPUT_PATH="${ROCE_RHC_RESULT:-${OUTPUT_ROOT}/rhc_direct_tate.rds}"
CUTOFFS="${ROCE_RHC_CUTOFFS:-1,1.5,2,2.5,3}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
RHC_SENSITIVITY_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/summarize_rhc_cutoff_sensitivity.sh" \
  "${PROJECT_ROOT}/scripts/slurm/summarize_rhc_cutoff_sensitivity.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_RHC_SENSITIVITY_WORKFLOW_FINGERPRINT="${RHC_SENSITIVITY_WORKFLOW_FINGERPRINT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/summarize_rhc_cutoff_sensitivity.R \
  "${INPUT_PATH}" "${OUTPUT_ROOT}" "${CUTOFFS}"
