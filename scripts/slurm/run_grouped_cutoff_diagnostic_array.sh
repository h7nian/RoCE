#!/bin/bash
#SBATCH --job-name=roce_group_cutoff
#SBATCH --time=06:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=20
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a_group_cutoff.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a_group_cutoff.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_PATH="${ROCE_MANIFEST:?ROCE_MANIFEST is required}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:?ROCE_OUTPUT_ROOT is required}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
roce_resolve_nuisance_cv_threads 5 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"
roce_require_nuisance_cv_cpu_capacity \
  "${SLURM_CPUS_PER_TASK}" "${NUISANCE_CV_THREADS}"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}"
export ROCE_PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
export ROCE_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/run_grouped_cutoff_diagnostic_array.sh" \
  "${PROJECT_ROOT}/scripts/slurm/run_grouped_cutoff_diagnostic_task.R" \
  "${PROJECT_ROOT}/scripts/slurm/grouped_cutoff_diagnostic_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/resource_topology.R" \
  "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"
export ROCE_MANIFEST_FINGERPRINT="$(sha256sum "${MANIFEST_PATH}" | awk '{print $1}')"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/run_grouped_cutoff_diagnostic_task.R \
  "${MANIFEST_PATH}" "${SLURM_ARRAY_TASK_ID}"
