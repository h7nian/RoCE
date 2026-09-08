#!/bin/bash
#SBATCH --job-name=roce_rho_group
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=40
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a_rho_group.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a_rho_group.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
GROUP_MANIFEST="${ROCE_GROUP_MANIFEST:?ROCE_GROUP_MANIFEST is required}"
PRIMARY_MANIFEST="${ROCE_PRIMARY_MANIFEST:?ROCE_PRIMARY_MANIFEST is required}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw}"
SENSITIVITY_OUTPUT_ROOT="${ROCE_SENSITIVITY_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/reused_sensitivity_raw}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

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
export ROCE_SENSITIVITY_OUTPUT_ROOT="${SENSITIVITY_OUTPUT_ROOT}"
export ROCE_PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
export ROCE_WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
export ROCE_MANIFEST_FINGERPRINT="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
export ROCE_GROUP_MANIFEST_FINGERPRINT="$(sha256sum "${GROUP_MANIFEST}" | awk '{print $1}')"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${SENSITIVITY_OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

Rscript scripts/slurm/run_direct_tate_rho_group_task.R \
  "${GROUP_MANIFEST}" "${PRIMARY_MANIFEST}" "${SLURM_ARRAY_TASK_ID}"
