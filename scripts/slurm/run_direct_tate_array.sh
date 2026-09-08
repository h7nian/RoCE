#!/bin/bash
#SBATCH --job-name=roce_tate
# A full p=100 replicate fits both treatment arms plus all comparison methods.
# Empirical MSI smoke runs exceeded three hours, so leave enough headroom to
# avoid discarding a nearly complete atomic task output.
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_PATH="${ROCE_MANIFEST:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/manifest.csv}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw}"
SENSITIVITY_OUTPUT_ROOT="${ROCE_SENSITIVITY_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/reused_sensitivity_raw}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
roce_resolve_nuisance_cv_threads 1 5
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"
roce_require_nuisance_cv_cpu_capacity \
  "${SLURM_CPUS_PER_TASK}" "${NUISANCE_CV_THREADS}"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
MANIFEST_FINGERPRINT="$(sha256sum "${MANIFEST_PATH}" | awk '{print $1}')"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}"
export ROCE_SENSITIVITY_OUTPUT_ROOT="${SENSITIVITY_OUTPUT_ROOT}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_MANIFEST_FINGERPRINT="${MANIFEST_FINGERPRINT}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${SENSITIVITY_OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

echo "[package] library=${PROJECT_LIBRARY} fingerprint=${PACKAGE_FINGERPRINT}"
echo "[workflow] fingerprint=${WORKFLOW_FINGERPRINT} manifest=${MANIFEST_FINGERPRINT}"
echo "[resource] cpus=${SLURM_CPUS_PER_TASK} nuisance_cv_threads=${ROCE_NUISANCE_CV_THREADS}"

Rscript scripts/slurm/run_direct_tate_task.R \
  "${MANIFEST_PATH}" "${SLURM_ARRAY_TASK_ID}"
