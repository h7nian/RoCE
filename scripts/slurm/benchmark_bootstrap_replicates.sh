#!/bin/bash
#SBATCH --job-name=roce_bootbench
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_bootbench.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_bootbench.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
OUTPUT_ROOT="${ROCE_BOOTSTRAP_BENCHMARK_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/bootstrap_benchmark}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
BOOTSTRAP_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/benchmark_bootstrap_replicates.sh" \
  "${PROJECT_ROOT}/scripts/slurm/benchmark_bootstrap_replicates.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
ROCE_PROJECT_LIB="${PROJECT_LIBRARY}" \
  ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}" \
  ROCE_BOOTSTRAP_WORKFLOW_FINGERPRINT="${BOOTSTRAP_WORKFLOW_FINGERPRINT}" \
  Rscript scripts/slurm/benchmark_bootstrap_replicates.R \
  "${OUTPUT_ROOT}" 10
