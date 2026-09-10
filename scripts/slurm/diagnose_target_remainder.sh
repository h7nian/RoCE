#!/bin/bash
#SBATCH --job-name=roce_remainder
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_remainder.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_remainder.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/target_remainder}"
N_REPLICATIONS="${ROCE_REMAINDER_REPLICATIONS:-500}"
# Configuration under diagnosis: C(1) isolates regularization bias, C(3) adds
# propensity misspecification (HISTORY #0011).
CONFIG="${ROCE_REMAINDER_CONFIG:-C3}"
NUISANCE_LAMBDA_RULE="${ROCE_NUISANCE_LAMBDA_RULE:-min}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
REMAINDER_DIAGNOSTIC_WORKFLOW_FINGERPRINT="$(roce_files_fingerprint \
  "${PROJECT_ROOT}/scripts/slurm/diagnose_target_remainder.sh" \
  "${PROJECT_ROOT}/scripts/slurm/diagnose_target_remainder.R" \
  "${PROJECT_ROOT}/scripts/slurm/atomic_output.R" \
  "${PROJECT_ROOT}/scripts/slurm/direct_tate_task_helpers.R" \
  "${PROJECT_ROOT}/scripts/slurm/result_provenance.R" \
  "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_REMAINDER_DIAGNOSTIC_WORKFLOW_FINGERPRINT="${REMAINDER_DIAGNOSTIC_WORKFLOW_FINGERPRINT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
echo "[package] library=${PROJECT_LIBRARY} fingerprint=${PACKAGE_FINGERPRINT}"
echo "[workflow] fingerprint=${REMAINDER_DIAGNOSTIC_WORKFLOW_FINGERPRINT}"
Rscript scripts/slurm/diagnose_target_remainder.R \
  "${N_REPLICATIONS}" "${OUTPUT_ROOT}" "${NUISANCE_LAMBDA_RULE}" "${CONFIG}"
