#!/bin/bash
#SBATCH --job-name=roce_smoke_audit
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_smoke_audit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_smoke_audit.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
OUTPUT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/smoke_final_direct_tate_raw}"
SENSITIVITY_OUTPUT_ROOT="${ROCE_SENSITIVITY_OUTPUT_ROOT:-${OUTPUT_ROOT%/}_reused_sensitivity}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_EXPECT_NUISANCE_CV_THREADS="${ROCE_EXPECT_NUISANCE_CV_THREADS:-5}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_SENSITIVITY_OUTPUT_ROOT="${SENSITIVITY_OUTPUT_ROOT}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/audit_direct_tate_smoke.R "${OUTPUT_ROOT}"
