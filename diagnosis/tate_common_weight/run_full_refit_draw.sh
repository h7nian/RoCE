#!/bin/bash
#SBATCH --job-name=roce_nuisance_refit
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=5
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_nuisance_refit.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_nuisance_refit.err
set -euo pipefail

if [[ "$#" -ne 4 ]]; then
  echo "usage: run_full_refit_draw.sh BUNDLE_DIR RHO DRAW_ID OUTPUT_DIR" >&2
  exit 1
fi
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"
CHECK_GATE="${ROCE_PACKAGE_CHECK_GATE:?ROCE_PACKAGE_CHECK_GATE is required}"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
SOURCE_FINGERPRINT="$(roce_package_source_fingerprint "${PROJECT_ROOT}")"
TEST_GATE="${PROJECT_LIBRARY}/audit_tests_passed.txt"

require_refit_gate() {
  local gate="$1" key="$2" expected="$3" observed
  observed="$(awk -F= -v key="${key}" '$1 == key {print substr($0, length($1) + 2)}' "${gate}")"
  if [[ "${observed}" != "${expected}" ]]; then
    echo "full-refit gate mismatch: ${gate}, field ${key}" >&2
    exit 1
  fi
}
require_refit_gate "${TEST_GATE}" package_tests passed
require_refit_gate "${TEST_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_refit_gate "${TEST_GATE}" package_source_fingerprint "${SOURCE_FINGERPRINT}"
require_refit_gate "${CHECK_GATE}" r_cmd_check passed
require_refit_gate "${CHECK_GATE}" package_source_fingerprint "${SOURCE_FINGERPRINT}"
require_refit_gate "${TEST_GATE}" test_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
require_refit_gate "${CHECK_GATE}" check_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
# The R driver pins the grouped-CV v19 package and records the separate v18
# source-fit provenance. This cannot silently resume the archived row-CV run.
export ROCE_REFIT_CORES=1
export ROCE_NUISANCE_CV_THREADS=5
if [[ ! "${SLURM_CPUS_PER_TASK:-}" =~ ^[0-9]+$ ]] || (( SLURM_CPUS_PER_TASK < 5 )); then
  echo "full-refit requires a Slurm allocation of at least 5 CPUs" >&2
  exit 1
fi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd "${PROJECT_ROOT}"
Rscript diagnosis/tate_common_weight/run_full_refit_draw.R "$@"
