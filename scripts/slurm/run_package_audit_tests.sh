#!/bin/bash
#SBATCH --job-name=roce_tests
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_tests.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_tests.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
AUDIT_LIBRARY="${ROCE_AUDIT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_final_audit}"
AUDIT_SOURCE="${ROCE_AUDIT_SOURCE:-${AUDIT_LIBRARY}_source/RoCE}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
"${PROJECT_ROOT}/scripts/slurm/audit_roce_contracts.sh" "${PROJECT_ROOT}"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
SOURCE_FINGERPRINT_BEFORE="$(
  roce_package_source_fingerprint "${PROJECT_ROOT}"
)"
TEST_SUITE_FINGERPRINT_BEFORE="$({
  find "${PROJECT_ROOT}/tests" -type f -name '*.R' -print0 \
    | sort -z | xargs -0 sha256sum
} | sha256sum | awk '{print $1}')"
if [[ -e "${AUDIT_SOURCE}" ]]; then
  echo "Refusing to reuse an existing package-test source stage: ${AUDIT_SOURCE}" >&2
  exit 2
fi
mkdir -p "${AUDIT_LIBRARY}" "${AUDIT_SOURCE}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

# Install from an immutable allow-listed stage. R CMD INSTALL --preclean may
# remove and recreate src/*.o in its input tree, so compiling directly from the
# shared checkout can race with R CMD check or another audit job.
cp "${PROJECT_ROOT}/DESCRIPTION" \
   "${PROJECT_ROOT}/LICENSE" \
   "${PROJECT_ROOT}/NAMESPACE" \
   "${PROJECT_ROOT}/README.md" \
   "${AUDIT_SOURCE}/"
cp -a "${PROJECT_ROOT}/R" \
      "${PROJECT_ROOT}/inst" \
      "${PROJECT_ROOT}/man" \
      "${AUDIT_SOURCE}/"
mkdir -p "${AUDIT_SOURCE}/src"
cp "${PROJECT_ROOT}"/src/*.cpp \
   "${PROJECT_ROOT}"/src/*.h \
   "${PROJECT_ROOT}/src/Makevars" \
   "${AUDIT_SOURCE}/src/"

R CMD INSTALL --preclean --clean --library="${AUDIT_LIBRARY}" "${AUDIT_SOURCE}"
ROCE_AUDIT_LIB="${AUDIT_LIBRARY}" ROCE_TEST_INSTALLED=1 \
  ROCE_C2_CV_SCALE_AUDIT=1 NOT_CRAN=true Rscript -e '
  .libPaths(c(Sys.getenv("ROCE_AUDIT_LIB"), .libPaths()))
  library(RoCE)
  slurm_r_scripts <- list.files(
    "scripts/slurm", pattern = "[.]R$", full.names = TRUE
  )
  root_r_launchers <- c("main.R", "realdata.R")
  invisible(lapply(c(root_r_launchers, slurm_r_scripts), parse))
  testthat::test_dir(
    "tests/testthat",
    reporter = "summary",
    stop_on_failure = TRUE,
    stop_on_warning = FALSE
  )
'

PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${AUDIT_LIBRARY}")"
SOURCE_FINGERPRINT_AFTER="$(
  roce_package_source_fingerprint "${PROJECT_ROOT}"
)"
TEST_SUITE_FINGERPRINT_AFTER="$({
  find "${PROJECT_ROOT}/tests" -type f -name '*.R' -print0 \
    | sort -z | xargs -0 sha256sum
} | sha256sum | awk '{print $1}')"
if [[ "${SOURCE_FINGERPRINT_AFTER}" != "${SOURCE_FINGERPRINT_BEFORE}" ||
      "${TEST_SUITE_FINGERPRINT_AFTER}" != "${TEST_SUITE_FINGERPRINT_BEFORE}" ]]; then
  echo "package source or tests changed during the immutable audit; no gate written" >&2
  exit 1
fi
AUDIT_GATE="${AUDIT_LIBRARY}/audit_tests_passed.txt"
AUDIT_GATE_TEMPORARY="$(mktemp "${AUDIT_LIBRARY}/.audit_tests_passed.XXXXXX")"
{
  echo "package_tests=passed"
  echo "package_fingerprint=${PACKAGE_FINGERPRINT}"
  echo "package_source_fingerprint=${SOURCE_FINGERPRINT_BEFORE}"
  echo "test_suite_fingerprint=${TEST_SUITE_FINGERPRINT_BEFORE}"
  echo "test_driver_md5=$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
} > "${AUDIT_GATE_TEMPORARY}"
mv "${AUDIT_GATE_TEMPORARY}" "${AUDIT_GATE}"
echo "[done] immutable-source package tests passed: ${AUDIT_GATE}"
