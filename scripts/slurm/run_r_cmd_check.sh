#!/bin/bash
#SBATCH --job-name=roce_check
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_check.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_check.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
CHECK_ROOT="${ROCE_CHECK_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/package_check}"
PACKAGE_TARBALL="${CHECK_ROOT}/RoCE_0.1.0.tar.gz"
CHECK_DIRECTORY="${CHECK_ROOT}/RoCE.Rcheck"
SOURCE_STAGE="${CHECK_ROOT}/source/RoCE"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export _R_CHECK_FORCE_SUGGESTS_=false

"${PROJECT_ROOT}/scripts/slurm/audit_roce_contracts.sh" "${PROJECT_ROOT}"
source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
SOURCE_FINGERPRINT_BEFORE="$(
  roce_package_source_fingerprint "${PROJECT_ROOT}"
)"

if [[ -e "${PACKAGE_TARBALL}" || -e "${CHECK_DIRECTORY}" ]]; then
  echo "Refusing to mix package-check artifacts in a non-fresh directory: ${CHECK_ROOT}" >&2
  echo "Set ROCE_CHECK_ROOT to a new path so earlier audit results remain immutable." >&2
  exit 2
fi
mkdir -p "${CHECK_ROOT}" "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

# Stage exactly the package files that survive .Rbuildignore.  The repository
# also contains tens of thousands of archived diagnostic files; asking
# R CMD build to recursively enumerate those ignored files is needlessly slow
# on MSI's shared filesystem.  A fresh allow-listed stage is equivalent to the
# package source after .Rbuildignore and makes the audit reproducibly fast.
mkdir -p "${SOURCE_STAGE}"
cp "${PROJECT_ROOT}/DESCRIPTION" \
   "${PROJECT_ROOT}/LICENSE" \
   "${PROJECT_ROOT}/NAMESPACE" \
   "${PROJECT_ROOT}/README.md" \
   "${SOURCE_STAGE}/"
cp -a "${PROJECT_ROOT}/R" \
      "${PROJECT_ROOT}/inst" \
      "${PROJECT_ROOT}/man" \
      "${SOURCE_STAGE}/"
mkdir -p "${SOURCE_STAGE}/src"
cp "${PROJECT_ROOT}"/src/*.cpp \
   "${PROJECT_ROOT}"/src/*.h \
   "${PROJECT_ROOT}/src/Makevars" \
   "${SOURCE_STAGE}/src/"
mkdir -p "${SOURCE_STAGE}/tests/testthat"
cp "${PROJECT_ROOT}/tests/testthat.R" "${SOURCE_STAGE}/tests/"
cp "${PROJECT_ROOT}"/tests/testthat/*.R "${SOURCE_STAGE}/tests/testthat/"

# Build the staged exact package tree and keep every check artifact outside
# the project source. --no-manual avoids coupling this audit to the manuscript
# TeX toolchain, which has its own compilation job.
cd "${CHECK_ROOT}"
R CMD build --no-build-vignettes --no-manual "${SOURCE_STAGE}"
test -s "${PACKAGE_TARBALL}"

R CMD check --no-manual --no-build-vignettes "${PACKAGE_TARBALL}"

if ! grep -Fqx -- "Status: OK" "${CHECK_DIRECTORY}/00check.log"; then
  echo "R CMD check exited without an exact Status: OK line" >&2
  exit 1
fi
SOURCE_FINGERPRINT_AFTER="$(
  roce_package_source_fingerprint "${PROJECT_ROOT}"
)"
if [[ "${SOURCE_FINGERPRINT_AFTER}" != "${SOURCE_FINGERPRINT_BEFORE}" ]]; then
  echo "package source changed during R CMD check; no gate written" >&2
  exit 1
fi
CHECK_GATE="${CHECK_ROOT}/r_cmd_check_passed.txt"
CHECK_GATE_TEMPORARY="$(mktemp "${CHECK_ROOT}/.r_cmd_check_passed.XXXXXX")"
{
  echo "r_cmd_check=passed"
  echo "package_source_fingerprint=${SOURCE_FINGERPRINT_BEFORE}"
  echo "package_tarball_sha256=$(sha256sum "${PACKAGE_TARBALL}" | awk '{print $1}')"
  echo "check_driver_md5=$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"
} > "${CHECK_GATE_TEMPORARY}"
mv "${CHECK_GATE_TEMPORARY}" "${CHECK_GATE}"
echo "[done] R CMD check passed: ${CHECK_GATE}"
