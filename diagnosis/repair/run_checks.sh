#!/bin/bash
# Run from an immutable stage created by stage_source.py. Build products and
# logs stay beside that stage, not in the working repository or home cache.
set -euo pipefail
CHECK_ROOT="${1:?Pass the staged check directory}"
TEST_FILTER="${2:-}"
CHECK_PHASE="${3:-all}"
case "${CHECK_PHASE}" in
  all|build|test) ;;
  *) echo "Phase must be all, build, or test" >&2; exit 2 ;;
esac
case "${CHECK_ROOT}" in
  /scratch.global/zhan9381/FACE-HD/*) ;;
  *) echo "Check directory must be under the FACE-HD scratch root" >&2; exit 2 ;;
esac
test -f "${CHECK_ROOT}/source/DESCRIPTION"
test -f "${CHECK_ROOT}/source_manifest.json"
mkdir -p "${CHECK_ROOT}/tmp"
export TMPDIR="${CHECK_ROOT}/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_NUISANCE_CV_THREADS=1 ROCE_TEST_INSTALLED=1 NOT_CRAN=true
export ROCE_CHECK_ROOT="${CHECK_ROOT}" ROCE_TEST_FILTER="${TEST_FILTER}"
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
cd "${CHECK_ROOT}/source"
if [[ "${CHECK_PHASE}" != test ]]; then
  python3 - "${CHECK_ROOT}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys
root = Path(sys.argv[1])
for entry in json.loads((root / "source_manifest.json").read_text()):
    path = root / "source" / entry["path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
        sys.exit("Staged source changed before compilation: " + entry["path"])
PY
  test ! -e "${CHECK_ROOT}/Rlib"
  mkdir "${CHECK_ROOT}/Rlib"
  Rscript --vanilla -e 'Rcpp::compileAttributes(".", verbose=FALSE)' > "${CHECK_ROOT}/attributes.log" 2>&1
  R CMD INSTALL --preclean --no-multiarch --library="${CHECK_ROOT}/Rlib" . > "${CHECK_ROOT}/install.log" 2>&1
  touch "${CHECK_ROOT}/INSTALLED"
fi
if [[ "${CHECK_PHASE}" != build ]]; then
  test -f "${CHECK_ROOT}/INSTALLED"
  export R_LIBS="${CHECK_ROOT}/Rlib:${R_LIBS_USER}"
  Rscript --vanilla diagnosis/repair/run_tests.R > "${CHECK_ROOT}/tests.log" 2>&1
  touch "${CHECK_ROOT}/TESTS_PASSED"
fi
