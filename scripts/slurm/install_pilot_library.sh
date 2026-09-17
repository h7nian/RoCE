#!/bin/bash
# Install RoCE from an immutable git-archive stage into a fresh library, so the
# installed objects are compiled with R's own flags (-O2 -g0 on MSI) and never
# reuse working-tree objects left behind by devtools::load_all (-O0 -g).
#
# usage: install_pilot_library.sh LIBRARY_DIR [GIT_REF]   (default GIT_REF = HEAD)
set -euo pipefail

LIBRARY_DIR="$1"
GIT_REF="${2:-HEAD}"
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
cd "${PROJECT_ROOT}"

if [[ -e "${LIBRARY_DIR}" ]]; then
  echo "Refusing to reuse an existing library directory: ${LIBRARY_DIR}" >&2
  exit 2
fi
STAGE="$(mktemp -d "${TMPDIR:-/tmp}/roce_stage_XXXXXX")"
trap 'rm -rf "${STAGE}"' EXIT
git archive --format=tar "${GIT_REF}" | tar -x -C "${STAGE}"
COMMIT="$(git rev-parse "${GIT_REF}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || true
# RcppExports.* are generated files (git-ignored); the stage needs them before install.
Rscript -e "Rcpp::compileAttributes('${STAGE}', verbose = FALSE)" > /dev/null
mkdir -p "${LIBRARY_DIR}"
R CMD INSTALL --no-docs --no-multiarch --library="${LIBRARY_DIR}" "${STAGE}" > "${LIBRARY_DIR}.install.log" 2>&1

SHARED_OBJECT="${LIBRARY_DIR}/RoCE/libs/RoCE.so"
if [[ ! -f "${SHARED_OBJECT}" ]]; then
  echo "Install failed; see ${LIBRARY_DIR}.install.log" >&2
  exit 1
fi
if grep -q -- "-O0" "${LIBRARY_DIR}.install.log" || readelf --debug-dump=info "${SHARED_OBJECT}" 2>/dev/null | grep -q -- "-O0"; then
  echo "Installed objects were compiled with -O0; refusing to publish ${LIBRARY_DIR}" >&2
  exit 1
fi
echo "${COMMIT}" > "${LIBRARY_DIR}/ROCE_SOURCE_COMMIT"
echo "installed ${LIBRARY_DIR} from ${GIT_REF} (${COMMIT:0:8}); compile flags: $(grep -m1 -o -- '-O[0-9] -g[0-9]' "${LIBRARY_DIR}.install.log" || echo 'not shown')"
