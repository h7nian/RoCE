#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${1:-${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}}"
cd "${PROJECT_ROOT}"

required_files=(
  DESCRIPTION NAMESPACE R/RcppExports.R src/RcppExports.cpp
  inst/include/RoCE_types.h main.R main.sh main.cmd test.sh test.cmd
  realdata.R realdata.sh realdata.cmd
)
for required_file in "${required_files[@]}"; do
  if [[ ! -f "${required_file}" ]]; then
    echo "RoCE contract audit: missing ${required_file}" >&2
    exit 1
  fi
done

grep -Fqx 'Package: RoCE' DESCRIPTION
grep -Fqx 'useDynLib(RoCE, .registration = TRUE)' NAMESPACE
grep -Fq 'R_init_RoCE' src/RcppExports.cpp
grep -Fq '_RoCE_fit_unified_density_ratio_cpp' R/RcppExports.R
grep -Fq 'ROCE_TEST_INSTALLED' test.cmd
grep -Fq 'choices = c("roce", "face")' main.R

for launcher in \
  main.sh main.cmd test.sh test.cmd realdata.sh realdata.cmd \
  data/download_rhc.sh data/download_rhc.cmd; do
  bash -n "${launcher}"
done

legacy_pattern='FACE''HD|FACE''-HD|face''hd|face''_hd'
legacy_hits="$({
  rg -n "${legacy_pattern}" \
    DESCRIPTION NAMESPACE README.md R src inst tests scripts docs diagnosis man \
    main.R main.sh main.cmd test.sh test.cmd realdata.R realdata.sh realdata.cmd \
    data/download_rhc.sh data/download_rhc.cmd \
    --glob '!docs/JASA_template_2025/**' \
    --glob '!**/results/**' --glob '!**/logs/**' --glob '!**/Rlib*/**' \
    --glob '!*.pdf' --glob '!*.rds' --glob '!*.csv' --glob '!*.log' \
    --glob '!*.aux' --glob '!*.bbl' --glob '!*.blg' \
    --glob '!*.fls' --glob '!*.fdb_latexmk' || true
} )"
if [[ -n "${legacy_hits}" ]]; then
  echo "RoCE contract audit: active legacy identifiers remain:" >&2
  echo "${legacy_hits}" >&2
  exit 1
fi

echo "[done] RoCE package, launcher, and naming contracts passed"
