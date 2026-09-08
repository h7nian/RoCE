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
# grep, not rg: batch jobs have no ripgrep on PATH, and the previous `|| true`
# turned the whole scan into a silent no-op there.  --exclude-dir prunes the
# archived diagnostic trees (retired candidate libraries and their job logs)
# and the vendored JASA template; the active sources are scanned in full.
# References to the repository directory or its GitHub slug are filtered out
# below: the checkout itself carries the old project name, which is a path, not
# a legacy package identifier.  This script is excluded because it necessarily
# contains the pattern text.
legacy_hits="$({
  grep -rnE "${legacy_pattern}" \
    DESCRIPTION NAMESPACE README.md R src inst tests scripts docs diagnosis man \
    main.R main.sh main.cmd test.sh test.cmd realdata.R realdata.sh realdata.cmd \
    data/download_rhc.sh data/download_rhc.cmd \
    --exclude-dir=JASA_template_2025 --exclude-dir=results --exclude-dir=log \
    --exclude-dir=logs --exclude-dir='Rlib*' --exclude-dir=fix_validation \
    --exclude-dir=out --exclude-dir=archived_tests \
    --exclude='*.pdf' --exclude='*.rds' --exclude='*.csv' --exclude='*.log' \
    --exclude='*.aux' --exclude='*.bbl' --exclude='*.blg' \
    --exclude='*.fls' --exclude='*.fdb_latexmk' \
    --exclude='CurrentState.md' --exclude='HISTORY.md' \
    --exclude='audit_roce_contracts.sh' || true
} | grep -vE '(/|[A-Za-z0-9_-]+/)FACE-HD([/[:space:]"'"'"'`]|$)' || true)"
if [[ -n "${legacy_hits}" ]]; then
  echo "RoCE contract audit: active legacy identifiers remain:" >&2
  echo "${legacy_hits}" >&2
  exit 1
fi

echo "[done] RoCE package, launcher, and naming contracts passed"
