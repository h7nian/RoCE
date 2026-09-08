#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/fitting_radius_diagnostic/raw}"
CPUS_PER_TASK="${3:-8}"
MANIFEST_PATH="${MANIFEST_ROOT}/manifest_truncation_diagnostic.csv"

if [[ ! -f "${MANIFEST_PATH}" ]]; then
  echo "fitting-radius manifest not found: ${MANIFEST_PATH}" >&2
  exit 1
fi

"${PROJECT_ROOT}/scripts/slurm/submit_direct_tate.sh" \
  "${MANIFEST_PATH}" "${OUTPUT_ROOT}" "${CPUS_PER_TASK}" 8G
