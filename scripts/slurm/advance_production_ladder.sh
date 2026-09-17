#!/bin/bash
# Advance every production setting block by one checkpoint rung.
#
# HISTORY #0010. The grid has 12 setting blocks: 9 in the negative-transfer
# family (C1-C3 x K = 2/4/8) and 3 in the shared-shift family (C1 x K = 2/4/8).
# Each call re-invokes submit_rho_group_direct_tate.sh once per block with
# ROCE_SETTING scoping; the submitter truncates the batch at the next rung of
# the pre-registered ladder (1/5/10/25/50/100/200/300/400/500) and refuses a
# block whose previous checkpoint has not passed review, so running this
# repeatedly is safe and cannot outrun the gates.
#
# Usage: advance_production_ladder.sh [PRODUCTION_ROOT] [LIBRARY] [CHECK_GATE]
set -uo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
cd "${PROJECT_ROOT}"
RESULT_ROOT="${ROCE_RESULT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
PRODUCTION_ROOT="${1:-${RESULT_ROOT}/production_20260911_v2}"
PROJECT_LIBRARY="${2:-${RESULT_ROOT}/Rlib_production_20260911_v2}"
PACKAGE_CHECK_GATE="${3:-${RESULT_ROOT}/package_check_production_20260911_v2/r_cmd_check_passed.txt}"
CUTOFF_SELECTION_GATE="${ROCE_CUTOFF_SELECTION_GATE:-${RESULT_ROOT}/grouped_cutoff_pilot/cutoff_decision_n010/cutoff_selection_passed.txt}"
BATCH_SIZE="${ROCE_BATCH_SIZE:-500}"
MAX_CONCURRENT="${ROCE_MAX_CONCURRENT:-50}"
failed_blocks=0

for family in "negative_transfer:${PRODUCTION_ROOT}:C1 C2 C3" \
              "shared_shift:${PRODUCTION_ROOT}/shared_shift:C1"; do
  name="${family%%:*}"
  rest="${family#*:}"
  root="${rest%%:*}"
  configs="${rest#*:}"
  gate="${root}/rho_reuse_equivalence_final/audit/rho_reuse_equivalence_passed.txt"
  for config in ${configs}; do
    for k in 2 4 8; do
      if out="$(ROCE_SETTING="${config}:K${k}" \
        ROCE_BATCH_SIZE="${BATCH_SIZE}" \
        ROCE_MAX_CONCURRENT="${MAX_CONCURRENT}" \
        ROCE_PROJECT_LIB="${PROJECT_LIBRARY}" \
        ROCE_PACKAGE_CHECK_GATE="${PACKAGE_CHECK_GATE}" \
        ROCE_RHO_REUSE_GATE="${gate}" \
        ROCE_CUTOFF_SELECTION_GATE="${CUTOFF_SELECTION_GATE}" \
        scripts/slurm/submit_rho_group_direct_tate.sh "${root}" "${root}/raw" 2>&1)"; then
        status=0
      else
        status=$?
        failed_blocks=$((failed_blocks + 1))
      fi
      printf '[%s %s K%s] exit=%s\n%s\n' "${name}" "${config}" "${k}" "${status}" "${out:-no submitter output}"
    done
  done
done
if [[ "${failed_blocks}" -gt 0 ]]; then
  printf '[failed] %s of 12 setting blocks could not advance.\n' "${failed_blocks}" >&2
  exit 1
fi
