#!/bin/bash
# Prelaunch integration checks against an audited checkpoint; no worker runs.
# Usage: test_independent_inference_pilot_array.sh [PILOT_ROOT [PREVIOUS_N]]
# Run before the next batch: existing output/lock paths intentionally fail.
set -euo pipefail
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PILOT_ROOT="$(readlink -f "${1:-results/direct_tate_mc500_b5000/independent_inference_pilot_v19}")"
PREVIOUS_N="${2:-1}"
case "${PREVIOUS_N}" in
  1) NEXT_N=5 ;;
  5) NEXT_N=10 ;;
  10) NEXT_N=25 ;;
  25) NEXT_N=50 ;;
  50) NEXT_N=100 ;;
  *) echo "PREVIOUS_N must be one of 1,5,10,25,50" >&2; exit 1 ;;
esac
printf -v CHECKPOINT_LABEL 'n%03d' "${PREVIOUS_N}"
CHECKPOINT="${PILOT_ROOT}/summaries/${CHECKPOINT_LABEL}"
ADAPTER="${PROJECT_ROOT}/diagnosis/tate_common_weight/run_independent_inference_pilot_array.sh"
STUB_DIRECTORY="$(mktemp -d /tmp/roce_inference_array_test.XXXXXX)"
ln -s /bin/echo "${STUB_DIRECTORY}/bash"
cleanup_stub() {
  unlink "${STUB_DIRECTORY}/bash"
  rmdir "${STUB_DIRECTORY}"
}
trap cleanup_stub EXIT

/bin/bash -n "${ADAPTER}"
for task_id in "$((PREVIOUS_N + 1))" "${NEXT_N}"; do
  printf -v seed_label 'seed_%06d' "$((10000 + task_id))"
  output="$(env PATH="${STUB_DIRECTORY}:${PATH}" SLURM_ARRAY_TASK_ID="${task_id}" \
    /bin/bash "${ADAPTER}" "${PILOT_ROOT}/manifest.csv" "${PILOT_ROOT}" \
    "${CHECKPOINT}")"
  expected="${PROJECT_ROOT}/diagnosis/tate_common_weight/run_independent_inference_pilot.sh ${PILOT_ROOT}/manifest.csv ${task_id} ${PILOT_ROOT}/${seed_label}"
  [[ "$(tail -n 1 <<< "${output}")" == "${expected}" ]] || {
    echo "stubbed worker arguments did not match" >&2; exit 1;
  }
done
for task_id in 0 "${PREVIOUS_N}" "$((NEXT_N + 1))" -1 invalid; do
  if env PATH="${STUB_DIRECTORY}:${PATH}" SLURM_ARRAY_TASK_ID="${task_id}" \
      /bin/bash "${ADAPTER}" "${PILOT_ROOT}/manifest.csv" "${PILOT_ROOT}" \
      "${CHECKPOINT}" >/dev/null 2>&1; then
    echo "out-of-batch task was accepted: ${task_id}" >&2; exit 1
  fi
done
echo "array adapter checks passed for ${CHECKPOINT_LABEL}: 2 stubbed positive mappings, 5 rejected IDs; no fitting"
