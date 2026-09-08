#!/bin/bash
#SBATCH --job-name=roce_inference_checkpoint
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=40
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%A_%a_independent_inference.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%A_%a_independent_inference.err
set -euo pipefail

# Scheduling adapter only. The scientific runner and its fingerprint stay fixed.
if [[ "$#" -ne 3 ]]; then
  echo "usage: run_independent_inference_pilot_array.sh MANIFEST OUTPUT_ROOT PREVIOUS_CHECKPOINT" >&2
  exit 1
fi
if [[ ! "${SLURM_ARRAY_TASK_ID:-}" =~ ^([1-9]|[1-9][0-9]|100)$ ]]; then
  echo "SLURM_ARRAY_TASK_ID must be one integer in [1,100]" >&2
  exit 1
fi
TASK_ID="${SLURM_ARRAY_TASK_ID}"
PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST="$(readlink -f "$1")"
OUTPUT_ROOT="$(readlink -f "$2")"
CHECKPOINT="$(readlink -f "$3")"
if [[ ! -f "${MANIFEST}" || ! -d "${OUTPUT_ROOT}" || ! -d "${CHECKPOINT}" ]] ||
   [[ "${MANIFEST}" != "${OUTPUT_ROOT}/manifest.csv" ]]; then
  echo "manifest, existing output root or checkpoint path is invalid" >&2
  exit 1
fi
(cd "${CHECKPOINT}" && sha256sum -c sha256.txt)
METADATA="${CHECKPOINT}/metadata.txt"
require_checkpoint_field() {
  local key="$1" expected="$2" observed
  observed="$(awk -F= -v key="${key}" '$1 == key {print substr($0,length($1)+2)}' "${METADATA}")"
  [[ "${observed}" == "${expected}" ]] || {
    echo "checkpoint field mismatch: ${key}" >&2; exit 1;
  }
}
require_checkpoint_field independent_inference_summary complete
require_checkpoint_field inference_validated FALSE
require_checkpoint_field statistical_review_required TRUE
require_checkpoint_field manifest_fingerprint "$(sha256sum "${MANIFEST}" | awk '{print $1}')"
require_checkpoint_field summary_script_fingerprint \
  "$(sha256sum "${PROJECT_ROOT}/diagnosis/tate_common_weight/summarize_independent_inference_pilot.R" | awk '{print $1}')"
PREVIOUS_N="$(awk -F= '$1 == "n" {print $2}' "${METADATA}")"
case "${PREVIOUS_N}" in
  1) NEXT_N=5 ;;
  5) NEXT_N=10 ;;
  10) NEXT_N=25 ;;
  25) NEXT_N=50 ;;
  50) NEXT_N=100 ;;
  *) echo "checkpoint does not permit a further declared batch" >&2; exit 1 ;;
esac
printf -v CHECKPOINT_LABEL 'n%03d' "${PREVIOUS_N}"
if [[ "${CHECKPOINT}" != "${OUTPUT_ROOT}/summaries/${CHECKPOINT_LABEL}" ]] ||
   (( TASK_ID <= PREVIOUS_N || TASK_ID > NEXT_N )); then
  echo "task or checkpoint path is outside the next declared prefix batch" >&2
  exit 1
fi
SIM_ID=$((10000 + TASK_ID))
printf -v SEED_LABEL 'seed_%06d' "${SIM_ID}"
DESTINATION="${OUTPUT_ROOT}/${SEED_LABEL}"
if [[ -e "${DESTINATION}" || -e "${DESTINATION}.lock" ]]; then
  echo "task output or lock already exists: ${DESTINATION}" >&2
  exit 1
fi
exec bash "${PROJECT_ROOT}/diagnosis/tate_common_weight/run_independent_inference_pilot.sh" \
  "${MANIFEST}" "${TASK_ID}" "${DESTINATION}"
