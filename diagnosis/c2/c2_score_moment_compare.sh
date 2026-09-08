#!/bin/bash
# Submit cached-vs-candidate C2 score/moment summary comparison.

set -euo pipefail

DRY_RUN=false
BASE_ROOT="diagnosis/c2/score_moment_audit"
CANDIDATE_ROOT="diagnosis/c2/score_moment_audit_nocache"
OUTPUT_ROOT="diagnosis/c2/score_moment_compare"
TIME_LIMIT="00:10:00"
MEMORY="2g"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,36p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --base-root) BASE_ROOT="${2:?Missing base root}"; shift 2 ;;
        --candidate-root) CANDIDATE_ROOT="${2:?Missing candidate root}"; shift 2 ;;
        --output-root) OUTPUT_ROOT="${2:?Missing output root}"; shift 2 ;;
        --time) TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem) MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p "${OUTPUT_ROOT}/log" "${OUTPUT_ROOT}/summary"

SBATCH_CMD=(
    sbatch
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --partition="${PARTITION}"
    --output="${OUTPUT_ROOT}/log/%x_%j.out"
    --error="${OUTPUT_ROOT}/log/%x_%j.err"
    --export="ALL,C2_SCORE_COMPARE_BASE_ROOT=${BASE_ROOT},C2_SCORE_COMPARE_CANDIDATE_ROOT=${CANDIDATE_ROOT},C2_SCORE_COMPARE_OUTPUT_ROOT=${OUTPUT_ROOT}"
    diagnosis/c2/c2_score_moment_compare.cmd
)

echo "======================================================"
echo " RoCE C2 Score/Moment Compare Submission"
echo "======================================================"
echo " Repo:       ${REPO_DIR}"
echo " Base:       ${BASE_ROOT}"
echo " Candidate:  ${CANDIDATE_ROOT}"
echo " Output:     ${OUTPUT_ROOT}"
echo " Dry run:    ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
