#!/bin/bash
# Submit Slurm-side summary of C2 score/moment audit outputs.

set -euo pipefail

DRY_RUN=false
TIME_LIMIT="00:10:00"
MEMORY="2g"
CPUS="1"
OUTPUT_ROOT="diagnosis/c2/score_moment_audit"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,34p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --time) TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem) MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --cpus) CPUS="${2:?Missing cpus value}"; shift 2 ;;
        --output-root) OUTPUT_ROOT="${2:?Missing output root}"; shift 2 ;;
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
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --output="${OUTPUT_ROOT}/log/%x_%j.out"
    --error="${OUTPUT_ROOT}/log/%x_%j.err"
    --export="ALL,C2_SCORE_AUDIT_OUTPUT_ROOT=${OUTPUT_ROOT}"
    diagnosis/c2/c2_score_moment_summarize.cmd
)

echo "======================================================"
echo " FACE-HD C2 Score/Moment Summary Submission"
echo "======================================================"
echo " Repo:        ${REPO_DIR}"
echo " Output root: ${OUTPUT_ROOT}"
echo " Time:        ${TIME_LIMIT}"
echo " Mem:         ${MEMORY}"
echo " Cpus:        ${CPUS}"
echo " Partition:   ${PARTITION}"
echo " Dry run:     ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
