#!/bin/bash
# Submit the Slurm-side summarizer for diagnosis-only DR-CV scale patch outputs.

set -euo pipefail

DRY_RUN=false
DEPENDENCY=""
TIME_LIMIT="00:20:00"
MEMORY="4g"
CPUS="1"
OUTPUT_ROOT="diagnosis/c2/dr_cv_scale_patch_probe"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,40p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --dependency) DEPENDENCY="${2:?Missing dependency value}"; shift 2 ;;
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
    --export="ALL,C2_RBAL_OUTPUT_ROOT=${OUTPUT_ROOT}"
)

if [[ -n "${DEPENDENCY}" ]]; then
    SBATCH_CMD+=(--dependency="${DEPENDENCY}")
fi

SBATCH_CMD+=(diagnosis/c2/c2_dr_cv_scale_patch_summarize.cmd)

echo "======================================================"
echo " FACE-HD C2 DR-CV Scale Patch Summary Submission"
echo "======================================================"
echo " Repo:       ${REPO_DIR}"
echo " Dependency: ${DEPENDENCY:-none}"
echo " Time:       ${TIME_LIMIT}"
echo " Mem:        ${MEMORY}"
echo " Cpus:       ${CPUS}"
echo " Output root:${OUTPUT_ROOT}"
echo " Partition:  ${PARTITION}"
echo " Dry run:    ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
