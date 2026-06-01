#!/bin/bash
# Submit the FACE-HD method/reference alignment audit to Slurm.

set -euo pipefail

DRY_RUN=false
OUTPUT_ROOT="diagnosis/method_alignment/audit_results"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,40p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --output-root) OUTPUT_ROOT="${2:?Missing output root}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p "${OUTPUT_ROOT}"

SBATCH_CMD=(
    sbatch
    --partition="${PARTITION}"
    --output="${OUTPUT_ROOT}/%x_%j.out"
    --error="${OUTPUT_ROOT}/%x_%j.err"
    --export="ALL,METHOD_AUDIT_OUTPUT_ROOT=${OUTPUT_ROOT}"
    diagnosis/method_alignment/method_alignment_audit.cmd
)

echo "======================================================"
echo " FACE-HD Method Alignment Audit Submission"
echo "======================================================"
echo " Repo:        ${REPO_DIR}"
echo " Output root: ${OUTPUT_ROOT}"
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
