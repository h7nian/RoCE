#!/bin/bash
# Submit local method-reference PDF extraction to Slurm.

set -euo pipefail

DRY_RUN=false
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,32p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p diagnosis/method_alignment/audit_results

SBATCH_CMD=(
    sbatch
    --partition="${PARTITION}"
    diagnosis/method_alignment/extract_method_pdfs.cmd
)

echo "======================================================"
echo " RoCE PDF Extraction Submission"
echo "======================================================"
echo " Repo:      ${REPO_DIR}"
echo " Partition: ${PARTITION}"
echo " Dry run:   ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
