#!/bin/bash
# Submit the Slurm-side summarizer for C2 nuisance-refit probe outputs.

set -euo pipefail

DRY_RUN=false
DEPENDENCY=""
TIME_LIMIT="00:20:00"
MEMORY="4g"
CPUS="1"
PARTITION="amdsmall,agsmall,msismall,msilarge"

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
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p diagnosis/c2/nuisance_refit_probe/log diagnosis/c2/nuisance_refit_probe/summary

SBATCH_CMD=(
    sbatch
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
)

if [[ -n "${DEPENDENCY}" ]]; then
    SBATCH_CMD+=(--dependency="${DEPENDENCY}")
fi

SBATCH_CMD+=(diagnosis/c2/c2_nuisance_refit_summarize.cmd)

echo "======================================================"
echo " RoCE C2 Nuisance Refit Summary Submission"
echo "======================================================"
echo " Repo:       ${REPO_DIR}"
echo " Dependency: ${DEPENDENCY:-none}"
echo " Time:       ${TIME_LIMIT}"
echo " Mem:        ${MEMORY}"
echo " Cpus:       ${CPUS}"
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
