#!/bin/bash
# Submit fast p=10 C2 oracle-gamma convention diagnosis jobs.

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="00:30:00"
MEMORY="8g"
CPUS="1"
OUTPUT_ROOT="diagnosis/c2/oracle_gamma_probe"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,34p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --tasks) TASKS="${2:?Missing tasks list}"; shift 2 ;;
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

mkdir -p "${OUTPUT_ROOT}/log" "${OUTPUT_ROOT}/results"

if [[ -z "${TASKS}" ]]; then
    TASKS="1-8"
fi

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --output="${OUTPUT_ROOT}/log/%x_%A_%a.out"
    --error="${OUTPUT_ROOT}/log/%x_%A_%a.err"
    --export="ALL,C2_ORACLE_GAMMA_OUTPUT_ROOT=${OUTPUT_ROOT}"
    diagnosis/c2/c2_oracle_gamma_probe.cmd
)

echo "======================================================"
echo " RoCE C2 Oracle Gamma Probe Submission"
echo "======================================================"
echo " Repo:        ${REPO_DIR}"
echo " Tasks:       ${TASKS}"
echo " Output root: ${OUTPUT_ROOT}"
echo " Time/task:   ${TIME_LIMIT}"
echo " Mem/task:    ${MEMORY}"
echo " Cpus/task:   ${CPUS}"
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
