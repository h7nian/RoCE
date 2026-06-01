#!/bin/bash
# Submit C2 score/moment audit jobs to Slurm.

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="04:00:00"
MEMORY="24g"
CPUS="1"
N_CORES="1"
N_FOLDS="10"
NLAMBDA_INIT="100"
MODES="one_round,two_round"
USE_LAMBDA_CACHE="true"
OUTPUT_ROOT="diagnosis/c2/score_moment_audit"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,42p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --tasks) TASKS="${2:?Missing tasks list}"; shift 2 ;;
        --time) TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem) MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --cpus) CPUS="${2:?Missing cpus value}"; shift 2 ;;
        --ncores) N_CORES="${2:?Missing ncores value}"; shift 2 ;;
        --folds) N_FOLDS="${2:?Missing folds value}"; shift 2 ;;
        --nlambda) NLAMBDA_INIT="${2:?Missing nlambda value}"; shift 2 ;;
        --modes) MODES="${2:?Missing modes value}"; shift 2 ;;
        --use-lambda-cache) USE_LAMBDA_CACHE="${2:?Missing true/false value}"; shift 2 ;;
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

EXPORT_MODES="${MODES//,/:}"
EXPORT_VARS="C2_SCORE_AUDIT_FOLDS=${N_FOLDS},C2_SCORE_AUDIT_NLAMBDA_INIT=${NLAMBDA_INIT},C2_SCORE_AUDIT_CORES=${N_CORES},C2_SCORE_AUDIT_MODES=${EXPORT_MODES},C2_SCORE_AUDIT_USE_LAMBDA_CACHE=${USE_LAMBDA_CACHE},C2_SCORE_AUDIT_OUTPUT_ROOT=${OUTPUT_ROOT}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --output="${OUTPUT_ROOT}/log/%x_%A_%a.out"
    --error="${OUTPUT_ROOT}/log/%x_%A_%a.err"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_score_moment_audit.cmd
)

echo "======================================================"
echo " FACE-HD C2 Score/Moment Audit Submission"
echo "======================================================"
echo " Repo:        ${REPO_DIR}"
echo " Tasks:       ${TASKS}"
echo " K_f:         ${N_FOLDS}"
echo " nlambda:     ${NLAMBDA_INIT}"
echo " n_cores/R:   ${N_CORES}"
echo " Modes:       ${MODES}"
echo " Lambda cache:${USE_LAMBDA_CACHE}"
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
