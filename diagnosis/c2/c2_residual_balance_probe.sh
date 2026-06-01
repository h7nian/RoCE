#!/bin/bash
# ============================================================================
# FACE-HD C2 residual-balance probe submission script
# ----------------------------------------------------------------------------
# Submits diagnosis-only jobs that decompose source correction into required,
# observed, true-outcome-residual, and true-propensity residual balance terms.
# Estimator source is not modified.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="06:00:00"
MEMORY="24g"
CPUS="1"
N_CORES="1"
N_FOLDS="10"
NLAMBDA_INIT="100"
MODES="one_round,two_round"
OUTPUT_ROOT="diagnosis/c2/residual_balance_probe"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,34p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --tasks)   TASKS="${2:?Missing tasks list}"; shift 2 ;;
        --time)    TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem)     MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --cpus)    CPUS="${2:?Missing cpus value}"; shift 2 ;;
        --ncores)  N_CORES="${2:?Missing ncores value}"; shift 2 ;;
        --folds)   N_FOLDS="${2:?Missing folds value}"; shift 2 ;;
        --nlambda) NLAMBDA_INIT="${2:?Missing nlambda value}"; shift 2 ;;
        --modes)   MODES="${2:?Missing modes value}"; shift 2 ;;
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
    TASKS="1-5"
fi

EXPORT_MODES="${MODES//,/:}"
EXPORT_VARS="C2_RBAL_FOLDS=${N_FOLDS},C2_RBAL_NLAMBDA_INIT=${NLAMBDA_INIT},C2_RBAL_CORES=${N_CORES},C2_RBAL_MODES=${EXPORT_MODES},C2_RBAL_OUTPUT_ROOT=${OUTPUT_ROOT}"

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
    diagnosis/c2/c2_residual_balance_probe.cmd
)

echo "======================================================"
echo " FACE-HD C2 Residual Balance Probe Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " K_f:            ${N_FOLDS}"
echo " nlambda:        ${NLAMBDA_INIT}"
echo " n_cores/R:      ${N_CORES}"
echo " Modes:          ${MODES}"
echo " Output root:    ${OUTPUT_ROOT}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Cpus/task:      ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Logs:           ${OUTPUT_ROOT}/log/"
echo " Results:        ${OUTPUT_ROOT}/results/"
echo " Dry run:        ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
