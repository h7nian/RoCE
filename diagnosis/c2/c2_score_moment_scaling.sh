#!/bin/bash
# ============================================================================
# FACE-HD C2 score/moment n-SCALING submission script
# ----------------------------------------------------------------------------
# Reuses the (proven, unmodified) diagnosis/c2/c2_score_moment_audit.R on the
# failing C2 p=10 settings across n_total in {5000, 10000, 20000}. Diagnosis
# only; estimator source is NOT modified. Purpose: decide whether the source
# correction underfill (out_noise = E_s[w (Y - m_true)] ~ -required, observed
# delta ~ 0, bias ~ -required) is a finite-sample higher-order effect that
# vanishes as n grows, or a structural method issue that persists.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="12:00:00"
MEMORY="48g"
CPUS="4"
N_CORES="1"
N_FOLDS="10"
NLAMBDA_INIT="100"
MODES="one_round,two_round"
USE_LAMBDA_CACHE="true"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,16p' "$0"
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
        --no-cache) USE_LAMBDA_CACHE="false"; shift ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

export C2_SCORE_AUDIT_OUTPUT_ROOT="diagnosis/c2/score_moment_scaling"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"

if [[ -z "${TASKS}" ]]; then
    TASKS="1-12"
fi

EXPORT_MODES="${MODES//,/:}"
EXPORT_VARS="C2_SCORE_AUDIT_FOLDS=${N_FOLDS},C2_SCORE_AUDIT_NLAMBDA_INIT=${NLAMBDA_INIT},C2_SCORE_AUDIT_CORES=${N_CORES},C2_SCORE_AUDIT_MODES=${EXPORT_MODES},C2_SCORE_AUDIT_USE_LAMBDA_CACHE=${USE_LAMBDA_CACHE},C2_SCORE_AUDIT_OUTPUT_ROOT=${C2_SCORE_AUDIT_OUTPUT_ROOT}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_score_moment_scaling.cmd
)

echo "======================================================"
echo " FACE-HD C2 Score/Moment n-Scaling Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " K_f:            ${N_FOLDS}"
echo " nlambda:        ${NLAMBDA_INIT}"
echo " n_cores/R:      ${N_CORES}"
echo " Modes:          ${MODES}"
echo " Lambda cache:   ${USE_LAMBDA_CACHE}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Cpus/task:      ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Output root:    ${C2_SCORE_AUDIT_OUTPUT_ROOT}"
echo " Dry run:        ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
