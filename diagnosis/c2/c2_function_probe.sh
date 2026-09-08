#!/bin/bash
# ============================================================================
# RoCE C2 function-level probe submission script
# ----------------------------------------------------------------------------
# Submits one-seed, fold-level C2 probes to MSI.  These jobs are diagnostic
# only: they record target nuisance decompositions, source transport/oracle
# variants, and aggregation weights under the default K_f=10 setting.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="04:00:00"
MEMORY="24g"
CPUS="4"
N_CORES="1"
N_FOLDS="10"
NLAMBDA_INIT="100"
MODES="one_round,two_round"
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
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p diagnosis/c2/function_probe/log diagnosis/c2/function_probe/results

if [[ -z "${TASKS}" ]]; then
    TASKS="1-10"
fi

EXPORT_MODES="${MODES//,/:}"
EXPORT_VARS="C2_FN_FOLDS=${N_FOLDS},C2_FN_NLAMBDA_INIT=${NLAMBDA_INIT},C2_FN_CORES=${N_CORES},C2_FN_MODES=${EXPORT_MODES}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_function_probe.cmd
)

echo "======================================================"
echo " RoCE C2 Function Probe Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " K_f:            ${N_FOLDS}"
echo " nlambda:        ${NLAMBDA_INIT}"
echo " n_cores/R:      ${N_CORES}"
echo " Modes:          ${MODES}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Cpus/task:      ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Logs:           diagnosis/c2/function_probe/log/"
echo " Results:        diagnosis/c2/function_probe/results/"
echo " Dry run:        ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
