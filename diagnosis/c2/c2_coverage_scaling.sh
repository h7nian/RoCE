#!/bin/bash
# ============================================================================
# FACE-HD C2 coverage-scaling submission script
# ----------------------------------------------------------------------------
# Estimates actual CI coverage for the C2 p=10 settings across many seeds at
# small (n=5000) vs large (n=20000) sample size, to confirm whether the
# coverage dip is finite-sample (recovers with n). Reuses the proven,
# UNMODIFIED diagnosis/c2/c2_score_moment_audit.R. Diagnosis only.
#
# Cells (K, n): {(3,5000),(3,20000),(4,5000),(4,20000)}; SEEDS_PER_CELL each.
# Default 50 seeds/cell x 4 cells = 200 array tasks (--array 1-200).
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="04:00:00"
MEMORY="32g"
CPUS="4"
N_CORES="1"
N_FOLDS="10"
NLAMBDA_INIT="100"
MODES="one_round"
USE_LAMBDA_CACHE="true"
SEEDS_PER_CELL="50"
SEED_BASE="1000"
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
        --seeds)   SEEDS_PER_CELL="${2:?Missing seeds value}"; shift 2 ;;
        --seed-base) SEED_BASE="${2:?Missing seed-base value}"; shift 2 ;;
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

export C2_SCORE_AUDIT_OUTPUT_ROOT="diagnosis/c2/coverage_scaling"
mkdir -p "${C2_SCORE_AUDIT_OUTPUT_ROOT}/log" "${C2_SCORE_AUDIT_OUTPUT_ROOT}/results"

N_CELLS=4
if [[ -z "${TASKS}" ]]; then
    TASKS="1-$(( SEEDS_PER_CELL * N_CELLS ))"
fi

EXPORT_MODES="${MODES//,/:}"
EXPORT_VARS="C2_SCORE_AUDIT_FOLDS=${N_FOLDS},C2_SCORE_AUDIT_NLAMBDA_INIT=${NLAMBDA_INIT},C2_SCORE_AUDIT_CORES=${N_CORES},C2_SCORE_AUDIT_MODES=${EXPORT_MODES},C2_SCORE_AUDIT_USE_LAMBDA_CACHE=${USE_LAMBDA_CACHE},C2_SCORE_AUDIT_OUTPUT_ROOT=${C2_SCORE_AUDIT_OUTPUT_ROOT},C2_COV_SEEDS_PER_CELL=${SEEDS_PER_CELL},C2_COV_SEED_BASE=${SEED_BASE}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_coverage_scaling.cmd
)

echo "======================================================"
echo " FACE-HD C2 Coverage-Scaling Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " Seeds/cell:     ${SEEDS_PER_CELL}  (base ${SEED_BASE})"
echo " Cells:          (3,5000) (3,20000) (4,5000) (4,20000)"
echo " K_f:            ${N_FOLDS}"
echo " Modes:          ${MODES}"
echo " Lambda cache:   ${USE_LAMBDA_CACHE}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
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
