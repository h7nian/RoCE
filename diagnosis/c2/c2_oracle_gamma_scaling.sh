#!/bin/bash
# ============================================================================
# RoCE C2 oracle-gamma coverage-vs-n submission script
# ----------------------------------------------------------------------------
# Confirmatory: does the ORACLE (true density-ratio gamma) estimator keep
# ~0.95 coverage at n=20000 while the FITTED-gamma estimator degrades to ~0.84?
# If yes, the estimand/identification is sound and the C2 under-coverage is a
# gamma-hat ESTIMATION problem (density-ratio CV validation-loss scale in
# src/cv_utils.h). Reuses unmodified diagnosis/c2/c2_oracle_gamma_probe.R.
#
# Cells (K,n): {(3,5000),(3,20000),(4,5000),(4,20000)}; SEEDS_PER_CELL each.
# Default 60 seeds/cell x 4 cells = 240 array tasks.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
TIME_LIMIT="02:00:00"
MEMORY="24g"
CPUS="2"
SEEDS_PER_CELL="60"
SEED_BASE="2000"
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
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

export C2_ORACLE_GAMMA_OUTPUT_ROOT="diagnosis/c2/oracle_gamma_scaling"
mkdir -p "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/log" "${C2_ORACLE_GAMMA_OUTPUT_ROOT}/results"

N_CELLS=4
if [[ -z "${TASKS}" ]]; then
    TASKS="1-$(( SEEDS_PER_CELL * N_CELLS ))"
fi

EXPORT_VARS="C2_ORACLE_GAMMA_OUTPUT_ROOT=${C2_ORACLE_GAMMA_OUTPUT_ROOT},C2_ORC_SEEDS_PER_CELL=${SEEDS_PER_CELL},C2_ORC_SEED_BASE=${SEED_BASE}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_oracle_gamma_scaling.cmd
)

echo "======================================================"
echo " RoCE C2 Oracle-Gamma Coverage-vs-n Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " Seeds/cell:     ${SEEDS_PER_CELL}  (base ${SEED_BASE})"
echo " Cells:          (3,5000) (3,20000) (4,5000) (4,20000)"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Partition:      ${PARTITION}"
echo " Output root:    ${C2_ORACLE_GAMMA_OUTPUT_ROOT}"
echo " Dry run:        ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
