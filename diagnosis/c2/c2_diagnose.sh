#!/bin/bash
# ============================================================================
# FACE-HD C2 diagnosis submission script
# ----------------------------------------------------------------------------
# Submits a SLURM array for measurement-only C2 diagnostics.  The grid is
# intentionally small and hypothesis-driven:
#   * tasks 1-12: exact problematic C2 settings by estimand at default K_f=10
#   * tasks 13-14: larger-n C2 checks at default K_f=10
#
# This mirrors diagnosis/bias/*.sh: all experiments are submitted to MSI through
# sbatch, write logs under diagnosis/c2/log/, and persist CSV/RData outputs
# under diagnosis/c2/results/.  No estimator code is modified here.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
N_SIMS=100
TIME_LIMIT="06:00:00"
MEMORY="32g"
CPUS="16"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"
RESET_CHECKPOINT=false

usage() {
    sed -n '1,38p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --tasks)   TASKS="${2:?Missing tasks list}"; shift 2 ;;
        --n-sims)  N_SIMS="${2:?Missing n-sims value}"; shift 2 ;;
        --time)    TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem)     MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --cpus)    CPUS="${2:?Missing cpus value}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        --reset-checkpoint) RESET_CHECKPOINT=true; shift ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"

mkdir -p diagnosis/c2/log diagnosis/c2/results diagnosis/c2/checkpoints diagnosis/c2/tmp

if [[ "${RESET_CHECKPOINT}" == true ]]; then
    echo "[c2_diagnose] Wiping diagnosis/c2/checkpoints/ before submission..."
    rm -rf diagnosis/c2/checkpoints
    mkdir -p diagnosis/c2/checkpoints
fi

if [[ -z "${TASKS}" ]]; then
    TASKS="1-14"
fi

EXPORT_VARS="C2_DIAG_N_SIMS=${N_SIMS}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/c2/c2_diagnose.cmd
)

echo "======================================================"
echo " FACE-HD C2 Diagnosis Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " n_sims/task:    ${N_SIMS}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Cpus/task:      ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Logs:           diagnosis/c2/log/"
echo " Results:        diagnosis/c2/results/"
echo " Checkpoints:    diagnosis/c2/checkpoints/task_<id>/"
echo " Reset chkpt:    ${RESET_CHECKPOINT}"
echo " Dry run:        ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
