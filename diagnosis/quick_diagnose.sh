#!/bin/bash
# ============================================================================
# RoCE Quick Diagnosis Submission Script
# ============================================================================
#
# Submits a small MSI Slurm job that runs static and package-load diagnostics.
# It does not run Monte Carlo simulations or submit experiment arrays.
#
# Usage:
#   bash diagnosis/quick_diagnose.sh
#   bash diagnosis/quick_diagnose.sh --dry-run
#   bash diagnosis/quick_diagnose.sh --install-check
#   bash diagnosis/quick_diagnose.sh --reset-checkpoint
#   bash diagnosis/quick_diagnose.sh --no-requeue
#
# Options:
#   --dry-run            Print the sbatch command without submitting.
#   --install-check      Also run a temporary R CMD INSTALL in the Slurm job.
#   --reset-checkpoint   Wipe any existing checkpoint for the new job before
#                        running the first attempt. Has no effect on later
#                        requeues — only consumed once per submission.
#   --no-requeue         Disable automatic SLURM requeue on preemption (the
#                        .cmd file defaults to #SBATCH --requeue).
#   --time HH:MM:SS      Override Slurm walltime, default 00:30:00.
#   --mem MEM            Override Slurm memory, default 8g.
#   --partition PARTS    Override partition list, default
#                        "msismall,preempt,amdsmall,agsmall".
#
# Checkpoint / Requeue Behavior:
#   The diagnostic job is registered as --requeue with a 90s SIGUSR1 grace
#   period (see diagnosis/quick_diagnose.cmd). On preemption:
#     1. SLURM sends SIGUSR1, our handler records state and runs
#        `scontrol requeue ${SLURM_JOB_ID}`.
#     2. The job is relaunched under the SAME job id.
#     3. Each diagnostic step is keyed on its name in
#        diagnosis/tmp/checkpoint_diag_<job>.txt; previously-completed steps
#        are skipped, so the new attempt resumes where the prior one stopped.
#     4. On successful completion the checkpoint files are deleted.
#
# ============================================================================

set -euo pipefail

DRY_RUN=false
INSTALL_CHECK=false
RESET_CHECKPOINT=false
NO_REQUEUE=false
TIME_LIMIT="00:30:00"
MEMORY="8g"
PARTITION="msismall,preempt,amdsmall,agsmall"

usage() {
    sed -n '1,42p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --install-check) INSTALL_CHECK=true; shift ;;
        --reset-checkpoint) RESET_CHECKPOINT=true; shift ;;
        --no-requeue) NO_REQUEUE=true; shift ;;
        --time) TIME_LIMIT="${2:?Missing time value}"; shift 2 ;;
        --mem) MEMORY="${2:?Missing mem value}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition value}"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_DIR}"

mkdir -p diagnosis/log diagnosis/results diagnosis/tmp

INSTALL_FLAG="FALSE"
if [[ "${INSTALL_CHECK}" == true ]]; then
    INSTALL_FLAG="TRUE"
fi
RESET_FLAG="FALSE"
if [[ "${RESET_CHECKPOINT}" == true ]]; then
    RESET_FLAG="TRUE"
fi

SBATCH_CMD=(
    sbatch
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --partition="${PARTITION}"
    --export="ALL,DIAG_INSTALL_CHECK=${INSTALL_FLAG},DIAG_RESET_CHECKPOINT=${RESET_FLAG}"
)
if [[ "${NO_REQUEUE}" == true ]]; then
    SBATCH_CMD+=(--no-requeue)
fi
SBATCH_CMD+=(diagnosis/quick_diagnose.cmd)

REQUEUE_STATUS="enabled (#SBATCH --requeue in .cmd)"
if [[ "${NO_REQUEUE}" == true ]]; then
    REQUEUE_STATUS="disabled (--no-requeue passed to sbatch)"
fi

echo "======================================================"
echo " RoCE Quick Diagnosis Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Dry run:        ${DRY_RUN}"
echo " Install check:  ${INSTALL_FLAG}"
echo " Reset chkpt:    ${RESET_FLAG}"
echo " Requeue:        ${REQUEUE_STATUS}"
echo " Time:           ${TIME_LIMIT}"
echo " Memory:         ${MEMORY}"
echo " Partition:      ${PARTITION}"
echo " Logs:           diagnosis/log/"
echo " Reports:        diagnosis/results/"
echo " Checkpoints:    diagnosis/tmp/checkpoint_diag_<job>.txt"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
