#!/bin/bash
# ============================================================================
# RoCE Bias Diagnosis Submission Script
# ----------------------------------------------------------------------------
# Submits a SLURM array job that runs 11 small simulation settings, all in
# config C1 with binary outcome + mild transform + superpopulation estimand.
# The grid varies n / p / K / shift_strength one-axis-at-a-time around an
# anchor point (n=1000, K=3, p=20, shift=0.5).
#
# Goal: identify which axis drives the positive bias of one_round_crossfit /
# two_round_crossfit in C1. THE ALGORITHM IS NOT MODIFIED.
#
# Grid (idx -> task):
#   1  anchor          n=1000 K=3 p=20  ss=0.5
#   2  n_small         n= 500 K=3 p=20  ss=0.5    (bias vs n)
#   3  n_large         n=2000 K=3 p=20  ss=0.5    (bias vs n)
#   4  p_50            n=1000 K=3 p=50  ss=0.5    (bias vs p, high-dim)
#   5  p_100           n=1000 K=3 p=100 ss=0.5    (bias vs p, high-dim)
#   6  p_200           n=1000 K=3 p=200 ss=0.5    (bias vs p, high-dim)
#   7  K_4             n=1000 K=4 p=20  ss=0.5    (bias vs K)
#   8  K_5             n=1000 K=5 p=20  ss=0.5    (bias vs K)
#   9  shift_005       n=1000 K=3 p=20  ss=0.05   (near-zero shift, isolate outcome model)
#  10  shift_1         n=1000 K=3 p=20  ss=1.0    (stress shift)
#
# Revisions from the first run (job 9276320):
#   * Dropped K_2 (algorithm degenerate at K=2; see earlier failed run).
#   * Replaced shift_0 with shift_005 (validator rejects shift=0).
#   * Switched to non-preempt partitions to avoid the requeue death-spiral
#     observed on saffo-2tb where most attempts were killed before
#     CHECKPOINT_SIM_INTERVAL (20) sims completed, stranding all progress.
#   * Cut n_sims default 100 -> 50 so a single 12 h walltime suffices for
#     most tasks (still gives bias_mean MC SE ~0.005/sqrt(50) ~ 7e-4).
#   * run_setting.R overrides ckpt_config$sim_interval to 5 so a USR1 hit
#     between sims 5..50 saves real progress.
#
# Usage:
#   bash diagnosis/bias/bias_diagnose.sh                # submit all
#   bash diagnosis/bias/bias_diagnose.sh --dry-run      # show sbatch only
#   bash diagnosis/bias/bias_diagnose.sh --tasks 1,4,5  # subset
#   bash diagnosis/bias/bias_diagnose.sh --n-sims 50    # override n_sims
#   bash diagnosis/bias/bias_diagnose.sh --reset-checkpoint  # wipe stored checkpoints first
#
# Checkpoint / requeue:
#   The .cmd script sets `#SBATCH --requeue` and `#SBATCH --signal=B:USR1@120`
#   so SLURM can preempt and re-launch the SAME array task with the SAME
#   SLURM_JOB_ID. The R loop (run_simulation_study) saves a per-setting
#   checkpoint under diagnosis/bias/checkpoints/task_<id>/, and resumes from
#   the saved simulation index on the next attempt. Use --reset-checkpoint
#   to start a fresh sweep (deletes diagnosis/bias/checkpoints/).
# ============================================================================

set -euo pipefail

DRY_RUN=false
TASKS=""
N_SIMS=50
TIME_LIMIT="12:00:00"
MEMORY="32g"
CPUS="16"
# Non-preempt partitions only. preempt/saffo-2tb were killing requeues
# faster than CHECKPOINT_SIM_INTERVAL could save progress; switching to
# stable partitions trades queue wait for completion certainty.
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

usage() {
    sed -n '1,40p' "$0"
    exit 1
}

RESET_CHECKPOINT=false

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

mkdir -p diagnosis/bias/log diagnosis/bias/results diagnosis/bias/checkpoints diagnosis/bias/tmp

if [[ "${RESET_CHECKPOINT}" == true ]]; then
    echo "[bias_diagnose] Wiping diagnosis/bias/checkpoints/ before submission..."
    rm -rf diagnosis/bias/checkpoints
    mkdir -p diagnosis/bias/checkpoints
fi

# Default task range: 1..10 (the full grid declared in bias_diagnose.cmd).
if [[ -z "${TASKS}" ]]; then
    TASKS="1-10"
fi

EXPORT_VARS="DIAG_N_SIMS=${N_SIMS}"

SBATCH_CMD=(
    sbatch
    --array="${TASKS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/bias/bias_diagnose.cmd
)

echo "======================================================"
echo " RoCE Bias Diagnosis Submission"
echo "======================================================"
echo " Repo:           ${REPO_DIR}"
echo " Tasks:          ${TASKS}"
echo " n_sims/task:    ${N_SIMS}"
echo " Time/task:      ${TIME_LIMIT}"
echo " Mem/task:       ${MEMORY}"
echo " Cpus/task:      ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Requeue:        enabled (#SBATCH --requeue in .cmd)"
echo " Preempt signal: SIGUSR1 @ 120 s grace"
echo " Logs:           diagnosis/bias/log/"
echo " Results:        diagnosis/bias/results/"
echo " Checkpoints:    diagnosis/bias/checkpoints/task_<id>/"
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
