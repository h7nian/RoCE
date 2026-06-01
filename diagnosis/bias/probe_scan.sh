#!/bin/bash
# ============================================================================
# FACE-HD bias probe scan submission
# ----------------------------------------------------------------------------
# Submits a single-task SLURM job that runs diagnosis/bias/probe_scan.R
# over 6 (n, p) settings × 3 seeds to quantify the relationship between
# (p / n_calib_arm) and mu_pred_ts instability.
#
# Usage:
#   bash diagnosis/bias/probe_scan.sh
#   bash diagnosis/bias/probe_scan.sh --seeds 5
#   bash diagnosis/bias/probe_scan.sh --time 02:00:00
#   bash diagnosis/bias/probe_scan.sh --dry-run
# ============================================================================

set -euo pipefail

DRY_RUN=false
N_SEEDS=3
TIME_LIMIT="01:00:00"
MEMORY="16g"
CPUS="4"
PARTITION="msismall,amdsmall,agsmall"

usage() {
    sed -n '1,15p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY_RUN=true; shift ;;
        --seeds)     N_SEEDS="${2:?Missing seeds}"; shift 2 ;;
        --time)      TIME_LIMIT="${2:?Missing time}"; shift 2 ;;
        --mem)       MEMORY="${2:?Missing mem}"; shift 2 ;;
        --cpus)      CPUS="${2:?Missing cpus}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition}"; shift 2 ;;
        -h|--help)   usage ;;
        *) echo "Error: unknown option '$1'"; usage ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
cd "${REPO_DIR}"
mkdir -p diagnosis/bias/log diagnosis/bias/probe_out

EXPORT_VARS="PROBE_SCAN_SEEDS=${N_SEEDS}"

SBATCH_CMD=(
    sbatch
    --job-name="FACE-probe-scan"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/bias/probe_scan.cmd
)

echo "======================================================"
echo " FACE-HD Bias Probe Scan Submission"
echo "======================================================"
echo " Settings:    6 (n, p) combos at K=3, C1, ss=0.5"
echo "                (1000,10) (1000,20) (1000,50)"
echo "                (2000,20) (5000,10) (5000,20)"
echo " Seeds/cell:  ${N_SEEDS}"
echo " Walltime:    ${TIME_LIMIT}"
echo " Mem:         ${MEMORY}"
echo " Cpus:        ${CPUS}"
echo " Partition:   ${PARTITION}"
echo " Logs:        diagnosis/bias/log/FACE-probe-scan_*.out"
echo " Dump:        diagnosis/bias/probe_out/probe_scan_*.csv"
echo " Dry run:     ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
