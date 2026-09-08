#!/bin/bash
# ============================================================================
# RoCE bias-diagnosis fast probe submission
# ----------------------------------------------------------------------------
# Submits a SHORT single-task SLURM job that runs diagnosis/bias/probe.R
# (3 phases: A data sanity, B per-function probe, C end-to-end seed sweep).
#
# Designed for FAST iteration — not an array, not 12 h walltime. Default
# 30 min on a small partition so it lands quickly.
#
# Usage:
#   bash diagnosis/bias/probe.sh                          # phase ALL, anchor setting
#   bash diagnosis/bias/probe.sh --phase B                # only per-function probe
#   bash diagnosis/bias/probe.sh --phase C --seeds 20     # 20 seeds end-to-end
#   bash diagnosis/bias/probe.sh --n 2000 --p 50          # change setting
#   bash diagnosis/bias/probe.sh --dry-run
# ============================================================================

set -euo pipefail

DRY_RUN=false
PHASE="ALL"
N_TOTAL=1000
K=3
P=20
N_SEEDS=10
TIME_LIMIT="00:30:00"
MEMORY="16g"
CPUS="8"
PARTITION="msismall,amdsmall,agsmall"

usage() {
    sed -n '1,30p' "$0"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY_RUN=true; shift ;;
        --phase)     PHASE="${2:?Missing phase}"; shift 2 ;;
        --n)         N_TOTAL="${2:?Missing n}"; shift 2 ;;
        --K)         K="${2:?Missing K}"; shift 2 ;;
        --p)         P="${2:?Missing p}"; shift 2 ;;
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

EXPORT_VARS="PROBE_PHASE=${PHASE},PROBE_N=${N_TOTAL},PROBE_K=${K},PROBE_P=${P},PROBE_SEEDS=${N_SEEDS}"

SBATCH_CMD=(
    sbatch
    --job-name="FACE-probe-${PHASE}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/bias/probe.cmd
)

echo "======================================================"
echo " RoCE Bias Probe Submission"
echo "======================================================"
echo " Phase:       ${PHASE}"
echo " Setting:     n=${N_TOTAL} K=${K} p=${P} seeds=${N_SEEDS}"
echo " Walltime:    ${TIME_LIMIT}"
echo " Mem:         ${MEMORY}"
echo " Cpus:        ${CPUS}"
echo " Partition:   ${PARTITION}"
echo " Logs:        diagnosis/bias/log/FACE-probe-${PHASE}_*.out"
echo " Dump:        diagnosis/bias/probe_out/"
echo " Dry run:     ${DRY_RUN}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
