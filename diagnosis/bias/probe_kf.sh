#!/bin/bash
# ============================================================================
# RoCE K_f probe scan submission
# ----------------------------------------------------------------------------
# Submits a single-task SLURM job that runs diagnosis/bias/probe_kf.R,
# sweeping K_f ∈ {3, 5, 7, 10} at two settings (catastrophic / production).
# ============================================================================

set -euo pipefail

DRY_RUN=false
N_SEEDS=3
TIME_LIMIT="01:30:00"
MEMORY="16g"
CPUS="4"
PARTITION="msismall,amdsmall,agsmall"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY_RUN=true; shift ;;
        --seeds)     N_SEEDS="${2:?Missing seeds}"; shift 2 ;;
        --time)      TIME_LIMIT="${2:?Missing time}"; shift 2 ;;
        --mem)       MEMORY="${2:?Missing mem}"; shift 2 ;;
        --cpus)      CPUS="${2:?Missing cpus}"; shift 2 ;;
        --partition) PARTITION="${2:?Missing partition}"; shift 2 ;;
        -h|--help)   sed -n '1,12p' "$0"; exit 1 ;;
        *) echo "Error: unknown option '$1'"; exit 1 ;;
    esac
done

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p diagnosis/bias/log diagnosis/bias/probe_out

EXPORT_VARS="PROBE_KF_SEEDS=${N_SEEDS}"

SBATCH_CMD=(
    sbatch
    --job-name="FACE-probe-kf"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/bias/probe_kf.cmd
)

echo "======================================================"
echo " RoCE K_f Probe Scan"
echo "======================================================"
echo " K_f values:  3, 5, 7, 10"
echo " Settings:    catastrophic (n=1000,p=20), production (n=5000,p=10)"
echo " Seeds/cell:  ${N_SEEDS}"
echo " Walltime:    ${TIME_LIMIT}"
echo " Partition:   ${PARTITION}"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
