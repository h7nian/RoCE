#!/bin/bash
# ============================================================================
# Submits a short SLURM job that runs diagnosis/bias/probe_plot.R to
# render PNG plots from the probe_scan / probe_kf CSV outputs.
# ============================================================================

set -euo pipefail

DRY_RUN=false
TIME_LIMIT="00:15:00"
MEMORY="4g"
PARTITION="msismall,amdsmall,agsmall"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY_RUN=true; shift ;;
        --time)      TIME_LIMIT="${2:?}"; shift 2 ;;
        --mem)       MEMORY="${2:?}"; shift 2 ;;
        --partition) PARTITION="${2:?}"; shift 2 ;;
        -h|--help)   sed -n '1,8p' "$0"; exit 1 ;;
        *) echo "Error: unknown option '$1'"; exit 1 ;;
    esac
done

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p diagnosis/bias/log diagnosis/bias/probe_out/plots

SBATCH_CMD=(
    sbatch
    --job-name="FACE-probe-plot"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task=2
    --partition="${PARTITION}"
    --export="ALL"
    diagnosis/bias/probe_plot.cmd
)

echo "[probe_plot] partition=${PARTITION}  time=${TIME_LIMIT}  mem=${MEMORY}"
if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
