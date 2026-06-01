#!/bin/bash
# ============================================================================
# FACE-HD K_f=3 production-scale verification submission
# ----------------------------------------------------------------------------
# Submits a single-task SLURM job that re-runs the exact production setting
#   n=5000, K=3, p=10, C1, superpopulation, mild, ss=0.5
# with n_folds = 3 (default) instead of 10, n_sims = 200.
#
# Expected: K_f=3 → one_round / two_round coverage moves from
#   0.886 / 0.878 (production at K_f=10) to ~0.94+.
#
# Usage:
#   bash diagnosis/bias/verify_kf3.sh
#   bash diagnosis/bias/verify_kf3.sh --kf 5         # also try K_f=5
#   bash diagnosis/bias/verify_kf3.sh --nsims 500
#   bash diagnosis/bias/verify_kf3.sh --dry-run
# ============================================================================

set -euo pipefail

DRY_RUN=false
N_FOLDS=3
N_SIMS=200
CONFIG="C1"
TIME_LIMIT="02:00:00"
MEMORY="32g"
CPUS="16"
PARTITION="msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)   DRY_RUN=true; shift ;;
        --kf)        N_FOLDS="${2:?Missing kf}"; shift 2 ;;
        --nsims)     N_SIMS="${2:?Missing nsims}"; shift 2 ;;
        --config)    CONFIG="${2:?Missing config}"; shift 2 ;;
        --time)      TIME_LIMIT="${2:?}"; shift 2 ;;
        --mem)       MEMORY="${2:?}"; shift 2 ;;
        --cpus)      CPUS="${2:?}"; shift 2 ;;
        --partition) PARTITION="${2:?}"; shift 2 ;;
        -h|--help)   sed -n '1,20p' "$0"; exit 1 ;;
        *) echo "Error: unknown option '$1'"; exit 1 ;;
    esac
done

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
mkdir -p diagnosis/bias/log diagnosis/bias/verify_out diagnosis/bias/checkpoints

EXPORT_VARS="VERIFY_N_FOLDS=${N_FOLDS},VERIFY_N_SIMS=${N_SIMS},VERIFY_CONFIG=${CONFIG}"

SBATCH_CMD=(
    sbatch
    --job-name="FACE-verify-${CONFIG}-kf${N_FOLDS}"
    --time="${TIME_LIMIT}"
    --mem="${MEMORY}"
    --cpus-per-task="${CPUS}"
    --partition="${PARTITION}"
    --export="ALL,${EXPORT_VARS}"
    diagnosis/bias/verify_kf3.cmd
)

echo "======================================================"
echo " FACE-HD ${CONFIG} K_f=${N_FOLDS} Production Verification"
echo "======================================================"
echo " Pinned setting: n=5000 K=3 p=10 ${CONFIG} ss=0.5"
echo " K_f:            ${N_FOLDS}      (production was 10)"
echo " n_sims:         ${N_SIMS}"
echo " Walltime:       ${TIME_LIMIT}"
echo " Cpus:           ${CPUS}"
echo " Partition:      ${PARTITION}"
echo " Output:         diagnosis/bias/verify_out/verify_${CONFIG}_n5000_K3_p10_kf${N_FOLDS}_summary.csv"
echo " Compare with:   results/n5000_K3_p10_${CONFIG}_..._kf10_..._summary.csv"
echo "======================================================"

if [[ "${DRY_RUN}" == true ]]; then
    printf '[DRY RUN] '
    printf '%q ' "${SBATCH_CMD[@]}"
    printf '\n'
else
    "${SBATCH_CMD[@]}"
fi
