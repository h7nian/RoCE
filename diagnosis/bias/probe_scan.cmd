#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=diagnosis/bias/log/%x_%j.out
#SBATCH --error=diagnosis/bias/log/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

# ============================================================================
# RoCE bias probe scan SLURM batch
# ----------------------------------------------------------------------------
# Runs diagnosis/bias/probe_scan.R, scanning (n, p) at K=3 to quantify
# how mu_pred_ts spread depends on p / n_calib_arm.
# ============================================================================

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/bias/log diagnosis/bias/probe_out

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

N_SEEDS="${PROBE_SCAN_SEEDS:-3}"
JOB_TAG="${SLURM_JOB_ID:-manual}"

echo "======================================================"
echo " RoCE Bias Probe Scan"
echo "======================================================"
echo " Seeds/cell: ${N_SEEDS}"
echo " Job:        ${SLURM_JOB_ID}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo "======================================================"

Rscript --vanilla diagnosis/bias/probe_scan.R "${JOB_TAG}" "${N_SEEDS}"

echo
echo "======================================================"
echo " probe_scan done at $(date)"
echo "======================================================"
