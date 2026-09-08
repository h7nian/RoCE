#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=diagnosis/bias/log/%x_%j.out
#SBATCH --error=diagnosis/bias/log/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

# ============================================================================
# RoCE bias probe SLURM batch
# ----------------------------------------------------------------------------
# Runs diagnosis/bias/probe.R for ONE phase (A / B / C / ALL) at the
# setting supplied by probe.sh via PROBE_PHASE/N/K/P/SEEDS env vars.
# Intentionally NO requeue / NO signal — this is short-running and we want
# any failure to fail fast so we can iterate.
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

PHASE="${PROBE_PHASE:-ALL}"
N_TOTAL="${PROBE_N:-1000}"
K="${PROBE_K:-3}"
P="${PROBE_P:-20}"
N_SEEDS="${PROBE_SEEDS:-10}"

echo "======================================================"
echo " RoCE Bias Probe | phase=${PHASE}"
echo "======================================================"
echo " Setting:    n=${N_TOTAL}  K=${K}  p=${P}  seeds=${N_SEEDS}"
echo " Job:        ${SLURM_JOB_ID}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo "======================================================"

Rscript --vanilla diagnosis/bias/probe.R \
    "${PHASE}" "${N_TOTAL}" "${K}" "${P}" "${N_SEEDS}"

echo
echo "======================================================"
echo " Probe ${PHASE} done at $(date)"
echo "======================================================"
