#!/bin/bash
# NOTE: Despite the .cmd suffix, this is a Bash SLURM submission script.
#SBATCH --output=diagnosis/bias/log/%x_%j.out
#SBATCH --error=diagnosis/bias/log/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=5

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

N_SEEDS="${PROBE_KF_SEEDS:-3}"
JOB_TAG="${SLURM_JOB_ID:-manual}"

echo "======================================================"
echo " FACE-HD K_f Probe Scan"
echo "======================================================"
echo " Seeds/cell: ${N_SEEDS}"
echo " Job:        ${SLURM_JOB_ID}"
echo " Start:      $(date)"
echo " Cpus:       ${SLURM_CPUS_PER_TASK:-?}"
echo "======================================================"

Rscript --vanilla diagnosis/bias/probe_kf.R "${JOB_TAG}" "${N_SEEDS}"

echo
echo "======================================================"
echo " probe_kf done at $(date)"
echo "======================================================"
