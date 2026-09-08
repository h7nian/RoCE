#!/bin/bash
#SBATCH --job-name=FACE-method-audit
#SBATCH --output=diagnosis/method_alignment/audit_results/%x_%j.out
#SBATCH --error=diagnosis/method_alignment/audit_results/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8g
#SBATCH --time=00:30:00
#SBATCH -p msismall,msilarge,msilong,amdsmall,agsmall,amdlarge,amd512,amd2tb
#SBATCH --nice=5

set -euo pipefail

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
export METHOD_AUDIT_OUTPUT_ROOT="${METHOD_AUDIT_OUTPUT_ROOT:-diagnosis/method_alignment/audit_results}"
mkdir -p "${METHOD_AUDIT_OUTPUT_ROOT}"

if command -v module >/dev/null 2>&1; then
    module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R
fi

export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1

echo "======================================================"
echo " RoCE Method Alignment Audit"
echo "======================================================"
echo " Repo:       $(pwd)"
echo " Output:     ${METHOD_AUDIT_OUTPUT_ROOT}"
echo " Start:      $(date)"
echo " Host:       $(hostname)"
echo " Job:        ${SLURM_JOB_ID:-NA}"
echo "======================================================"

Rscript --vanilla diagnosis/method_alignment/method_alignment_audit.R .

echo "======================================================"
echo " Method alignment audit completed at $(date)"
echo "======================================================"
