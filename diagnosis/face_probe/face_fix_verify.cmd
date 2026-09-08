#!/bin/bash
#SBATCH --job-name=RoCE-fixverify
#SBATCH --output=diagnosis/face_probe/log/%x_%j.out
#SBATCH --error=diagnosis/face_probe/log/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16g
#SBATCH --time=2:00:00
#SBATCH -p preempt,msismall,amdsmall,agsmall,amdlarge,msilarge
#SBATCH --nice=6
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/face_probe/log
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
echo "=== fix verify (load_all current source) | $(date) ==="
Rscript --vanilla diagnosis/face_probe/face_fix_verify.R
echo "=== done $(date) ==="
