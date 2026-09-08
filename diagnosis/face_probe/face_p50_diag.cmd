#!/bin/bash
#SBATCH --job-name=RoCE-p50diag
#SBATCH --output=diagnosis/face_probe/log/%x_%j.out
#SBATCH --error=diagnosis/face_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12g
#SBATCH --time=0:40:00
#SBATCH -p preempt,msismall,amdsmall,agsmall,amdlarge
#SBATCH --nice=5
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/face_probe/log
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
echo "=== p50 face-DGP smoke (current source via load_all) | $(date) ==="
Rscript --vanilla diagnosis/face_probe/face_p50_diag.R
echo "=== done $(date) ==="
