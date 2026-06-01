#!/bin/bash
#SBATCH --job-name=FACE-facerepro
#SBATCH --output=diagnosis/face_probe/log/%x_%j.out
#SBATCH --error=diagnosis/face_probe/log/%x_%j.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/face_probe/log
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
echo "=== face repro | $(date) ==="
Rscript --vanilla diagnosis/face_probe/repro.R
echo "=== done $(date) ==="
