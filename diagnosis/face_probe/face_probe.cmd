#!/bin/bash
#SBATCH --job-name=FACE-faceprobe
#SBATCH --output=diagnosis/face_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/face_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --nice=10
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/face_probe/log diagnosis/face_probe/results
if command -v module >/dev/null 2>&1; then module load R/4.2.2-gcc-8.2.0-vp7tyde 2>/dev/null || module load R; fi
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
KFS=(5 10)
KF="${KFS[$(( SLURM_ARRAY_TASK_ID - 1 ))]}"
echo "=== face_probe | task ${SLURM_ARRAY_TASK_ID} | kf=${KF} | $(date) ==="
Rscript --vanilla diagnosis/face_probe/face_probe.R "${KF}"
echo "=== done $(date) ==="
