#!/bin/bash
#SBATCH --job-name=FACE-facedev
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
# task -> (kf, ate_deviation, n_deviated): negative-transfer gradient at K=2 (<=2 bad sources)
KF=(5 5 5 10 10 10)
DEV=(1 2 3 1 2 3)
ND=(1 2 2 1 2 2)
i=$(( SLURM_ARRAY_TASK_ID - 1 ))
echo "=== face_dev | task ${SLURM_ARRAY_TASK_ID} | kf=${KF[$i]} ate_dev=${DEV[$i]} n_dev=${ND[$i]} | $(date) ==="
Rscript --vanilla diagnosis/face_probe/face_probe.R "${KF[$i]}" "${DEV[$i]}" "${ND[$i]}"
echo "=== done $(date) ==="
