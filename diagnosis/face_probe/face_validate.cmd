#!/bin/bash
#SBATCH --job-name=RoCE-val
#SBATCH --output=diagnosis/face_probe/log/%x_%A_%a.out
#SBATCH --error=diagnosis/face_probe/log/%x_%A_%a.err
#SBATCH --open-mode=append
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24g
#SBATCH --time=12:00:00
#SBATCH -p preempt,msismall,amdsmall,agsmall,amdlarge,msilarge,msilong
#SBATCH --nice=8
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p diagnosis/face_probe/log diagnosis/face_probe/validation
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export ROCE_VAL_NSIMS="${ROCE_VAL_NSIMS:-50}"

CELLS="diagnosis/face_probe/validation/cells.txt"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CELLS")
echo "=== val task ${SLURM_ARRAY_TASK_ID}: [$line] nsims=${ROCE_VAL_NSIMS} | $(date) ==="
# shellcheck disable=SC2086
Rscript --vanilla diagnosis/face_probe/face_validate.R $line
echo "=== done $(date) ==="
