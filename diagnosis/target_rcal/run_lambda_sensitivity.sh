#!/bin/bash
#SBATCH --job-name=roce_target_lambda
#SBATCH --time=00:20:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:?}"
Rscript diagnosis/target_rcal/run_lambda_sensitivity.R "${SLURM_ARRAY_TASK_ID:?}" \
  "${ROCE_DIAGNOSTIC_ROOT:?}"
