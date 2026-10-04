#!/bin/bash
#SBATCH --job-name=roce_rcal_audit
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
set -euo pipefail
PROJECT_ROOT="${ROCE_PROJECT_ROOT:?Set ROCE_PROJECT_ROOT}"
RESULT_ROOT="${ROCE_DIAGNOSTIC_ROOT:?Set ROCE_DIAGNOSTIC_ROOT}"
export ROCE_PROJECT_LIB="${ROCE_PROJECT_LIB:?Set ROCE_PROJECT_LIB}"
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_NUISANCE_CV_THREADS=1
cd "${PROJECT_ROOT}"
Rscript diagnosis/target_rcal/run_diagnostic.R "${ROCE_DIAGNOSTIC_MODE:?}" \
  "${ROCE_DIAGNOSTIC_CONFIG:?}" "${SLURM_ARRAY_TASK_ID:-1}" "${RESULT_ROOT}"
