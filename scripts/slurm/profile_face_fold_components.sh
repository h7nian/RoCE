#!/bin/bash
#SBATCH --job-name=roce_foldprof
#SBATCH --partition=preempt
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/foldprof_%j.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/foldprof_%j.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
PROFILE_ROOT="${ROCE_PROFILE_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/profile_face_fold_components}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PROFILE_ROOT="${PROFILE_ROOT}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "${PROJECT_ROOT}"
mkdir -p "${PROFILE_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript scripts/slurm/profile_face_fold_components.R
