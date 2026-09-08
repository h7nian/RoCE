#!/bin/bash
#SBATCH --job-name=roce_profile
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-2%2
#SBATCH --output=results/direct_tate_mc500_b5000/logs/profile_%A_%a.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/profile_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
PROFILE_ROOT="${ROCE_PROFILE_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/profile_p100}"
STAGES=(face_mu1 naive_tate dr_tate)
STAGE="${STAGES[${SLURM_ARRAY_TASK_ID}]}"

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
Rscript scripts/slurm/profile_p100_stages.R "${STAGE}"
