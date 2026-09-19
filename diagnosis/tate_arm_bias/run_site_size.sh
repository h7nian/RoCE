#!/bin/bash
#SBATCH --job-name=roce_site_size
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=20
#SBATCH --output=diagnosis/tate_arm_bias/logs/%A_%a.out
#SBATCH --error=diagnosis/tate_arm_bias/logs/%A_%a.err
set -euo pipefail
COMPARISON_LIBRARY="${ROCE_CMP_LIB:?}"
CACHE="${ROCE_CACHE_FLAG:?}"
SITE_N="${ROCE_CMP_NSITE:-1000}"
CONFIG="${ROCE_CMP_CONFIG:-C1}"
K="${ROCE_CMP_K:-2}"
SEED="${SLURM_ARRAY_TASK_ID:?}"
module load R/4.2.2-gcc-8.2.0-vp7tyde
export ROCE_PROJECT_LIB="${COMPARISON_LIBRARY}"
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_NUISANCE_CV_THREADS=5
mkdir -p diagnosis/tate_arm_bias/logs diagnosis/tate_arm_bias/out
Rscript diagnosis/tate_arm_bias/compare_site_size.R "${CONFIG}" "${K}" "${SEED}" "${CACHE}" \
  > "diagnosis/tate_arm_bias/out/${CONFIG}_K${K}_cache${CACHE}_n${SITE_N}_seed${SEED}.txt"
