#!/bin/bash
#SBATCH --job-name=roce_roundmode
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=20
#SBATCH --output=diagnosis/tate_arm_bias/logs/%A_%a.out
#SBATCH --error=diagnosis/tate_arm_bias/logs/%A_%a.err
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${ROCE_CMP_LIB:?}"
export ROCE_NUISANCE_CV_THREADS=5
mkdir -p diagnosis/tate_arm_bias/logs diagnosis/tate_arm_bias/out
Rscript diagnosis/tate_arm_bias/compare_round_mode.R "${ROCE_CMP_CONFIG:?}" "${ROCE_CMP_K:?}" \
  "${SLURM_ARRAY_TASK_ID}" "${ROCE_ROUND_MODE:?}" \
  > "diagnosis/tate_arm_bias/out/${ROCE_ROUND_MODE}_${ROCE_CMP_CONFIG}_K${ROCE_CMP_K}_seed${SLURM_ARRAY_TASK_ID}.txt"
