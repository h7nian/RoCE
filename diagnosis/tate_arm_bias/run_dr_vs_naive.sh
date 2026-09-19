#!/bin/bash
#SBATCH --job-name=roce_drnaive
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=diagnosis/tate_arm_bias/logs/%A_%a.out
#SBATCH --error=diagnosis/tate_arm_bias/logs/%A_%a.err
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export ROCE_PROJECT_LIB="${ROCE_CMP_LIB:?}"
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
mkdir -p diagnosis/tate_arm_bias/logs diagnosis/tate_arm_bias/out
Rscript diagnosis/tate_arm_bias/dr_vs_naive.R "${ROCE_CMP_CONFIG:-C3}" \
  "${SLURM_ARRAY_TASK_ID}" "${ROCE_CMP_NSITE:-1000}" \
  > "diagnosis/tate_arm_bias/out/${ROCE_CMP_CONFIG:-C3}_n${ROCE_CMP_NSITE:-1000}_seed${SLURM_ARRAY_TASK_ID}.txt"
