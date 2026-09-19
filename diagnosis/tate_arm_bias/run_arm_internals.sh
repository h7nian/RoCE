#!/bin/bash
#SBATCH --job-name=roce_arminternals
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=20
#SBATCH --output=diagnosis/tate_arm_bias/logs/%A_%a.out
#SBATCH --error=diagnosis/tate_arm_bias/logs/%A_%a.err
set -euo pipefail
module load R/4.2.2-gcc-8.2.0-vp7tyde
export ROCE_PROJECT_LIB="${ROCE_CMP_LIB:?}"
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_NUISANCE_CV_THREADS=5
mkdir -p diagnosis/tate_arm_bias/logs diagnosis/tate_arm_bias/out
Rscript diagnosis/tate_arm_bias/arm_internals.R "${ROCE_CMP_CONFIG:-C3}" "${ROCE_CMP_K:-2}" \
  "${SLURM_ARRAY_TASK_ID}" > "diagnosis/tate_arm_bias/out/seed${SLURM_ARRAY_TASK_ID}.txt"
