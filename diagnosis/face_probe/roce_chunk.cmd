#!/bin/bash
#SBATCH --output=diagnosis/face_probe/validation/chunks/hdchunk_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --job-name=hd_chunk
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")
IFS=$'\t' read -r CK_K CK_CFG CK_P CK_RHO CK_NDEV CK_SIM_START CK_SIM_END <<< "${line}"
export CK_K CK_CFG CK_P CK_RHO CK_NDEV CK_SIM_START CK_SIM_END CK_NSITE=1000 CK_NFOLDS=5 CK_ESTIMATE_ATE=TRUE
Rscript diagnosis/face_probe/face_sim_chunk.R
