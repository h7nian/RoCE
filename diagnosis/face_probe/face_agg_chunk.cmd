#!/bin/bash
#SBATCH --output=diagnosis/face_probe/validation/chunks/aggchunk_%A_%a.out
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=12g
#SBATCH --job-name=aggdiag
#SBATCH --nice=5
#SBATCH --requeue
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")
IFS=$'\t' read -r CK_RHO CK_NDEV CK_SIM_START CK_SIM_END <<< "${line}"
export CK_RHO CK_NDEV CK_SIM_START CK_SIM_END CK_K=2 CK_P=10 CK_NSITE=1000 CK_NFOLDS=5
Rscript diagnosis/face_probe/face_agg_diag_chunk.R
