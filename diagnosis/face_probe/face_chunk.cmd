#!/bin/bash
#SBATCH --output=diagnosis/face_probe/validation/chunks/chunk_%A_%a.out
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12g
#SBATCH --job-name=negT_chunk
#SBATCH --nice=5
#SBATCH --requeue
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="${HOME}/Rlibs"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:-$(git rev-parse --show-toplevel)}"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${COMBO_FILE}")
IFS=$'\t' read -r CK_K CK_CFG CK_P CK_RHO CK_NDEV CK_SIM_START CK_SIM_END <<< "${line}"
export CK_K CK_CFG CK_P CK_RHO CK_NDEV CK_SIM_START CK_SIM_END
export CK_NSITE=1000 CK_NFOLDS=5 CK_ESTIMATE_ATE=TRUE
echo "[task ${SLURM_ARRAY_TASK_ID}] K=${CK_K} ${CK_CFG} p=${CK_P} rho=${CK_RHO} nd=${CK_NDEV} sims=${CK_SIM_START}-${CK_SIM_END}"
Rscript diagnosis/face_probe/face_sim_chunk.R
