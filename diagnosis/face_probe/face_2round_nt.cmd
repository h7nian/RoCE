#!/bin/bash
#SBATCH --job-name=2r_nt
#SBATCH --time=02:00:00
#SBATCH --mem=8g
#SBATCH --cpus-per-task=4
#SBATCH --output=diagnosis/face_probe/validation/chunks2r_nt/log_%A_%a.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=$HOME/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" diagnosis/face_probe/em_tasks_nt.txt)
read RHO P S0 S1 <<< "$line"
export CK_RHO=$RHO CK_NDEV=1 CK_EM=0 CK_P=$P CK_SIM_START=$S0 CK_SIM_END=$S1 CK_K=4 CK_CFG=C1 CK_NFOLDS=5 CK_OUTCOME=binary
Rscript diagnosis/face_probe/face_2round_chunk.R
