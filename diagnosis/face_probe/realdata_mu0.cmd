#!/bin/bash
#SBATCH --job-name=rhc_mu0
#SBATCH --time=02:30:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=5
#SBATCH --output=diagnosis/face_probe/validation/rhc_mu0_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS=$HOME/Rlibs_em:$HOME/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript realdata.R 5 death30 10 0 mu0wald ninsclas 5 Private
