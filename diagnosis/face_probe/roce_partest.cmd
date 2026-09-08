#!/bin/bash
#SBATCH --job-name=hd_par
#SBATCH --time=01:10:00
#SBATCH --mem=24g
#SBATCH --cpus-per-task=16
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/roce_partest_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/roce_partest.R
