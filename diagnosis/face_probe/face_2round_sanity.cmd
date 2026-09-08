#!/bin/bash
#SBATCH --job-name=2r_sanity
#SBATCH --time=01:30:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=4
#SBATCH --output=diagnosis/face_probe/validation/face_2round_sanity_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=$HOME/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/face_2round_sanity.R
