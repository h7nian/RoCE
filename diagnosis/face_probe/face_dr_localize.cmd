#!/bin/bash
#SBATCH --job-name=face_dr_loc
#SBATCH --time=01:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/face_dr_localize_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_dr_localize.R
