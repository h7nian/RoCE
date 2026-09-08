#!/bin/bash
#SBATCH --job-name=forest
#SBATCH --time=00:10:00
#SBATCH --mem=4g
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/forest_replot_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/forest_replot.R
