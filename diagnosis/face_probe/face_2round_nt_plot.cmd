#!/bin/bash
#SBATCH --job-name=2rnt_plot
#SBATCH --time=00:15:00
#SBATCH --mem=6g
#SBATCH --output=diagnosis/face_probe/validation/face_2round_nt_plot_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=$HOME/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/face_2round_nt_plot.R
