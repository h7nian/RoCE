#!/bin/bash
#SBATCH --job-name=restyle
#SBATCH --time=00:15:00
#SBATCH --mem=6g
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/face_restyle_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
unset ROCE_FIG_P; Rscript diagnosis/face_probe/face_fig1_plot.R
for P in 50 100; do for cfg in C1 C2 C3; do export ROCE_FIG_CONFIG=$cfg ROCE_FIG_P=$P; Rscript diagnosis/face_probe/face_fig1_plot.R; done; done
