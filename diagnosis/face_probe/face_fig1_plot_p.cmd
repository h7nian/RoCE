#!/bin/bash
#SBATCH --job-name=fig1_p
#SBATCH --time=00:15:00
#SBATCH --mem=6g
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/face_fig1_plot_p_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_negT_aggregate.R
export ROCE_FIG_P=50;  Rscript diagnosis/face_probe/face_fig1_plot.R
export ROCE_FIG_P=100; Rscript diagnosis/face_probe/face_fig1_plot.R
