#!/bin/bash
#SBATCH --job-name=fig2_w
#SBATCH --time=02:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=24
#SBATCH --output=diagnosis/face_probe/validation/face_fig2_weights_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
export ROCE_NSIMS=200
export ROCE_PROBE_K=2
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_fig2_weights.R
