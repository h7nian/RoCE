#!/bin/bash
#SBATCH --job-name=face_bias_se
#SBATCH --time=03:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=diagnosis/face_probe/validation/face_bias_se_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
export ROCE_NSIMS=100
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_bias_se.R
