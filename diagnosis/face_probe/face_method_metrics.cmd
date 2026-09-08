#!/bin/bash
#SBATCH --job-name=face_metrics
#SBATCH --time=01:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=32
#SBATCH --output=diagnosis/face_probe/validation/face_method_metrics_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
export ROCE_NSIMS=200
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_method_metrics.R
