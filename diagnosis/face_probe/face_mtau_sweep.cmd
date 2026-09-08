#!/bin/bash
#SBATCH --job-name=face_mtau
#SBATCH --time=01:30:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=24
#SBATCH --output=diagnosis/face_probe/validation/face_mtau_sweep_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
export ROCE_NSIMS=150
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_mtau_sweep.R
