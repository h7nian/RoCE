#!/bin/bash
#SBATCH --job-name=aggdbg
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8g
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/face_aggdbg_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/face_aggdbg.R
