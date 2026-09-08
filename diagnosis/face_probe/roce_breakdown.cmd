#!/bin/bash
#SBATCH --job-name=hd_break
#SBATCH --time=01:30:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=4
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/roce_breakdown_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/roce_breakdown.R
