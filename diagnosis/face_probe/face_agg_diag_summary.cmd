#!/bin/bash
#SBATCH --job-name=aggsum
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/face_agg_diag_summary_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde; export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"; Rscript diagnosis/face_probe/face_agg_diag_summary.R
