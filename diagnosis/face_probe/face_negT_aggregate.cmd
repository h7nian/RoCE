#!/bin/bash
#SBATCH --job-name=negT_agg
#SBATCH --time=00:20:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --nice=5
#SBATCH --output=diagnosis/face_probe/validation/face_negT_aggregate_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_negT_aggregate.R
