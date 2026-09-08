#!/bin/bash
#SBATCH --job-name=2rc_plot
#SBATCH --time=00:15:00
#SBATCH --mem=6g
#SBATCH --cpus-per-task=1
#SBATCH --output=diagnosis/face_probe/validation/face_2round_cont_plot_%j.out
set -euo pipefail

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_2round_cont_plot.R
