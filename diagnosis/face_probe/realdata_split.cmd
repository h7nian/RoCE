#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --mem=16g
#SBATCH --cpus-per-task=5
#SBATCH --output=diagnosis/face_probe/validation/rhc_split_%x_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS=$HOME/Rlibs_em:$HOME/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
echo "=== split: K=$1 site_var=$2 ==="
Rscript realdata.R "$1" death30 10 1 "spl$2" "$2" 5
