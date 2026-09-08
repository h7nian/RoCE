#!/bin/bash
#SBATCH --job-name=face_tests_dgp
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=diagnosis/face_probe/validation/face_tests_dgp_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
./test.sh --source --filter 'data_generation|per-site' --reporter progress
