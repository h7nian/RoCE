#!/bin/bash
#SBATCH --job-name=negT_K8
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --output=diagnosis/face_probe/validation/face_negT_C1_K8_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=~/Rlibs
export ROCE_NSIMS=300
export ROCE_PROBE_CONFIG=C1
export ROCE_PROBE_P=10
export ROCE_PROBE_K=8
export ROCE_PROBE_NSITE=1000
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
Rscript diagnosis/face_probe/face_negtransfer_grid.R
