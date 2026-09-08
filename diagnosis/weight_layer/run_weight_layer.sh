#!/bin/bash
# HISTORY: 2026-09-07 #0001 Weight-layer influence function for the soft-threshold aggregation rule
# Task:    Finite-difference derivative checks on the seed 10013 v19 bundle, then a
#          replay over the 100 saved v19 seeds (600 fits); writes
#          diagnosis/out/weight_layer/v1. Expected wall-time about 1-2 h, single core.

#SBATCH --job-name=weight_layer
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --partition=msismall
#SBATCH --output=diagnosis/logs/%x-%j.out
#SBATCH --error=diagnosis/logs/%x-%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@60

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${PROJECT_ROOT}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
# Package library byte-identical to the current workspace source (built 2026-09-07).
export R_LIBS_USER="${PROJECT_ROOT}/results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/lib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

TASK="weight_layer"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"

Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}"
