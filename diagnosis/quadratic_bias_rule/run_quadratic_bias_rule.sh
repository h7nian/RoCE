#!/bin/bash
# HISTORY: 2026-09-07 #0007 Land the smooth quadratic-bias weight rule as the sensitivity estimator
# Task:    Replay the 100 saved v19 seeds with the in-package quadratic-bias rule and
#          compare with the Stage-1 candidate rows (results/.../quadratic_*_n100_v1);
#          writes diagnosis/out/quadratic_bias_rule/v1. About 5 min, single core.

#SBATCH --job-name=quadratic_bias_rule
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --partition=msismall
#SBATCH --output=diagnosis/logs/%x-%j.out
#SBATCH --error=diagnosis/logs/%x-%j.err
#SBATCH --requeue
#SBATCH --signal=B:USR1@60

set -euo pipefail

# sbatch copies the script to a spool directory, so locate the project by the
# submission directory (submit from the project root).
PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"
test -f DESCRIPTION || { echo "submit from the project root (DESCRIPTION not found in ${PROJECT_ROOT})" >&2; exit 1; }

module load R/4.2.2-gcc-8.2.0-vp7tyde
# Uses the package installed by ./test.sh into R_LIBS_USER (~/.Renviron).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

TASK="quadratic_bias_rule"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"

Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}"
