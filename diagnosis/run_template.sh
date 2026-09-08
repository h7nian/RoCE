#!/bin/bash
# diagnosis/<task>/run_<task>.sh — sbatch wrapper for a Stage-1 prototype
# (iterate-with-discipline Rules 16, 16a, 17, 19; R-project amendment A2/A4
# recorded in diagnosis/HISTORY.md #0001).
#
# Copy this file to diagnosis/<task>/run_<task>.sh, fill every <FILL>, and
# submit with `sbatch diagnosis/<task>/run_<task>.sh` from the project root.
# Logs land in diagnosis/logs/, results in diagnosis/out/<task>/.

# HISTORY: YYYY-MM-DD <one-line title of the iteration this serves>
# Task:    <1-2 lines: inputs, outputs, expected wall-time>

#SBATCH --job-name=<task>
#SBATCH --time=<FILL>
#SBATCH --mem=<FILL>
#SBATCH --cpus-per-task=<FILL>
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
# ~/.Renviron pins R_LIBS_USER; put any candidate library in R_LIBS instead.
export R_LIBS_USER=/users/0/zhan9381/Rlibs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

TASK="<task>"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"

Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}"
