#!/bin/bash
# HISTORY: 2026-09-07 #0006 Land the common-basis DGP with X-dagger misspecification
# Task:    4-cell array (C1, C2, C3, C4 at the pre-registered strength), seed 20001,
#          p = 100, K = 2, rho = 0: in-package DGP versus the #0002 prototype cells
#          (standardization, calibration, truth) and a full refit comparison; writes
#          diagnosis/out/dgp_common_basis/package_<config>/. About 1.5 h per cell on
#          10 cores (2 source workers x 5 nuisance CV threads).

#SBATCH --job-name=dgp_common_basis
#SBATCH --array=1-4
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=10
#SBATCH --partition=msismall
#SBATCH --output=diagnosis/logs/%x-%A_%a.out
#SBATCH --error=diagnosis/logs/%x-%A_%a.err
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
export ROCE_NUISANCE_CV_THREADS=5

TASK="dgp_common_basis"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"

Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}" "${SLURM_ARRAY_TASK_ID}"
