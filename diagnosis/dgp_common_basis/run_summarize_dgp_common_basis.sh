#!/bin/bash
# HISTORY: 2026-09-07 #0002 Common working basis and X-dagger misspecification: strength pilot
# Task:    Apply the pre-registered omega selection rule to the 13 finished pilot
#          cells; writes diagnosis/out/dgp_common_basis/selection/. A few minutes.

#SBATCH --job-name=dgp_common_basis_summary
#SBATCH --time=00:30:00
#SBATCH --mem=8G
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
# ~/.Renviron overrides R_LIBS_USER, so the candidate library goes in R_LIBS.
export R_LIBS="${PROJECT_ROOT}/results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/lib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

mkdir -p "diagnosis/logs"
Rscript diagnosis/dgp_common_basis/summarize_dgp_common_basis.R
