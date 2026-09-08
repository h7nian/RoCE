#!/bin/bash
# HISTORY: 2026-09-07 #0002 Common working basis and X-dagger misspecification: strength pilot
# Task:    13-cell array (C1; C2/C3/C4 x omega in {0.25,0.5,0.75,1}), one seed each,
#          p = 100, K = 2, rho = 0, full one-round TATE fit; writes
#          diagnosis/out/dgp_common_basis/<config>_omega<omega>/. About 30 min per cell
#          on 10 cores (2 source workers x 5 nuisance CV threads).

#SBATCH --job-name=dgp_common_basis
#SBATCH --array=1-13
#SBATCH --time=02:00:00
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
# Package library byte-identical to the current workspace source (built 2026-09-07).
# ~/.Renviron overrides R_LIBS_USER, so the candidate library goes in R_LIBS.
export R_LIBS="${PROJECT_ROOT}/results/direct_tate_mc500_b5000/outcome_cv_scale_candidate_v1/lib"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export ROCE_NUISANCE_CV_THREADS=5

TASK="dgp_common_basis"
mkdir -p "diagnosis/logs" "diagnosis/out/${TASK}"

Rscript "diagnosis/${TASK}/${TASK}.R" "diagnosis/out/${TASK}" "${SLURM_ARRAY_TASK_ID}"
