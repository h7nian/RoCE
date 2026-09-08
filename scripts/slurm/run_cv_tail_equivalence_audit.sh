#!/bin/bash
#SBATCH --job-name=roce_cv_tail_qc
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_cv_tail_qc.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_cv_tail_qc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
BASELINE_RESULT="${ROCE_CV_TAIL_BASELINE_RESULT:?baseline result is required}"
CANDIDATE_RESULT="${ROCE_CV_TAIL_CANDIDATE_RESULT:?candidate result is required}"
BASELINE_SENSITIVITY="${ROCE_CV_TAIL_BASELINE_SENSITIVITY:?baseline sensitivity is required}"
CANDIDATE_SENSITIVITY="${ROCE_CV_TAIL_CANDIDATE_SENSITIVITY:?candidate sensitivity is required}"
AUDIT_OUTPUT="${ROCE_CV_TAIL_AUDIT_OUTPUT:?audit output is required}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
Rscript scripts/slurm/audit_cv_tail_equivalence.R \
  "${BASELINE_RESULT}" "${CANDIDATE_RESULT}" \
  "${BASELINE_SENSITIVITY}" "${CANDIDATE_SENSITIVITY}" \
  "${AUDIT_OUTPUT}"
