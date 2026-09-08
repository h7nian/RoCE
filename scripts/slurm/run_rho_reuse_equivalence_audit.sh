#!/bin/bash
#SBATCH --job-name=roce_rho_reuse_qc
#SBATCH --time=00:15:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rho_reuse_qc.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rho_reuse_qc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
INDEPENDENT_RESULT="${ROCE_INDEPENDENT_RESULT:?ROCE_INDEPENDENT_RESULT is required}"
GROUPED_RESULT="${ROCE_GROUPED_RESULT:?ROCE_GROUPED_RESULT is required}"
AUDIT_OUTPUT="${ROCE_REUSE_AUDIT_OUTPUT:?ROCE_REUSE_AUDIT_OUTPUT is required}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
cd "${PROJECT_ROOT}"

Rscript scripts/slurm/audit_rho_reuse_equivalence.R \
  "${INDEPENDENT_RESULT}" "${GROUPED_RESULT}" "${AUDIT_OUTPUT}"
