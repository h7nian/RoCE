#!/bin/bash
#SBATCH --job-name=roce_rhc_figure
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rhc_figure.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rhc_figure.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
RESULT_ROOT="${ROCE_OUTPUT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/rhc_supported_primary}"

if [[ ! -s "${RESULT_ROOT}/rhc_direct_tate_methods.csv" ]]; then
  echo "audited RHC method table is missing: ${RESULT_ROOT}/rhc_direct_tate_methods.csv" >&2
  exit 1
fi
mapfile -t RHC_AUDIT_GATES < <(
  find "${RESULT_ROOT}/audits" -type f \
    -name 'rhc_direct_tate_audit_passed.txt' -print 2>/dev/null | sort
)
if [[ "${#RHC_AUDIT_GATES[@]}" -ne 1 ]] ||
   ! grep -Fqx -- 'rhc_direct_tate_audit=passed' "${RHC_AUDIT_GATES[0]:-}"; then
  echo "exactly one passed RHC audit is required before rendering: ${RESULT_ROOT}" >&2
  exit 1
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

Rscript scripts/slurm/render_rhc_direct_tate_figure.R \
  "${RESULT_ROOT}/rhc_direct_tate_methods.csv" \
  "${RESULT_ROOT}/rhc_direct_tate_forest_highlighted.pdf" \
  "${PROJECT_ROOT}/docs/figures/rhc_direct_tate_forest.pdf" \
  "${PROJECT_ROOT}/overleaf/figures/rhc_direct_tate_forest.pdf"
