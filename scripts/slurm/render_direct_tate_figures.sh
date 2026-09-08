#!/bin/bash
#SBATCH --job-name=roce_figures
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_figures.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_figures.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${ROCE_RESULT_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
MANIFEST_ROOT="${ROCE_MANIFEST_ROOT:-${RESULT_ROOT}}"
SUMMARY_ROOT="${RESULT_ROOT}/summary"
LEGACY_FIGURE_ROOT="${PROJECT_ROOT}/results/direct_tate_v1/legacy_figures"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${RESULT_ROOT}/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
PRIMARY_MANIFEST="${MANIFEST_ROOT}/manifest_main.csv"
if [[ ! -f "${PRIMARY_MANIFEST}" ]]; then
  echo "primary simulation manifest is missing: ${PRIMARY_MANIFEST}" >&2
  exit 1
fi
MANIFEST_FINGERPRINT="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_MANIFEST_FINGERPRINT="${MANIFEST_FINGERPRINT}"
export ROCE_RESULTS_DIR="${SUMMARY_ROOT}"
export ROCE_FIG_P=100
export ROCE_FIG_ESTIMAND_SCOPE=direct

cd "${PROJECT_ROOT}"
mkdir -p "${RESULT_ROOT}/logs" "${LEGACY_FIGURE_ROOT}"

Rscript scripts/slurm/aggregate_direct_tate.R "${RESULT_ROOT}" 500
Rscript scripts/slurm/diagnose_direct_tate_settings.R \
  "${RESULT_ROOT}" 500 "${PRIMARY_MANIFEST}" TRUE

for config in C1 C2 C3; do
  export ROCE_FIG_CONFIG="${config}"
  Rscript diagnosis/face_probe/face_fig1_plot.R

  generated="${SUMMARY_ROOT}/fig1_negtransfer_${config}_p100.pdf"
  legacy_manuscript_name="sim_negtransfer_${config}_p100.pdf"
  manuscript_name="sim_negtransfer_${config}_p100_direct_tate.pdf"
  legacy_name="sim_negtransfer_${config}_p100_armwise_legacy.pdf"
  test -s "${generated}"

  # Preserve the pre-TATE p=100 panel exactly once. TATE output
  # receives its own filename, so rendering never overwrites an earlier result.
  if [[ ! -s "${LEGACY_FIGURE_ROOT}/${legacy_name}" ]]; then
    test -s "${PROJECT_ROOT}/docs/figures/${legacy_manuscript_name}"
    install -m 0644 "${PROJECT_ROOT}/docs/figures/${legacy_manuscript_name}" \
      "${LEGACY_FIGURE_ROOT}/${legacy_name}"
    install -m 0644 "${PROJECT_ROOT}/docs/figures/${legacy_manuscript_name}" \
      "${PROJECT_ROOT}/docs/figures/${legacy_name}"
    install -m 0644 "${PROJECT_ROOT}/docs/figures/${legacy_manuscript_name}" \
      "${PROJECT_ROOT}/overleaf/figures/${legacy_name}"
  fi

  install -m 0644 "${generated}" \
    "${PROJECT_ROOT}/docs/figures/${manuscript_name}"
  install -m 0644 "${generated}" \
    "${PROJECT_ROOT}/overleaf/figures/${manuscript_name}"
done

echo "Rendered TATE p=100 figures with enlarged manuscript fonts."
