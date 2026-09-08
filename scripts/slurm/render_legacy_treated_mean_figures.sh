#!/bin/bash
#SBATCH --job-name=roce_legacy_figs
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_v1/logs/%j_legacy_figures.out
#SBATCH --error=results/direct_tate_v1/logs/%j_legacy_figures.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
SUMMARY_ROOT="${PROJECT_ROOT}/diagnosis/face_probe/validation"
ARCHIVE_ROOT="${PROJECT_ROOT}/results/direct_tate_v1/legacy_figures"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_RESULTS_DIR="${SUMMARY_ROOT}"
export ROCE_FIG_ESTIMAND_SCOPE=legacy

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_v1/logs" "${ARCHIVE_ROOT}"

for dimension in 50 100; do
  export ROCE_FIG_P="${dimension}"
  for config in C1 C2 C3; do
    export ROCE_FIG_CONFIG="${config}"
    Rscript diagnosis/face_probe/face_fig1_plot.R

    generated="${SUMMARY_ROOT}/fig1_negtransfer_${config}_p${dimension}.pdf"
    manuscript_name="sim_negtransfer_${config}_p${dimension}.pdf"
    archive_name="sim_negtransfer_${config}_p${dimension}_treated_mean_legacy.pdf"
    test -s "${generated}"

    install -m 0644 "${generated}" "${ARCHIVE_ROOT}/${archive_name}"
    install -m 0644 "${generated}" \
      "${PROJECT_ROOT}/docs/figures/${manuscript_name}"
    install -m 0644 "${generated}" \
      "${PROJECT_ROOT}/overleaf/figures/${manuscript_name}"
  done
done

echo "Rendered six estimand-consistent legacy treated-mean figures."
