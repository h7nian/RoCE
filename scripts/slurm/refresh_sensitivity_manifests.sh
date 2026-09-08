#!/bin/bash
#SBATCH --job-name=roce_sens_manifest
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_sens_manifest.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_sens_manifest.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${ROCE_MANIFEST_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
ARCHIVE_ROOT="${ROCE_MANIFEST_ARCHIVE_ROOT:-${MANIFEST_ROOT}/manifest_archive_pre_group_20260814}"
CUTOFF_MANIFEST="${MANIFEST_ROOT}/manifest_cutoff_diagnostic.csv"
TRUNCATION_MANIFEST="${MANIFEST_ROOT}/manifest_truncation_diagnostic.csv"

if [[ -e "${ARCHIVE_ROOT}" ]]; then
  echo "refusing to reuse sensitivity-manifest archive: ${ARCHIVE_ROOT}" >&2
  exit 2
fi
if [[ ! -f "${CUTOFF_MANIFEST}" || ! -f "${TRUNCATION_MANIFEST}" ]]; then
  echo "both legacy sensitivity manifests must exist before archival." >&2
  exit 1
fi
if find "${MANIFEST_ROOT}" -maxdepth 3 -type f \
    \( -path '*cutoff*/*task_*.csv' -o \
       -path '*truncation*/*task_*.csv' \) -print -quit | grep -q .; then
  echo "sensitivity task outputs exist; refusing to remap their task IDs." >&2
  exit 1
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
cd "${PROJECT_ROOT}"
mkdir -p "${ARCHIVE_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
mv "${CUTOFF_MANIFEST}" "${ARCHIVE_ROOT}/"
mv "${TRUNCATION_MANIFEST}" "${ARCHIVE_ROOT}/"

Rscript scripts/slurm/build_direct_tate_cutoff_manifest.R \
  "${CUTOFF_MANIFEST}" 500
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${TRUNCATION_MANIFEST}" truncation 500
Rscript scripts/slurm/audit_mc500_manifests.R \
  "${MANIFEST_ROOT}" 500 5000

echo "[done] archived legacy sensitivity manifests under ${ARCHIVE_ROOT}"
