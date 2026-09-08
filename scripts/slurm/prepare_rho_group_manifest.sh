#!/bin/bash
#SBATCH --job-name=roce_rho_manifest
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_rho_manifest.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_rho_manifest.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${ROCE_MANIFEST_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
PRIMARY_MANIFEST="${MANIFEST_ROOT}/manifest_main.csv"
GROUP_MANIFEST="${MANIFEST_ROOT}/manifest_main_rho_groups.csv"

if [[ ! -f "${PRIMARY_MANIFEST}" ]]; then
  echo "primary manifest is missing: ${PRIMARY_MANIFEST}" >&2
  exit 1
fi
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"

if [[ -e "${GROUP_MANIFEST}" ]]; then
  echo "[reuse] auditing existing grouped manifest: ${GROUP_MANIFEST}"
else
  Rscript scripts/slurm/build_direct_tate_rho_group_manifest.R \
    "${PRIMARY_MANIFEST}" "${GROUP_MANIFEST}"
fi
Rscript scripts/slurm/audit_mc500_manifests.R \
  "${MANIFEST_ROOT}" 500 5000
