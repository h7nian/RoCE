#!/bin/bash
#SBATCH --job-name=roce_manifest
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_manifest.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_manifest.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${ROCE_MANIFEST_ROOT:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
N_REPLICATIONS="${ROCE_N_REPLICATIONS:-500}"
CUTOFF_SELECTION_GATE="${ROCE_CUTOFF_SELECTION_GATE:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot/cutoff_decision_n010/cutoff_selection_passed.txt}"

if [[ ! -f "${CUTOFF_SELECTION_GATE}" ]]; then
  echo "audited cutoff-selection gate is missing: ${CUTOFF_SELECTION_GATE}" >&2
  exit 1
fi
read_gate_value() {
  local key="$1"
  local value
  value="$(awk -F= -v key="${key}" '$1 == key { print substr($0, length(key) + 2) }' \
    "${CUTOFF_SELECTION_GATE}")"
  if [[ -z "${value}" || "$(printf '%s\n' "${value}" | wc -l)" -ne 1 ]]; then
    echo "cutoff-selection gate has no unique ${key}." >&2
    exit 1
  fi
  printf '%s\n' "${value}"
}
if [[ "$(read_gate_value cutoff_selection)" != "passed" ]]; then
  echo "cutoff-selection gate did not pass." >&2
  exit 1
fi
ROCE_PRIMARY_CUTOFF="$(read_gate_value selected_cutoff)"
if [[ ! "${ROCE_PRIMARY_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v cutoff="${ROCE_PRIMARY_CUTOFF}" \
      'BEGIN { exit !(cutoff > 0) }'; then
  echo "selected cutoff is not a positive number." >&2
  exit 1
fi
ROCE_CUTOFF_SELECTION_FINGERPRINT="$(sha256sum \
  "${CUTOFF_SELECTION_GATE}" | awk '{print $1}')"
export ROCE_PRIMARY_CUTOFF
export ROCE_CUTOFF_SELECTION_FINGERPRINT
export ROCE_N_FOLDS="${ROCE_N_FOLDS:-5}"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"

cd "${PROJECT_ROOT}"
mkdir -p "${MANIFEST_ROOT}" "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
echo "[cutoff] selected=${ROCE_PRIMARY_CUTOFF} gate_sha256=${ROCE_CUTOFF_SELECTION_FINGERPRINT}"

Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${MANIFEST_ROOT}/manifest_main.csv" main "${N_REPLICATIONS}"
Rscript scripts/slurm/split_direct_tate_manifest_by_k.R \
  "${MANIFEST_ROOT}/manifest_main.csv" "${MANIFEST_ROOT}"
Rscript scripts/slurm/build_direct_tate_rho_group_manifest.R \
  "${MANIFEST_ROOT}/manifest_main.csv" \
  "${MANIFEST_ROOT}/manifest_main_rho_groups.csv"
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${MANIFEST_ROOT}/manifest_smoke_single.csv" smoke 1
Rscript scripts/slurm/build_direct_tate_cutoff_manifest.R \
  "${MANIFEST_ROOT}/manifest_cutoff_diagnostic.csv" "${N_REPLICATIONS}"
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${MANIFEST_ROOT}/manifest_truncation_diagnostic.csv" truncation \
  "${N_REPLICATIONS}"
Rscript scripts/slurm/audit_mc500_manifests.R \
  "${MANIFEST_ROOT}" "${N_REPLICATIONS}" 5000 main

# Shared-shift family (HISTORY #0009): its own root with the same file names,
# so the grouped submission and checkpoint tooling apply unchanged.
SHARED_SHIFT_ROOT="${MANIFEST_ROOT}/shared_shift"
mkdir -p "${SHARED_SHIFT_ROOT}"
Rscript scripts/slurm/build_direct_tate_manifest.R \
  "${SHARED_SHIFT_ROOT}/manifest_main.csv" shared_shift "${N_REPLICATIONS}"
Rscript scripts/slurm/split_direct_tate_manifest_by_k.R \
  "${SHARED_SHIFT_ROOT}/manifest_main.csv" "${SHARED_SHIFT_ROOT}"
Rscript scripts/slurm/build_direct_tate_rho_group_manifest.R \
  "${SHARED_SHIFT_ROOT}/manifest_main.csv" \
  "${SHARED_SHIFT_ROOT}/manifest_main_rho_groups.csv"
Rscript scripts/slurm/audit_mc500_manifests.R \
  "${SHARED_SHIFT_ROOT}" "${N_REPLICATIONS}" 5000 shared_shift
