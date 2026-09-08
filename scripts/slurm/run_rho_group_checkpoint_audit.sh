#!/bin/bash
#SBATCH --job-name=roce_group_qc
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_group_qc.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_group_qc.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
RESULT_ROOT="${ROCE_RESULT_ROOT:?ROCE_RESULT_ROOT is required}"
PRIMARY_MANIFEST="${ROCE_PRIMARY_MANIFEST:?ROCE_PRIMARY_MANIFEST is required}"
CONFIGURATION="${ROCE_CHECKPOINT_CONFIG:?ROCE_CHECKPOINT_CONFIG is required}"
SOURCE_COUNT="${ROCE_CHECKPOINT_K:?ROCE_CHECKPOINT_K is required}"
EXPECTED_REPLICATIONS="${ROCE_CHECKPOINT_REPLICATIONS:?ROCE_CHECKPOINT_REPLICATIONS is required}"
PRIMARY_CUTOFF="${ROCE_PRIMARY_CUTOFF:?ROCE_PRIMARY_CUTOFF is required}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:?ROCE_PROJECT_LIB is required}"

if [[ ! "${CONFIGURATION}" =~ ^C[123]$ ]] ||
   [[ ! "${SOURCE_COUNT}" =~ ^(2|4|8)$ ]] ||
   [[ ! "${EXPECTED_REPLICATIONS}" =~ ^[1-9][0-9]*$ ]] ||
   [[ "${EXPECTED_REPLICATIONS}" -gt 500 ]]; then
  echo "invalid grouped-checkpoint configuration." >&2
  exit 1
fi
if [[ ! "${PRIMARY_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v cutoff="${PRIMARY_CUTOFF}" 'BEGIN { exit !(cutoff > 0) }'; then
  echo "ROCE_PRIMARY_CUTOFF must be one positive number." >&2
  exit 1
fi

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
MANIFEST_FINGERPRINT="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
if [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  EXPECTED_THREADS=2
else
  EXPECTED_THREADS=5
fi

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"
export ROCE_MANIFEST_FINGERPRINT="${MANIFEST_FINGERPRINT}"
export ROCE_EXPECT_NUISANCE_CV_THREADS="${EXPECTED_THREADS}"
cd "${PROJECT_ROOT}"

for rho in 0 0.5 1 1.5 2 2.5; do
  Rscript scripts/slurm/audit_direct_tate_checkpoint.R \
    "${RESULT_ROOT}" "${PRIMARY_MANIFEST}" "${CONFIGURATION}" \
    "${SOURCE_COUNT}" "${rho}" "${EXPECTED_REPLICATIONS}" \
    "${PRIMARY_CUTOFF}"
done

GATE_ROOT="${RESULT_ROOT}/group_checkpoint_gates"
GATE_PATH="${GATE_ROOT}/$(printf '%s_K%d_n%03d_passed.txt' \
  "${CONFIGURATION}" "${SOURCE_COUNT}" "${EXPECTED_REPLICATIONS}")"
if [[ -e "${GATE_PATH}" ]]; then
  echo "refusing to overwrite grouped-checkpoint gate: ${GATE_PATH}" >&2
  exit 2
fi
mkdir -p "${GATE_ROOT}"
GATE_TEMPORARY="$(mktemp "${GATE_ROOT}/.group_checkpoint.XXXXXX")"
{
  echo "rho_group_checkpoint=passed"
  echo "config=${CONFIGURATION}"
  echo "K=${SOURCE_COUNT}"
  echo "expected_replications=${EXPECTED_REPLICATIONS}"
  echo "primary_cutoff=${PRIMARY_CUTOFF}"
  echo "rho_settings_audited=6"
  echo "package_fingerprint=${PACKAGE_FINGERPRINT}"
  echo "workflow_fingerprint=${WORKFLOW_FINGERPRINT}"
  echo "manifest_fingerprint=${MANIFEST_FINGERPRINT}"
  echo "checkpoint_audit_driver_md5=$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_direct_tate_checkpoint.R" | awk '{print $1}')"
  echo "simulation_qc_policy_md5=$(md5sum "${PROJECT_ROOT}/scripts/slurm/simulation_qc_policy.R" | awk '{print $1}')"
  echo "group_checkpoint_driver_md5=$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_rho_group_checkpoint_audit.sh" | awk '{print $1}')"
} > "${GATE_TEMPORARY}"
mv "${GATE_TEMPORARY}" "${GATE_PATH}"
echo "[pass] all six rho settings passed cumulative implementation QC: ${GATE_PATH}"
