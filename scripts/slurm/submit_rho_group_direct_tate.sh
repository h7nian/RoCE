#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
MANIFEST_ROOT="${1:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000}"
OUTPUT_ROOT="${2:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/raw}"
OUTPUT_ROOT="$(readlink -m -- "${OUTPUT_ROOT}")"
RESULT_ROOT="$(dirname "${OUTPUT_ROOT}")"
PRIMARY_MANIFEST="${MANIFEST_ROOT}/manifest_main.csv"
GROUP_MANIFEST="${MANIFEST_ROOT}/manifest_main_rho_groups.csv"
MANIFEST_GATE="${ROCE_MANIFEST_GATE:-${MANIFEST_ROOT}/manifest_audit_passed.txt}"
# The reuse-equivalence gate lives under the family's result root: each family
# (negative transfer with treated-arm reuse, shared shift with both-arm reuse)
# must have audited its own grouped reuse (HISTORY #0009).
REUSE_GATE="${ROCE_RHO_REUSE_GATE:-${RESULT_ROOT}/rho_reuse_equivalence_final/rho_reuse_equivalence_passed.txt}"
CUTOFF_SELECTION_GATE="${ROCE_CUTOFF_SELECTION_GATE:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/grouped_cutoff_pilot/cutoff_decision_n010/cutoff_selection_passed.txt}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"
PACKAGE_TEST_GATE="${ROCE_PACKAGE_TEST_GATE:-${PROJECT_LIBRARY}/audit_tests_passed.txt}"
PACKAGE_CHECK_GATE="${ROCE_PACKAGE_CHECK_GATE:-}"
SLURM_PARTITION="${ROCE_SLURM_PARTITION:-msismall}"
SENSITIVITY_OUTPUT_ROOT="${ROCE_SENSITIVITY_OUTPUT_ROOT:-${RESULT_ROOT}/reused_sensitivity_raw}"
SENSITIVITY_OUTPUT_ROOT="$(readlink -m -- "${SENSITIVITY_OUTPUT_ROOT}")"
if [[ "$(basename "${OUTPUT_ROOT}")" != "raw" ]]; then
  echo "grouped production OUTPUT_ROOT must be a canonical root/raw directory." >&2
  exit 1
fi
if [[ ! "${SLURM_PARTITION}" =~ ^[A-Za-z0-9_-]+$ ]]; then
  echo "ROCE_SLURM_PARTITION contains invalid characters." >&2
  exit 1
fi

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
source "${PROJECT_ROOT}/scripts/slurm/resource_topology.sh"

for required_path in \
  "${PRIMARY_MANIFEST}" "${GROUP_MANIFEST}" "${MANIFEST_GATE}" \
  "${REUSE_GATE}" "${CUTOFF_SELECTION_GATE}" "${PACKAGE_TEST_GATE}"; do
  if [[ ! -f "${required_path}" ]]; then
    echo "required production gate/input is missing: ${required_path}" >&2
    exit 1
  fi
done
if [[ -z "${PACKAGE_CHECK_GATE}" || ! -f "${PACKAGE_CHECK_GATE}" ]]; then
  echo "ROCE_PACKAGE_CHECK_GATE must name the final R CMD check gate." >&2
  exit 1
fi

roce_require_package_library "${PROJECT_LIBRARY}"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
PACKAGE_SOURCE_FINGERPRINT="$(roce_package_source_fingerprint "${PROJECT_ROOT}")"
TEST_SUITE_FINGERPRINT="$({
  find "${PROJECT_ROOT}/tests" -type f -name '*.R' -print0 \
    | sort -z | xargs -0 sha256sum
} | sha256sum | awk '{print $1}')"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"
PRIMARY_MANIFEST_SHA256="$(sha256sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
GROUP_MANIFEST_SHA256="$(sha256sum "${GROUP_MANIFEST}" | awk '{print $1}')"

require_gate_value() {
  local gate_path="$1"
  local key="$2"
  local expected="$3"
  if ! grep -Fqx -- "${key}=${expected}" "${gate_path}"; then
    echo "gate mismatch in ${gate_path}: expected ${key}=${expected}" >&2
    exit 1
  fi
}

require_gate_value "${MANIFEST_GATE}" manifest_audit passed
require_gate_value "${MANIFEST_GATE}" replications_per_setting 500
require_gate_value "${MANIFEST_GATE}" bootstrap_draws 5000
require_gate_value "${MANIFEST_GATE}" p 100
CUTOFF_SELECTION_FINGERPRINT="$(sha256sum \
  "${CUTOFF_SELECTION_GATE}" | awk '{print $1}')"
SELECTED_CUTOFF="$(awk -F= '$1 == "selected_cutoff" {
  print substr($0, length($1) + 2)
}' "${CUTOFF_SELECTION_GATE}")"
if [[ "$(awk -F= '$1 == "cutoff_selection" {
       print substr($0, length($1) + 2)
     }' "${CUTOFF_SELECTION_GATE}")" != "passed" ]] ||
   [[ ! "${SELECTED_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "cutoff-selection gate is malformed or did not pass." >&2
  exit 1
fi
require_gate_value \
  "${MANIFEST_GATE}" primary_cutoff "${SELECTED_CUTOFF}"
require_gate_value \
  "${MANIFEST_GATE}" cutoff_selection_fingerprint \
  "${CUTOFF_SELECTION_FINGERPRINT}"
require_gate_value \
  "${MANIFEST_GATE}" manifest_main_md5 \
  "$(md5sum "${PRIMARY_MANIFEST}" | awk '{print $1}')"
require_gate_value \
  "${MANIFEST_GATE}" manifest_main_rho_groups_md5 \
  "$(md5sum "${GROUP_MANIFEST}" | awk '{print $1}')"
require_gate_value \
  "${MANIFEST_GATE}" audit_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_mc500_manifests.R" | awk '{print $1}')"
require_gate_value \
  "${MANIFEST_GATE}" rho_group_builder_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/build_direct_tate_rho_group_manifest.R" | awk '{print $1}')"

require_gate_value "${PACKAGE_TEST_GATE}" package_tests passed
require_gate_value \
  "${PACKAGE_TEST_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" package_source_fingerprint \
  "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" test_suite_fingerprint "${TEST_SUITE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_TEST_GATE}" test_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_package_audit_tests.sh" | awk '{print $1}')"
require_gate_value "${PACKAGE_CHECK_GATE}" r_cmd_check passed
require_gate_value \
  "${PACKAGE_CHECK_GATE}" package_source_fingerprint \
  "${PACKAGE_SOURCE_FINGERPRINT}"
require_gate_value \
  "${PACKAGE_CHECK_GATE}" check_driver_md5 \
  "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_r_cmd_check.sh" | awk '{print $1}')"

# Family from the primary manifest (deviation_mechanism is its last column).
FAMILY_MECHANISM="$(awk -F, 'NR == 2 { value = $NF; gsub(/"/, "", value); print value; exit }' "${PRIMARY_MANIFEST}")"
case "${FAMILY_MECHANISM}" in
  treated_arm) REUSE_CONFIG="C3" ;;
  both_arms) REUSE_CONFIG="C1" ;;
  *) echo "primary manifest has an unknown deviation_mechanism: ${FAMILY_MECHANISM}" >&2; exit 1 ;;
esac
require_gate_value "${REUSE_GATE}" rho_reuse_equivalence passed
require_gate_value "${REUSE_GATE}" config "${REUSE_CONFIG}"
require_gate_value "${REUSE_GATE}" deviation_mechanism "${FAMILY_MECHANISM}"
require_gate_value "${REUSE_GATE}" p 100
require_gate_value "${REUSE_GATE}" K 4
require_gate_value "${REUSE_GATE}" rho 2.5
require_gate_value "${REUSE_GATE}" nlambda_init 100
require_gate_value "${REUSE_GATE}" n_bootstrap 5000
require_gate_value "${REUSE_GATE}" primary_cutoff "${SELECTED_CUTOFF}"
require_gate_value "${REUSE_GATE}" positive_rho_workers 5
require_gate_value "${REUSE_GATE}" positive_rho_backend psock
require_gate_value \
  "${REUSE_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
require_gate_value \
  "${REUSE_GATE}" workflow_fingerprint "${WORKFLOW_FINGERPRINT}"
require_gate_value \
  "${REUSE_GATE}" manifest_fingerprint "${PRIMARY_MANIFEST_SHA256}"
require_gate_value \
  "${REUSE_GATE}" group_manifest_fingerprint "${GROUP_MANIFEST_SHA256}"
for gate_file_pair in \
  "reuse_audit_driver_md5:scripts/slurm/audit_rho_reuse_equivalence.R" \
  "rho_group_task_driver_md5:scripts/slurm/run_direct_tate_rho_group_task.R" \
  "rho_group_array_driver_md5:scripts/slurm/run_direct_tate_rho_group_array.sh" \
  "task_helpers_md5:scripts/slurm/direct_tate_task_helpers.R" \
  "checkpoint_audit_driver_md5:scripts/slurm/audit_direct_tate_checkpoint.R" \
  "simulation_qc_policy_md5:scripts/slurm/simulation_qc_policy.R" \
  "group_checkpoint_driver_md5:scripts/slurm/run_rho_group_checkpoint_audit.sh" \
  "rho_group_submit_driver_md5:scripts/slurm/submit_rho_group_direct_tate.sh"; do
  gate_key="${gate_file_pair%%:*}"
  relative_path="${gate_file_pair#*:}"
  require_gate_value \
    "${REUSE_GATE}" "${gate_key}" \
    "$(md5sum "${PROJECT_ROOT}/${relative_path}" | awk '{print $1}')"
done

N_GROUPS="$(awk 'END { print NR - 1 }' "${GROUP_MANIFEST}")"
# Every grouped job covers the six rho values of one seed-setting, so the
# primary manifest must be exactly six times the group manifest (4500 groups
# for the negative-transfer family, 1500 for the shared-shift family).
N_PRIMARY="$(awk 'END { print NR - 1 }' "${PRIMARY_MANIFEST}")"
if [[ "${N_GROUPS}" -lt 1 || $((N_GROUPS * 6)) -ne "${N_PRIMARY}" ]]; then
  echo "group manifest (${N_GROUPS} rows) must partition the primary manifest (${N_PRIMARY} rows) into six-rho groups." >&2
  exit 1
fi
COMMIT_ROOT="${OUTPUT_ROOT}/rho_group_commits"
# Optional per-setting scope (HISTORY #0010): ROCE_SETTING="C1:K2" restricts
# this call to the grouped rows of one config/K block, so independent settings
# can run in parallel, each walking its own checkpoint ladder.
SETTING="${ROCE_SETTING:-}"
SETTING_CONFIG=""
SETTING_K=""
if [[ -n "${SETTING}" ]]; then
  if [[ ! "${SETTING}" =~ ^(C[123]):K(2|4|8)$ ]]; then
    echo "ROCE_SETTING must look like C1:K2." >&2
    exit 1
  fi
  SETTING_CONFIG="${BASH_REMATCH[1]}"
  SETTING_K="${BASH_REMATCH[2]}"
fi
in_setting() {
  # Row fields: group_task_id, experiment, sim_id, config, p, K, ...
  local config="$1" k="$2"
  [[ -z "${SETTING}" || ( "${config}" == "${SETTING_CONFIG}" && "${k}" == "${SETTING_K}" ) ]]
}
FIRST_MISSING=""
SETTING_LAST_ROW=""
while IFS=, read -r group_task_id _ _ config _ k _; do
  group_task_id="${group_task_id//\"/}"
  config="${config//\"/}"
  k="${k//\"/}"
  in_setting "${config}" "${k}" || continue
  SETTING_LAST_ROW="${group_task_id}"
  if [[ -z "${FIRST_MISSING}" &&
        ! -s "${COMMIT_ROOT}/$(printf 'group_%06d_committed.txt' "${group_task_id}")" ]]; then
    FIRST_MISSING="${group_task_id}"
  fi
done < <(tail -n +2 "${GROUP_MANIFEST}")
if [[ -n "${SETTING}" && -z "${SETTING_LAST_ROW}" ]]; then
  echo "setting ${SETTING} has no rows in the grouped manifest." >&2
  exit 1
fi
if [[ -z "${FIRST_MISSING}" ]]; then
  echo "all rho-group jobs${SETTING:+ of ${SETTING}} are already committed"
  exit 0
fi

BATCH_START="${ROCE_BATCH_START:-${FIRST_MISSING}}"
BATCH_SIZE="${ROCE_BATCH_SIZE:-1}"
MAX_CONCURRENT="${ROCE_MAX_CONCURRENT:-1}"
if [[ ! "${BATCH_START}" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "${BATCH_SIZE}" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "${MAX_CONCURRENT}" =~ ^[1-9][0-9]*$ ]]; then
  echo "batch start, size, and concurrency must be positive integers." >&2
  exit 1
fi
# Caps (HISTORY #0010): one call never exceeds one setting (500 grouped jobs)
# and 50 concurrent array tasks; the checkpoint ladder below bounds a call
# further.
if [[ "${BATCH_SIZE}" -gt 500 || "${MAX_CONCURRENT}" -gt 50 ]]; then
  echo "cap is 500 grouped jobs per call and concurrency 50." >&2
  exit 1
fi
if [[ "${BATCH_START}" -gt "${N_GROUPS}" ]]; then
  echo "batch start exceeds the ${N_GROUPS}-row grouped manifest." >&2
  exit 1
fi
if [[ "${BATCH_START}" -ne "${FIRST_MISSING}" ]]; then
  echo "production batches must start at the first uncommitted grouped row${SETTING:+ of ${SETTING}} (${FIRST_MISSING})." >&2
  exit 1
fi
BATCH_END=$((BATCH_START + BATCH_SIZE - 1))
if [[ "${BATCH_END}" -gt "${N_GROUPS}" ]]; then
  BATCH_END="${N_GROUPS}"
fi
if [[ -n "${SETTING}" && "${BATCH_END}" -gt "${SETTING_LAST_ROW}" ]]; then
  BATCH_END="${SETTING_LAST_ROW}"
fi

read_group_field() {
  local row_number="$1"
  local column_number="$2"
  awk -F, -v row="$((row_number + 1))" -v column="${column_number}" '
    NR == row { value = $column; gsub(/"/, "", value); print value; exit }
  ' "${GROUP_MANIFEST}"
}
SOURCE_COUNT="$(read_group_field "${BATCH_START}" 6)"
PRIMARY_CUTOFF="$(read_group_field "${BATCH_START}" 7)"
CONFIGURATION="$(read_group_field "${BATCH_START}" 4)"
START_SIM_ID="$(read_group_field "${BATCH_START}" 3)"
while [[ "${BATCH_END}" -gt "${BATCH_START}" ]]; do
  end_source_count="$(read_group_field "${BATCH_END}" 6)"
  end_configuration="$(read_group_field "${BATCH_END}" 4)"
  if [[ "${end_source_count}" == "${SOURCE_COUNT}" &&
        "${end_configuration}" == "${CONFIGURATION}" ]]; then
    break
  fi
  BATCH_END=$((BATCH_END - 1))
done
# Stop at the next predeclared cumulative QC checkpoint even when a caller
# requests a larger batch. This keeps compute bounded while avoiding an
# expensive six-setting re-audit after every individual replication.
CHECKPOINTS=(1 5 10 25 50 100 200 300 400 500)
NEXT_CHECKPOINT=""
for checkpoint in "${CHECKPOINTS[@]}"; do
  if [[ "${checkpoint}" -ge "${START_SIM_ID}" ]]; then
    NEXT_CHECKPOINT="${checkpoint}"
    break
  fi
done
if [[ -z "${NEXT_CHECKPOINT}" ]]; then
  echo "could not resolve the next cumulative QC checkpoint." >&2
  exit 1
fi
CHECKPOINT_END_ROW=$((BATCH_START + NEXT_CHECKPOINT - START_SIM_ID))
if [[ "${BATCH_END}" -gt "${CHECKPOINT_END_ROW}" ]]; then
  BATCH_END="${CHECKPOINT_END_ROW}"
fi
END_SIM_ID="$(read_group_field "${BATCH_END}" 3)"
if [[ ! "${START_SIM_ID}" =~ ^[1-9][0-9]*$ ]] ||
   [[ ! "${END_SIM_ID}" =~ ^[1-9][0-9]*$ ]] ||
   [[ "${END_SIM_ID}" -gt 500 ]]; then
  echo "could not resolve the cumulative replication count." >&2
  exit 1
fi
if [[ ! "${SOURCE_COUNT}" =~ ^(2|4|8)$ ]]; then
  echo "could not resolve K from grouped manifest row ${BATCH_START}." >&2
  exit 1
fi
if [[ ! "${PRIMARY_CUTOFF}" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
   ! awk -v cutoff="${PRIMARY_CUTOFF}" 'BEGIN { exit !(cutoff > 0) }'; then
  echo "could not resolve a positive cutoff from grouped manifest row ${BATCH_START}." >&2
  exit 1
fi
if [[ "${PRIMARY_CUTOFF}" != "${SELECTED_CUTOFF}" ]]; then
  echo "grouped manifest cutoff does not match the audited selected cutoff." >&2
  exit 1
fi

# A new batch cannot outrun review of the most recent predeclared cumulative
# checkpoint. When walking the whole manifest, starting a new config/K block
# requires the previous block's n=500 gate; under ROCE_SETTING the settings
# are independent and only the setting's own ladder applies.
PREVIOUS_REVIEW_N=""
PREVIOUS_CONFIG="${CONFIGURATION}"
PREVIOUS_K="${SOURCE_COUNT}"
if [[ "${START_SIM_ID}" -gt 1 ]]; then
  for checkpoint in "${CHECKPOINTS[@]}"; do
    if [[ "${checkpoint}" -lt "${START_SIM_ID}" ]]; then
      PREVIOUS_REVIEW_N="${checkpoint}"
    fi
  done
elif [[ -z "${SETTING}" && "${BATCH_START}" -gt 1 ]]; then
  PREVIOUS_ROW=$((BATCH_START - 1))
  PREVIOUS_CONFIG="$(read_group_field "${PREVIOUS_ROW}" 4)"
  PREVIOUS_K="$(read_group_field "${PREVIOUS_ROW}" 6)"
  PREVIOUS_REVIEW_N="$(read_group_field "${PREVIOUS_ROW}" 3)"
fi
if [[ -n "${PREVIOUS_REVIEW_N}" ]]; then
  PREVIOUS_GATE="${RESULT_ROOT}/group_checkpoint_gates/$(printf \
    '%s_K%d_n%03d_passed.txt' \
    "${PREVIOUS_CONFIG}" "${PREVIOUS_K}" "${PREVIOUS_REVIEW_N}")"
  if [[ ! -f "${PREVIOUS_GATE}" ]]; then
    echo "prior grouped checkpoint has not passed review: ${PREVIOUS_GATE}" >&2
    exit 1
  fi
  require_gate_value "${PREVIOUS_GATE}" rho_group_checkpoint passed
  require_gate_value "${PREVIOUS_GATE}" config "${PREVIOUS_CONFIG}"
  require_gate_value "${PREVIOUS_GATE}" K "${PREVIOUS_K}"
  require_gate_value \
    "${PREVIOUS_GATE}" expected_replications "${PREVIOUS_REVIEW_N}"
  require_gate_value \
    "${PREVIOUS_GATE}" primary_cutoff "${PRIMARY_CUTOFF}"
  require_gate_value \
    "${PREVIOUS_GATE}" package_fingerprint "${PACKAGE_FINGERPRINT}"
  require_gate_value \
    "${PREVIOUS_GATE}" workflow_fingerprint "${WORKFLOW_FINGERPRINT}"
  require_gate_value \
    "${PREVIOUS_GATE}" manifest_fingerprint "${PRIMARY_MANIFEST_SHA256}"
  require_gate_value \
    "${PREVIOUS_GATE}" checkpoint_audit_driver_md5 \
    "$(md5sum "${PROJECT_ROOT}/scripts/slurm/audit_direct_tate_checkpoint.R" | awk '{print $1}')"
  require_gate_value \
    "${PREVIOUS_GATE}" simulation_qc_policy_md5 \
    "$(md5sum "${PROJECT_ROOT}/scripts/slurm/simulation_qc_policy.R" | awk '{print $1}')"
  require_gate_value \
    "${PREVIOUS_GATE}" group_checkpoint_driver_md5 \
    "$(md5sum "${PROJECT_ROOT}/scripts/slurm/run_rho_group_checkpoint_audit.sh" | awk '{print $1}')"
fi

RUN_CHECKPOINT_AUDIT=0
for checkpoint in "${CHECKPOINTS[@]}"; do
  if [[ "${END_SIM_ID}" -eq "${checkpoint}" ]]; then
    RUN_CHECKPOINT_AUDIT=1
    break
  fi
done

if [[ -n "${ROCE_NUISANCE_CV_THREADS:-}" ]]; then
  roce_resolve_nuisance_cv_threads 1 5
elif [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  ROCE_NUISANCE_CV_THREADS=2
  export ROCE_NUISANCE_CV_THREADS
  roce_resolve_nuisance_cv_threads 2 5
else
  roce_resolve_nuisance_cv_threads 5 5
fi
NUISANCE_CV_THREADS="${ROCE_NUISANCE_CV_THREADS_RESOLVED}"
CPUS_PER_TASK=$((2 * SOURCE_COUNT * NUISANCE_CV_THREADS))
DEFAULT_POSITIVE_RHO_WORKERS=$((CPUS_PER_TASK / NUISANCE_CV_THREADS))
if [[ "${DEFAULT_POSITIVE_RHO_WORKERS}" -gt 5 ]]; then
  DEFAULT_POSITIVE_RHO_WORKERS=5
fi
POSITIVE_RHO_WORKERS="${ROCE_POSITIVE_RHO_WORKERS:-${DEFAULT_POSITIVE_RHO_WORKERS}}"
if [[ ! "${POSITIVE_RHO_WORKERS}" =~ ^[1-9][0-9]*$ ]] ||
   [[ "${POSITIVE_RHO_WORKERS}" -gt "${DEFAULT_POSITIVE_RHO_WORKERS}" ]]; then
  echo "ROCE_POSITIVE_RHO_WORKERS must be an integer from 1 through ${DEFAULT_POSITIVE_RHO_WORKERS}." >&2
  exit 1
fi
if [[ "${SOURCE_COUNT}" -eq 8 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-32G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-24:00:00}"
elif [[ "${SOURCE_COUNT}" -eq 4 ]]; then
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-16G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-18:00:00}"
else
  MEMORY_PER_TASK="${ROCE_MEMORY_PER_TASK:-12G}"
  TIME_PER_TASK="${ROCE_TIME_PER_TASK:-12:00:00}"
fi
ARRAY_SPEC="${BATCH_START}-${BATCH_END}%${MAX_CONCURRENT}"

if [[ "${ROCE_SUBMIT_DRY_RUN:-0}" == "1" ]]; then
  echo "dry run: group rows ${BATCH_START}-${BATCH_END} of ${N_GROUPS}"
  echo "dry run: C=${CONFIGURATION}, K=${SOURCE_COUNT}, array=${ARRAY_SPEC}"
  echo "dry run: primary cutoff=${PRIMARY_CUTOFF}"
  echo "dry run: partition=${SLURM_PARTITION}"
  echo "dry run: ${CPUS_PER_TASK} CPUs, ${NUISANCE_CV_THREADS} CV threads, ${MEMORY_PER_TASK}, ${TIME_PER_TASK}"
  echo "dry run: ${POSITIVE_RHO_WORKERS} concurrent positive-rho PSOCK workers"
  echo "dry run: cumulative n=${END_SIM_ID}, checkpoint audit=${RUN_CHECKPOINT_AUDIT}"
  echo "dry run: tested package ${PROJECT_LIBRARY}"
  exit 0
fi

cd "${PROJECT_ROOT}"
mkdir -p "${OUTPUT_ROOT}" "${SENSITIVITY_OUTPUT_ROOT}" \
  "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_NUISANCE_CV_THREADS="${NUISANCE_CV_THREADS}"
export ROCE_POSITIVE_RHO_WORKERS="${POSITIVE_RHO_WORKERS}"
GROUP_JOB="$(
  ROCE_GROUP_MANIFEST="${GROUP_MANIFEST}" \
  ROCE_PRIMARY_MANIFEST="${PRIMARY_MANIFEST}" \
  ROCE_OUTPUT_ROOT="${OUTPUT_ROOT}" \
  ROCE_SENSITIVITY_OUTPUT_ROOT="${SENSITIVITY_OUTPUT_ROOT}" \
  sbatch --parsable --partition="${SLURM_PARTITION}" \
    --array="${ARRAY_SPEC}" \
    --cpus-per-task="${CPUS_PER_TASK}" \
    --mem="${MEMORY_PER_TASK}" \
    --time="${TIME_PER_TASK}" \
    scripts/slurm/run_direct_tate_rho_group_array.sh
)"
CHECKPOINT_JOB=""
if [[ "${RUN_CHECKPOINT_AUDIT}" -eq 1 ]]; then
  CHECKPOINT_JOB="$(
    ROCE_RESULT_ROOT="${RESULT_ROOT}" \
    ROCE_PRIMARY_MANIFEST="${PRIMARY_MANIFEST}" \
    ROCE_CHECKPOINT_CONFIG="${CONFIGURATION}" \
    ROCE_CHECKPOINT_K="${SOURCE_COUNT}" \
    ROCE_CHECKPOINT_REPLICATIONS="${END_SIM_ID}" \
    ROCE_PRIMARY_CUTOFF="${PRIMARY_CUTOFF}" \
    sbatch --parsable --partition="${SLURM_PARTITION}" \
      --dependency="afterok:${GROUP_JOB}" \
      scripts/slurm/run_rho_group_checkpoint_audit.sh
  )"
fi

echo "submitted rho-group rows ${BATCH_START}-${BATCH_END} (${ARRAY_SPEC}): ${GROUP_JOB}"
if [[ -n "${CHECKPOINT_JOB}" ]]; then
  echo "submitted dependent six-setting checkpoint audit at n=${END_SIM_ID}: ${CHECKPOINT_JOB}"
else
  echo "next six-setting cumulative audit is scheduled at n=${NEXT_CHECKPOINT}"
fi
echo "setting C=${CONFIGURATION}, K=${SOURCE_COUNT}, cutoff=${PRIMARY_CUTOFF}; package=${PROJECT_LIBRARY}"
