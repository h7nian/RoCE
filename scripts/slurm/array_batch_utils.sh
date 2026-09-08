#!/bin/bash

# Shared validation for bounded SLURM array submissions. Callers set
# ROCE_ARRAY_START, ROCE_BATCH_SIZE, and ROCE_MAX_CONCURRENT as needed.
resolve_roce_array_start_from_outputs() {
  local manifest_path="$1"
  local output_root="$2"

  ROCE_ALL_OUTPUTS_COMPLETE=0
  if [[ -n "${ROCE_ARRAY_START:-}" ]]; then
    return 0
  fi

  # The Slurm array index is a manifest row, whereas split manifests retain a
  # global task_id for filenames. Keep those two identifiers distinct.
  local next_missing_row=""
  local manifest_row
  local manifest_task_id
  local task_output
  while IFS=, read -r manifest_row manifest_task_id; do
    manifest_task_id="${manifest_task_id//\"/}"
    if [[ ! "${manifest_task_id}" =~ ^[1-9][0-9]*$ ]]; then
      echo "invalid task_id in manifest row ${manifest_row}: ${manifest_task_id}" >&2
      return 1
    fi
    task_output="${output_root}/$(printf 'task_%06d.csv' "${manifest_task_id}")"
    if [[ ! -s "${task_output}" ]]; then
      next_missing_row="${manifest_row}"
      break
    fi
  done < <(
    awk -F, 'NR > 1 { task_id = $1; gsub(/"/, "", task_id); print NR - 1 "," task_id }' \
      "${manifest_path}"
  )

  if [[ -z "${next_missing_row}" ]]; then
    ROCE_ALL_OUTPUTS_COMPLETE=1
    return 0
  fi
  export ROCE_ARRAY_START="${next_missing_row}"
}

configure_roce_array_batch() {
  local total_tasks="$1"
  local default_batch_size="${2:-100}"
  local default_max_concurrent="${3:-20}"

  ROCE_BATCH_START="${ROCE_ARRAY_START:-1}"
  ROCE_BATCH_SIZE_RESOLVED="${ROCE_BATCH_SIZE:-${default_batch_size}}"
  ROCE_MAX_CONCURRENT_RESOLVED="${ROCE_MAX_CONCURRENT:-${default_max_concurrent}}"

  for value_name in \
    ROCE_BATCH_START \
    ROCE_BATCH_SIZE_RESOLVED \
    ROCE_MAX_CONCURRENT_RESOLVED; do
    local value="${!value_name}"
    if [[ ! "${value}" =~ ^[1-9][0-9]*$ ]]; then
      echo "${value_name} must be a positive integer: ${value}" >&2
      return 1
    fi
  done
  if [[ "${ROCE_BATCH_START}" -gt "${total_tasks}" ]]; then
    echo "array start ${ROCE_BATCH_START} exceeds ${total_tasks} tasks." >&2
    return 1
  fi

  local candidate_end=$((ROCE_BATCH_START + ROCE_BATCH_SIZE_RESOLVED - 1))
  if [[ "${candidate_end}" -lt "${total_tasks}" ]]; then
    ROCE_BATCH_END="${candidate_end}"
  else
    ROCE_BATCH_END="${total_tasks}"
  fi
  ROCE_ARRAY_SPEC="${ROCE_BATCH_START}-${ROCE_BATCH_END}%${ROCE_MAX_CONCURRENT_RESOLVED}"
}

# MSI guardrail for expensive p=100 jobs. Keep this separate from the default
# values above: defaults choose the usual submission size, whereas this helper
# rejects an accidental broad override before sbatch is called.
enforce_roce_array_safety_cap() {
  local maximum_batch_size="${1:-25}"
  local maximum_concurrent="${2:-5}"

  if [[ "${ROCE_BATCH_SIZE_RESOLVED}" -gt "${maximum_batch_size}" ]]; then
    echo "refusing batch size ${ROCE_BATCH_SIZE_RESOLVED}; hard cap is ${maximum_batch_size}" >&2
    return 1
  fi
  if [[ "${ROCE_MAX_CONCURRENT_RESOLVED}" -gt "${maximum_concurrent}" ]]; then
    echo "refusing concurrency ${ROCE_MAX_CONCURRENT_RESOLVED}; hard cap is ${maximum_concurrent}" >&2
    return 1
  fi
}

# Prevent a bounded array from spilling into the next statistical setting.
# Manifest rows are grouped by experiment/configuration/p/K/rho/cutoff, while
# sim_id varies within a setting.  The first eight columns are common to the
# main and cutoff-diagnostic manifests; the final field is named either cutoff
# or cutoffs but has the same role in the grouping key.
cap_roce_array_batch_to_manifest_setting() {
  local manifest_path="$1"
  local proposed_end="${ROCE_BATCH_END}"
  local capped_end

  capped_end="$(awk -F, \
    -v start="${ROCE_BATCH_START}" \
    -v proposed_end="${proposed_end}" '
      function clean(value) {
        gsub(/^"|"$/, "", value)
        return value
      }
      NR == start + 1 {
        setting_key = clean($2) "|" clean($4) "|" clean($5) "|" \
          clean($6) "|" clean($7) "|" clean($8)
        next
      }
      NR > start + 1 && NR <= proposed_end + 1 {
        candidate_key = clean($2) "|" clean($4) "|" clean($5) "|" \
          clean($6) "|" clean($7) "|" clean($8)
        if (candidate_key != setting_key) {
          print NR - 2
          found_boundary = 1
          exit
        }
      }
      END {
        if (!found_boundary) print proposed_end
      }
    ' "${manifest_path}")"

  if [[ ! "${capped_end}" =~ ^[1-9][0-9]*$ ]] ||
      [[ "${capped_end}" -lt "${ROCE_BATCH_START}" ]] ||
      [[ "${capped_end}" -gt "${proposed_end}" ]]; then
    echo "failed to resolve a valid within-setting array end: ${capped_end}" >&2
    return 1
  fi

  ROCE_BATCH_END="${capped_end}"
  ROCE_ARRAY_SPEC="${ROCE_BATCH_START}-${ROCE_BATCH_END}%${ROCE_MAX_CONCURRENT_RESOLVED}"
}
