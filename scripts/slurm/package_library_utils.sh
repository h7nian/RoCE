#!/bin/bash

roce_require_package_library() {
  local project_library="$1"
  local package_root="${project_library}/RoCE"
  local required_file
  for required_file in \
    "${package_root}/DESCRIPTION" \
    "${package_root}/libs/RoCE.so" \
    "${package_root}/R/RoCE.rdb" \
    "${package_root}/R/RoCE.rdx"; do
    if [[ ! -f "${required_file}" ]]; then
      echo "tested RoCE installation is incomplete: ${required_file}" >&2
      return 1
    fi
  done
}

roce_package_fingerprint() {
  local project_library="$1"
  local package_root="${project_library}/RoCE"
  roce_require_package_library "${project_library}"
  {
    sha256sum "${package_root}/DESCRIPTION"
    sha256sum "${package_root}/libs/RoCE.so"
    sha256sum "${package_root}/R/RoCE.rdb"
    sha256sum "${package_root}/R/RoCE.rdx"
  } | awk '{print $1}' | sha256sum | awk '{print $1}'
}

roce_resolve_package_library() {
  local project_library="$1"
  roce_require_package_library "${project_library}"
  readlink -f "${project_library}"
}

roce_files_fingerprint() {
  if [[ "$#" -lt 1 ]]; then
    echo "roce_files_fingerprint requires at least one file" >&2
    return 1
  fi
  local fingerprint_file
  for fingerprint_file in "$@"; do
    if [[ ! -f "${fingerprint_file}" ]]; then
      echo "fingerprint input is missing: ${fingerprint_file}" >&2
      return 1
    fi
  done
  sha256sum "$@" | awk '{print $1}' | sha256sum | awk '{print $1}'
}

roce_package_source_fingerprint() {
  local project_root="$1"
  if [[ ! -d "${project_root}/R" || ! -d "${project_root}/src" ||
        ! -d "${project_root}/tests" ]]; then
    echo "package source tree is incomplete: ${project_root}" >&2
    return 1
  fi
  (
    cd "${project_root}"
    find DESCRIPTION LICENSE NAMESPACE README.md R src inst man tests \
      -type f ! -name '*.o' ! -name '*.so' ! -name '*.dll' -print0 \
      | sort -z | xargs -0 sha256sum
  ) | sha256sum | awk '{print $1}'
}

roce_simulation_workflow_fingerprint() {
  local project_root="$1"
  roce_files_fingerprint \
    "${project_root}/scripts/slurm/run_direct_tate_array.sh" \
    "${project_root}/scripts/slurm/run_direct_tate_task.R" \
    "${project_root}/scripts/slurm/run_direct_tate_rho_group_array.sh" \
    "${project_root}/scripts/slurm/run_direct_tate_rho_group_task.R" \
    "${project_root}/scripts/slurm/direct_tate_task_helpers.R" \
    "${project_root}/scripts/slurm/direct_tate_sensitivity_rows.R" \
    "${project_root}/scripts/slurm/result_provenance.R" \
    "${project_root}/scripts/slurm/resource_topology.R" \
    "${project_root}/scripts/slurm/resource_topology.sh" \
    "${project_root}/scripts/slurm/package_library_utils.sh"
}

roce_cutoff_workflow_fingerprint() {
  local project_root="$1"
  roce_files_fingerprint \
    "${project_root}/scripts/slurm/run_direct_tate_cutoff_array.sh" \
    "${project_root}/scripts/slurm/run_direct_tate_cutoff_task.R" \
    "${project_root}/scripts/slurm/direct_tate_task_helpers.R" \
    "${project_root}/scripts/slurm/result_provenance.R" \
    "${project_root}/scripts/slurm/resource_topology.R" \
    "${project_root}/scripts/slurm/resource_topology.sh" \
    "${project_root}/scripts/slurm/package_library_utils.sh"
}
