#!/bin/bash

roce_resolve_nuisance_cv_threads() {
  local default_threads="${1:-1}"
  local maximum_threads="${2:-5}"
  local requested_threads="${ROCE_NUISANCE_CV_THREADS:-${default_threads}}"

  if [[ ! "${default_threads}" =~ ^[1-9][0-9]*$ ]] ||
      [[ ! "${maximum_threads}" =~ ^[1-9][0-9]*$ ]] ||
      [[ "${default_threads}" -gt "${maximum_threads}" ]]; then
    echo "invalid internal nuisance-CV thread bounds" >&2
    return 1
  fi
  if [[ ! "${requested_threads}" =~ ^[1-9][0-9]*$ ]] ||
      [[ "${requested_threads}" -gt "${maximum_threads}" ]]; then
    echo "ROCE_NUISANCE_CV_THREADS must be an integer from 1 through ${maximum_threads}." >&2
    return 1
  fi

  ROCE_NUISANCE_CV_THREADS_RESOLVED="${requested_threads}"
}

roce_require_nuisance_cv_cpu_capacity() {
  local allocated_cores="$1"
  local nuisance_cv_threads="$2"

  if [[ ! "${allocated_cores}" =~ ^[1-9][0-9]*$ ]]; then
    echo "allocated CPU count must be a positive integer: ${allocated_cores}" >&2
    return 1
  fi
  if [[ ! "${nuisance_cv_threads}" =~ ^[1-9][0-9]*$ ]]; then
    echo "nuisance-CV thread count must be a positive integer: ${nuisance_cv_threads}" >&2
    return 1
  fi
  if [[ "${allocated_cores}" -lt "${nuisance_cv_threads}" ]]; then
    echo "allocated CPUs must be at least the nuisance-CV thread count." >&2
    return 1
  fi
}
