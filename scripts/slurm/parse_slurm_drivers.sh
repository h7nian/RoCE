#!/bin/bash
#SBATCH --job-name=roce_parse
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=results/direct_tate_mc500_b5000/logs/%j_parse.out
#SBATCH --error=results/direct_tate_mc500_b5000/logs/%j_parse.err

set -euo pipefail

PROJECT_ROOT="${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
PROJECT_LIBRARY="${ROCE_PROJECT_LIB:-${PROJECT_ROOT}/results/direct_tate_mc500_b5000/Rlib_current}"

source "${PROJECT_ROOT}/scripts/slurm/package_library_utils.sh"
PROJECT_LIBRARY="$(roce_resolve_package_library "${PROJECT_LIBRARY}")"
PACKAGE_FINGERPRINT="$(roce_package_fingerprint "${PROJECT_LIBRARY}")"
WORKFLOW_FINGERPRINT="$(roce_simulation_workflow_fingerprint "${PROJECT_ROOT}")"

module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER="/users/0/zhan9381/Rlibs"
export ROCE_PROJECT_LIB="${PROJECT_LIBRARY}"
export ROCE_PACKAGE_FINGERPRINT="${PACKAGE_FINGERPRINT}"
export ROCE_WORKFLOW_FINGERPRINT="${WORKFLOW_FINGERPRINT}"

cd "${PROJECT_ROOT}"
mkdir -p "${PROJECT_ROOT}/results/direct_tate_mc500_b5000/logs"
Rscript -e '
  scripts <- list.files(
    "scripts/slurm", pattern = "[.]R$", full.names = TRUE
  )
  invisible(lapply(scripts, parse))
  cat(sprintf("parsed %d Slurm R drivers\n", length(scripts)))

  source("scripts/slurm/result_provenance.R")
  sha_fixture <- tempfile("roce_sha256_")
  writeLines("RoCE", sha_fixture)
  sha_fixture_digest <- roce_sha256_file(sha_fixture)
  stopifnot(grepl("^[0-9a-f]{64}$", sha_fixture_digest))
  unlink(sha_fixture)
  tested_fingerprint <- Sys.getenv("ROCE_PACKAGE_FINGERPRINT")
  stopifnot(identical(
    roce_sha256_environment("ROCE_WORKFLOW_FINGERPRINT"),
    Sys.getenv("ROCE_WORKFLOW_FINGERPRINT")
  ))
  valid <- data.frame(
    package_library = rep("/tmp/roce-tested-library", 2L),
    package_fingerprint = rep(tested_fingerprint, 2L)
  )
  stopifnot(roce_result_provenance(
    valid, "/tmp/roce-tested-library"
  )$passed)
  stopifnot(!roce_result_provenance(
    valid[c("package_library")], "/tmp/roce-tested-library"
  )$passed)
  mixed <- valid
  mixed$package_fingerprint[[2L]] <- paste(rep("b", 64L), collapse = "")
  stopifnot(!roce_result_provenance(
    mixed, "/tmp/roce-tested-library"
  )$passed)
  .libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  runtime <- roce_runtime_package_provenance(
    Sys.getenv("ROCE_PROJECT_LIB")
  )
  stopifnot(
    identical(runtime$library, Sys.getenv("ROCE_PROJECT_LIB")),
    identical(runtime$fingerprint, Sys.getenv("ROCE_PACKAGE_FINGERPRINT"))
  )

  source("scripts/slurm/resource_topology.R")
  full_plan <- roce_slurm_resource_plan(
    source_count = 4L, n_folds = 5L,
    allocated_cores = 40L, nuisance_cv_threads = 5L
  )
  fallback_plan <- roce_slurm_resource_plan(
    source_count = 4L, n_folds = 5L,
    allocated_cores = 20L, nuisance_cv_threads = 5L
  )
  stopifnot(
    identical(full_plan$source_workers, 4L),
    isTRUE(full_plan$parallel_treatment_arms),
    identical(full_plan$fully_parallel_cores, 40L),
    identical(fallback_plan$source_workers, 4L),
    !fallback_plan$parallel_treatment_arms
  )
  expect_resource_error <- function(code) {
    inherits(try(code(), silent = TRUE), "try-error")
  }
  stopifnot(
    expect_resource_error(function() {
      roce_slurm_resource_plan(4L, 5L, 4L, 5L)
    }),
    expect_resource_error(function() {
      roce_slurm_resource_plan(4L, 5L, 40L, 6L)
    })
  )
  resource_rows <- data.frame(
    K = rep(4L, 2L), n_folds = rep(5L, 2L),
    allocated_cores = rep(40L, 2L),
    nuisance_cv_threads = rep(5L, 2L),
    source_workers_per_arm = rep(4L, 2L),
    parallel_treatment_arms = rep(TRUE, 2L),
    fully_parallel_cores = rep(40L, 2L)
  )
  stopifnot(roce_validate_result_resource_metadata(resource_rows, 5L))
  resource_rows$source_workers_per_arm[[1L]] <- 3L
  stopifnot(expect_resource_error(function() {
    roce_validate_result_resource_metadata(resource_rows, 5L)
  }))
  # Restore the valid fixture before using it to construct the production
  # topology below. The preceding mutation exists only to exercise fail-fast
  # validation and must not leak into the next independent check.
  resource_rows$source_workers_per_arm[[1L]] <- 4L
  production_rows <- rbind(
    transform(resource_rows[1L, ], K = 4L),
    data.frame(
      K = 8L, n_folds = 5L, allocated_cores = 32L,
      nuisance_cv_threads = 2L, source_workers_per_arm = 8L,
      parallel_treatment_arms = TRUE, fully_parallel_cores = 32L
    )
  )
  stopifnot(roce_validate_production_resource_metadata(production_rows))
  production_rows$nuisance_cv_threads[[2L]] <- 5L
  stopifnot(expect_resource_error(function() {
    roce_validate_production_resource_metadata(production_rows)
  }))
  cat("resource-topology checks passed\n")

  source("scripts/slurm/direct_tate_sensitivity_rows.R")
  if (exists(
      ".face_heterogeneity_type", envir = asNamespace("RoCE"),
      inherits = FALSE
  )) {
    mock_task <- data.frame(
      sim_id = 1L, K = 4L, p = 100L, config = "C3", rho = 0,
      task_id = 1L, experiment = "single_task_smoke", n_site = 1000L,
      n_folds = 5L
    )
    mock_fit <- list(estimate = 0.2, se = 0.03, M_tau = 5)
    mock_row <- roce_make_tate_result_row(
      task = mock_task, fit = mock_fit, method = "target_only_ate",
      truth = 0.21, cutoff = 2, M_tau_inference = 5,
      primary_cutoff = 2, n_total = 5000L, nlambda_init = 100L,
      n_bootstrap = 5000L, outcome_family = "binomial"
    )
    stopifnot(
      identical(mock_row$primary_experiment, "single_task_smoke"),
      identical(mock_row$experiment, "c3_reused_sensitivity"),
      identical(mock_row$primary_cutoff, 2),
      abs(mock_row$bias - (mock_row$estimate - mock_row$truth)) <= 1e-15
    )
  } else {
    message(
      "legacy tested package lacks .face_heterogeneity_type; ",
      "sidecar behavior is deferred to the fresh-package test gate"
    )
  }
  cat("result and sidecar provenance checks passed\n")
'
