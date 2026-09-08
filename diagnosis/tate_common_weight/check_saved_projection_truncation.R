#!/usr/bin/env Rscript

# Exercise derivative branches using copies of saved coefficients. No refit,
# resampling, new data or claim that perturbed coefficients solve any fit.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: check_saved_projection_truncation.R V19_LIBRARY MATRIX_PROBE")
  lib <- normalizePath(args[1], mustWork = TRUE)
  directory <- normalizePath(args[2], mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"),
    .weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")
  lines <- readLines(file.path(directory, "sha256.txt"))
  expected <- substr(lines[substring(lines, 67L) == "matrix_probe.rds"], 1L, 64L)
  stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, "matrix_probe.rds")) == expected)
  source_env <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", source_env)
  evaluate <- source_env$.evaluate_saved_projection
  state <- readRDS(file.path(directory, "matrix_probe.rds"))$state
  cases <- c(initial_outcome = "alpha_initial_k2_2",
             initial_tilting = "gamma_initial_k2_2", final_tilting = "gamma_final")
  checks <- lapply(names(cases), function(label) {
    theta <- state$theta
    theta[state$indices[[cases[[label]]]][1L]] <- theta[state$indices[[cases[[label]]]][1L]]+6
    analytical <- evaluate(state, theta, TRUE)
    stopifnot(analytical$truncation_counts[[label]] > 0L)
    errors <- vapply(1:8, function(k) {
      direction <- cos(seq_along(theta)*k); direction <- direction/sqrt(sum(direction^2))
      epsilon <- 1e-6
      plus <- evaluate(state, theta+epsilon*direction)
      minus <- evaluate(state, theta-epsilon*direction)
      c(max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(analytical$jacobian %*% direction))),
        abs((plus$score-minus$score)/(2*epsilon)-sum(analytical$score_gradient*direction)))
    }, numeric(2L))
    native <- RoCE:::calculate_correction_term_cpp(state$source$Z[, -1, drop = FALSE],
      state$source$A, state$source$Y, theta[state$indices$gamma_final],
      theta[state$indices$alpha_final], state$source$W[, -1, drop = FALSE], 5, 1L, 1L, 1L)
    native_score <- mean(RoCE:::predict_glm_cpp(state$target$W[, -1, drop = FALSE],
      theta[state$indices$alpha_final], 1L, 1L))+native$delta_ts
    stopifnot(max(errors) < 1e-7, abs(native_score-analytical$score) < 1e-12)
    data.frame(branch = label, activated_records = analytical$truncation_counts[[label]],
      jacobian_error = max(errors[1, ]), score_gradient_error = max(errors[2, ]),
      native_score_error = abs(native_score-analytical$score),
      minimum_boundary_distance = analytical$boundary_distance)
  })
  print(do.call(rbind, checks), row.names = FALSE, digits = 12)
  cat("Coefficient-perturbation derivative checks only; fitted states and all outputs remain unchanged.\n")
}

if (sys.nframe() == 0L) main()
