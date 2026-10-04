#!/usr/bin/env Rscript
# Exact 1000/site finite-support fixture for the feature-balance counterexample.
# Calls production nuisance primitives at lambda=0. This is not a cross-fitted
# FACE Monte Carlo replicate and does not estimate production coverage.

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: check_package_population.R NEW_OUTPUT_DIRECTORY")
  parent <- normalizePath(dirname(args[[1L]]), mustWork = TRUE)
  if (!startsWith(parent, "/scratch.global/")) stop("Output parent must be in scratch.global")
  if (file.exists(args[[1L]])) stop("Output exists; choose a new directory")
  project_library <- Sys.getenv("ROCE_PROJECT_LIB")
  if (!nzchar(project_library)) stop("Set ROCE_PROJECT_LIB to the frozen package library")
  .libPaths(c(project_library, .libPaths()))
  loadNamespace("RoCE")

  covariate <- c(-1, 0, 1)
  site_mass <- c(.3, .4, .3)
  propensity <- c(.25, .75, .25)
  outcome_mean <- c(.8, .2, .8)
  rows <- list()
  for (index in seq_along(covariate)) for (arm in 0:1) for (response in 0:1) {
    arm_probability <- if (arm == 1L) propensity[index] else 1 - propensity[index]
    response_probability <- if (response == 1L) outcome_mean[index] else 1 - outcome_mean[index]
    count <- 1000 * site_mass[index] * arm_probability * response_probability
    stopifnot(abs(count - round(count)) < 1e-10)
    rows[[length(rows) + 1L]] <- data.frame(
      x = rep(covariate[index], round(count)), a = arm, y = response
    )
  }
  source <- do.call(rbind, rows)
  # Two separate empirical measures with the same exact cell probabilities.
  target <- source
  stopifnot(nrow(source) == 1000L, nrow(target) == 1000L)
  source_design <- matrix(source$x, ncol = 1L)
  target_design <- matrix(target$x, ncol = 1L)
  design <- cbind(1, source_design)
  target_features <- cbind(1, target_design)
  summaries <- list()
  for (arm in 0:1) {
    initial_outcome <- RoCE:::fit_general_glm_cpp(
      X = source_design, Y = source$y, weights = as.numeric(source$a == arm),
      family_int = 1L, link_int = 1L, lambda = 0, max_iter = 1000L,
      tol = 1e-10, warm_start = numeric()
    )
    initial_tilt <- RoCE:::fit_initial_density_ratio_cpp(
      Z_site = source_design, A_source = source$a, mean_phi = colMeans(target_features),
      lambda = 0, max_iter = 1000L, tol = 1e-10, A_val = arm, M_tau = 5,
      warm_start = numeric()
    )
    initial_prediction <- plogis(drop(target_features %*% initial_outcome$alpha))
    derivative_moment <- colMeans(target_features * initial_prediction * (1 - initial_prediction))
    final_tilt <- RoCE:::fit_unified_density_ratio_cpp(
      Z_site = source_design, A_source = source$a, mean_grad_psi = derivative_moment,
      alpha_init = initial_outcome$alpha, lambda = 0, max_iter = 1000L, tol = 1e-10,
      calibrated = TRUE, M_tau = 5, W_outcome = source_design, A_val = arm,
      family_int = 1L, link_int = 1L, warm_start = numeric()
    )
    final_outcome <- RoCE:::fit_unified_outcome_cpp(
      W_outcome = source_design, Y_source = source$y, A_source = source$a,
      gamma_s = initial_tilt$gamma, family_int = 1L, link_int = 1L,
      lambda = 0, max_iter = 1000L, tol = 1e-10, A_val = arm,
      calibrated = TRUE, M_tau = 5, Z_site = source_design, warm_start = numeric()
    )
    stopifnot(initial_outcome$converged, initial_tilt$converged,
              final_outcome$converged, final_tilt$converged)
    prediction <- plogis(drop(design %*% final_outcome$alpha))
    weight <- as.numeric(source$a == arm) * exp(-drop(design %*% final_tilt$gamma))
    balance <- colMeans(design * weight) - colMeans(target_features)
    estimate <- mean(prediction) + mean(weight * (source$y - prediction))
    correction <- RoCE:::calculate_correction_term_cpp(
      Z_site = source_design, A_source = source$a, Y_source = source$y,
      gamma_s = final_tilt$gamma, alpha_ts = final_outcome$alpha,
      W_outcome = source_design, M_tau = 5, family_int = 1L, link_int = 1L, A_val = arm
    )
    stopifnot(is.numeric(correction$delta_ts), length(correction$delta_ts) == 1L,
              is.finite(correction$delta_ts),
              abs(correction$delta_ts - mean(weight * (source$y - prediction))) < 1e-8)
    expected <- sum(site_mass * (if (arm == 1L) propensity else 1 - propensity) * outcome_mean) /
      sum(site_mass * (if (arm == 1L) propensity else 1 - propensity))
    stopifnot(max(abs(balance)) < 1e-8, abs(estimate - expected) < 1e-8)
    summaries[[length(summaries) + 1L]] <- data.frame(
      arm = arm, n_target = nrow(target), n_source = nrow(source),
      estimate = estimate, truth = sum(site_mass * outcome_mean),
      max_balance_error = max(abs(balance)),
      alpha_initial_final_difference = max(abs(initial_outcome$alpha - final_outcome$alpha)),
      gamma_initial_final_difference = max(abs(initial_tilt$gamma - final_tilt$gamma))
    )
  }
  summaries <- do.call(rbind, summaries)
  stopifnot(abs(diff(summaries$estimate) + .29090909090909) < 1e-8)
  dir.create(args[[1L]])
  write.csv(summaries, file.path(args[[1L]], "package_population_checks.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(args[[1L]], "session_info.txt"))
  writeLines(c(paste("library:", find.package("RoCE")),
               "scope: fixed-penalty production primitives; exact finite-support fixture; no Monte Carlo"),
             file.path(args[[1L]], "scope.txt"))
  print(summaries, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
