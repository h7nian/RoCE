#!/usr/bin/env Rscript

# Constructive population audit of an estimating-equation correction. This is
# not a high-dimensional implementation, a new simulation, or a CI adjustment.
# Original outcome and tilting bases remain different in C2/C3-style cases.
main <- function(args = commandArgs(trailingOnly = TRUE), arm_value = 1L,
                 target_treatment_probability = rep(.5, 3L),
                 return_details = FALSE, verbose = TRUE,
                 projection_solver = NULL) {
  if (length(args) != 1L) stop("usage: check_full_nuisance_orthogonalization.R V19_LIBRARY")
  lib <- normalizePath(args[[1L]], mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"))
  source("scripts/slurm/result_provenance.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  stopifnot(.weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")
  support <- c(-1, 0, 1)
  target_mass <- rep(1/3, 3)
  stopifnot(length(arm_value) == 1L, arm_value %in% c(0L, 1L),
            length(target_treatment_probability) == 3L,
            all(is.finite(target_treatment_probability)),
            all(target_treatment_probability > 0 & target_treatment_probability < 1))
  source_arm_mass <- if (arm_value == 1L) c(.08, .28, .12) else
    c(.3, .4, .3)-c(.08, .28, .12)
  target_arm_probability <- if (arm_value == 1L) target_treatment_probability else
    1-target_treatment_probability
  true_mean <- plogis(-.2 + .7*support + 1.2*support^2 - .8*(1-arm_value))
  source_x <- rep(support, c(90L, 120L, 90L))
  source_a <- c(rep(c(1, 0), c(24L, 66L)), rep(c(1, 0), c(84L, 36L)),
                rep(c(1, 0), c(36L, 54L)))
  source_mean <- rep(true_mean, c(90L, 120L, 90L))
  basis <- function(x, quadratic) {
    if (quadratic) cbind(x, x^2) else matrix(x, ncol = 1L)
  }
  numeric_jacobian <- function(fun, theta, step = 1e-5) {
    vapply(seq_along(theta), function(j) {
      shift <- numeric(length(theta)); shift[j] <- step
      (fun(theta + shift)-fun(theta-shift))/(2*step)
    }, numeric(length(fun(theta))))
  }
  cases <- list(C1_style = c(TRUE, TRUE), C2_style = c(FALSE, TRUE),
                C3_style = c(TRUE, FALSE))
  case_details <- list()
  result <- lapply(names(cases), function(label) {
    features <- cases[[label]]
    W0 <- basis(support, features[1]); Z0 <- basis(support, features[2])
    W <- cbind(1, W0); Z <- cbind(1, Z0)
    source_W <- basis(source_x, features[1]); source_Z <- basis(source_x, features[2])
    a <- ncol(W); g <- ncol(Z)
    indices <- list(alpha_initial = seq_len(a), gamma_initial = a+seq_len(g),
                    alpha_final = a+g+seq_len(a), gamma_final = 2*a+g+seq_len(g))
    # Fractional Y values here evaluate the exact expected Bernoulli loss;
    # they are not simulated responses or a replacement for a binary dataset.
    alpha_initial <- RoCE:::fit_general_glm_cpp(W0, true_mean,
      3*target_mass*target_arm_probability, 1L, 1L, 0, 10000L, 1e-12, numeric())
    gamma_initial <- RoCE:::fit_initial_density_ratio_cpp(source_Z, source_a,
      colMeans(Z), 0, 10000L, 1e-12, arm_value, numeric())
    alpha_final <- RoCE:::fit_unified_outcome_cpp(source_W, source_mean, source_a,
      gamma_initial$gamma, 1L, 1L, 0, 10000L, 1e-12, arm_value, TRUE, 5,
      source_Z, numeric())
    m_initial <- plogis(drop(W %*% alpha_initial$alpha))
    gamma_final <- RoCE:::fit_unified_density_ratio_cpp(source_Z, source_a,
      drop(crossprod(Z, target_mass*m_initial*(1-m_initial))),
      alpha_initial$alpha, 0, 10000L, 1e-12, TRUE, 5,
      source_W, arm_value, 1L, 1L, numeric())
    fits <- list(alpha_initial, gamma_initial, alpha_final, gamma_final)
    stopifnot(all(vapply(fits, function(x) isTRUE(x$converged), logical(1))))
    theta <- c(alpha_initial$alpha, gamma_initial$gamma,
               alpha_final$alpha, gamma_final$gamma)
    evaluate <- function(theta) {
      m0 <- plogis(drop(W %*% theta[indices$alpha_initial]))
      m <- plogis(drop(W %*% theta[indices$alpha_final]))
      w0 <- exp(-drop(Z %*% theta[indices$gamma_initial]))
      w <- exp(-drop(Z %*% theta[indices$gamma_final]))
      moments <- c(
        crossprod(W, target_arm_probability*target_mass*(true_mean-m0)),
        crossprod(Z, target_mass-source_arm_mass*w0),
        crossprod(W, source_arm_mass*w0*(true_mean-m)),
        crossprod(Z, (target_mass-source_arm_mass*w)*m0*(1-m0)))
      list(score = sum(target_mass*m + source_arm_mass*w*(true_mean-m)),
           moments = as.numeric(moments), m0 = m0, m = m, w0 = w0, w = w)
    }
    population <- evaluate(theta)
    stopifnot(max(abs(drop(W %*% alpha_initial$alpha))) < 5,
              max(abs(drop(W %*% alpha_final$alpha))) < 5,
              max(abs(drop(Z %*% gamma_final$gamma))) < 5,
              max(abs(drop(cbind(1, source_Z) %*% gamma_initial$gamma))) < 5,
              max(abs(population$moments)) < 1e-8)
    d0 <- population$m0*(1-population$m0)
    d <- population$m*(1-population$m)
    moment_jacobian <- matrix(0, length(theta), length(theta))
    block <- function(left, diagonal, right) crossprod(left, diagonal*right)
    moment_jacobian[indices$alpha_initial, indices$alpha_initial] <-
      -block(W, target_arm_probability*target_mass*d0, W)
    moment_jacobian[indices$gamma_initial, indices$gamma_initial] <- block(Z, source_arm_mass*population$w0, Z)
    moment_jacobian[indices$alpha_final, indices$alpha_final] <- -block(W, source_arm_mass*population$w0*d, W)
    moment_jacobian[indices$alpha_final, indices$gamma_initial] <-
      -block(W, source_arm_mass*population$w0*(true_mean-population$m), Z)
    moment_jacobian[indices$gamma_final, indices$gamma_final] <- block(Z, source_arm_mass*population$w*d0, Z)
    moment_jacobian[indices$gamma_final, indices$alpha_initial] <-
      block(Z, (target_mass-source_arm_mass*population$w)*d0*(1-2*population$m0), W)
    score_gradient <- numeric(length(theta))
    score_gradient[indices$alpha_final] <- crossprod(W, (target_mass-source_arm_mass*population$w)*d)
    score_gradient[indices$gamma_final] <- -crossprod(Z, source_arm_mass*population$w*(true_mean-population$m))
    jacobian_error <- max(abs(moment_jacobian - numeric_jacobian(function(t) evaluate(t)$moments, theta)))
    gradient_error <- max(abs(score_gradient - as.numeric(numeric_jacobian(function(t) evaluate(t)$score, theta))))
    stopifnot(rcond(moment_jacobian) > 1e-5, jacobian_error < 1e-8, gradient_error < 1e-8)
    # All four nuisance blocks, including initial fits, enter this system.
    # No coefficient is chosen by coverage or by the estimator's error.
    projection <- if (is.null(projection_solver)) solve(t(moment_jacobian), score_gradient) else
      projection_solver(moment_jacobian, score_gradient, indices)
    stopifnot(is.numeric(projection), length(projection) == length(theta),
              all(is.finite(projection)))
    corrected_score <- function(t) {
      state <- evaluate(t)
      state$score - sum(projection*state$moments)
    }
    corrected_gradient <- as.numeric(numeric_jacobian(corrected_score, theta))
    target_truth <- sum(target_mass*true_mean)
    stopifnot(max(abs(corrected_gradient)) < 1e-8,
              abs(population$score-target_truth) < 1e-8,
              abs(corrected_score(theta)-target_truth) < 1e-8)
    if (label != "C1_style") stopifnot(max(abs(score_gradient)) > 1e-3)
    if (return_details) case_details[[label]] <<- list(
      arm_value = arm_value, W = W, Z = Z, theta = theta, indices = indices,
      projection = projection, population_moments = population$moments,
      target_mean = target_truth, corrected_population_score = corrected_score(theta))
    data.frame(case = label, outcome_dimension = a, tilting_dimension = g,
      full_nuisance_dimension = length(theta),
      native_moment_residual = max(abs(population$moments)),
      analytic_jacobian_error = jacobian_error,
      uncorrected_gradient_max = max(abs(score_gradient)),
      corrected_gradient_max = max(abs(corrected_gradient)),
      projection_equation_error = max(abs(drop(t(moment_jacobian)%*%projection)-score_gradient)),
      reciprocal_condition = rcond(moment_jacobian),
      population_target_identity_error = abs(corrected_score(theta)-target_truth))
  })
  result <- do.call(rbind, result)
  if (return_details) attr(result, "case_details") <- case_details
  if (verbose) {
    print(result, row.names = FALSE, digits = 12)
    cat("No bootstrap, simulated replications, SE inflation or feature-union substitution.\n",
      "Finite-dimensional population construction only: high-dimensional projection rates,\n",
      "truncation, sample implementation, federated summaries, and adaptive aggregation\n",
      "remain unvalidated. This does not resolve the C1 weight-stability problem.\n", sep = "")
  }
  invisible(result)
}

if (sys.nframe() == 0L) main()
