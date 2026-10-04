#!/usr/bin/env Rscript
# Population audit of the treated-arm score at a source-specific biased mean.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: frozen_R_library new_scratch_output")
library_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
output <- arguments[[2L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/") || dir.exists(output) ||
    !dir.create(output, recursive = TRUE)) stop("Use a new FACE-HD scratch directory")

evaluate <- function(order) {
  quadrature <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1, quadrature$X)
  mass <- quadrature$weights
  base <- drop(quadrature$X %*% c(.25, .20, .10, .05))
  target_mean <- plogis(.5 + base)
  rows <- list()
  for (scenario in c("C1", "C3")) {
    strength <- RoCE:::.face_misspecification_strengths(scenario, .75)$propensity
    assignment <- RoCE:::.bounded_active_features(quadrature$X, strength)
    for (site in 1:2) {
      common <- site / 2 * c(.25, -.20, .15, -.10)
      contrast <- c(.35, -.25, .125, -.0625) / 2
      exponential <- exp(assignment %*% cbind(common - contrast, common + contrast))
      joint <- exponential[, 2L] / sum(mass * rowSums(exponential))
      for (rho in c(0, .5, 1, 2)) {
        source_mean <- plogis(.5 + rho + base)
        source_derivative <- source_mean * (1 - source_mean)
        for (initialization in c("target", "source")) {
          initial_mean <- if (initialization == "target") target_mean else source_mean
          initial_derivative <- initial_mean * (1 - initial_mean)
          fit <- glm.fit(design, 1 / joint, weights = mass * joint * initial_derivative,
            family = quasipoisson(), control = glm.control(epsilon = 1e-12, maxit = 100))
          calibration_gradient <- drop(crossprod(design,
            mass * initial_derivative * (joint * fit$fitted.values - 1)))
          score_gradient <- drop(crossprod(design,
            mass * source_derivative * (1 - joint * fit$fitted.values)))
          predictor_bound <- abs(fit$coefficients[1L]) +
            RoCE:::.bounded_covariate_spec()$bound * sum(abs(fit$coefficients[-1L]))
          stopifnot(isTRUE(fit$converged), max(abs(calibration_gradient)) < 1e-9, predictor_bound < 12)
          source_alpha <- c(.5 + rho, .25, .20, .10, .05)
          score <- function(alpha) {
            predicted <- plogis(drop(design %*% alpha))
            sum(mass * (predicted + joint * fit$fitted.values * (source_mean - predicted)))
          }
          finite_difference <- vapply(seq_along(source_alpha), function(index) {
            step <- numeric(length(source_alpha)); step[index] <- 1e-5
            (score(source_alpha + step) - score(source_alpha - step)) / 2e-5
          }, numeric(1L))
          derivative_error <- max(abs(score_gradient - finite_difference))
          stopifnot(derivative_error < 1e-8)
          if (initialization == "source" || scenario == "C1" || rho == 0) {
            stopifnot(max(abs(score_gradient)) < 1e-9)
          }
          rows[[length(rows) + 1L]] <- data.frame(config = scenario, site = site, rho = rho,
            initialization = initialization, source_target_mean = sum(mass * source_mean),
            target_mean = sum(mass * target_mean),
            candidate_limit_bias = score(source_alpha) - sum(mass * source_mean),
            calibration_gradient = max(abs(calibration_gradient)),
            outcome_score_gradient = max(abs(score_gradient)), derivative_check_error = derivative_error,
            weight_predictor_bound = predictor_bound)
        }
      }
    }
  }
  do.call(rbind, rows)
}
coarse <- evaluate(16L)
fine <- evaluate(24L)
fields <- c("source_target_mean", "candidate_limit_bias", "calibration_gradient", "outcome_score_gradient")
quadrature_error <- max(abs(as.matrix(coarse[fields]) - as.matrix(fine[fields])))
stopifnot(quadrature_error < 1e-8)
write.csv(fine, file.path(output, "population_gradients.csv"), row.names = FALSE)
writeLines(c(paste("Maximum16/24 quadrature error:", quadrature_error),
  "Audit of regularity around the SOURCE-specific target-population mean, not source validity for the target effect.",
  "No finite-sample estimator or default DGP was changed."), file.path(output, "checks.txt"))
cat("BIASED_SOURCE_POPULATION_AUDIT_PASSED\n")
