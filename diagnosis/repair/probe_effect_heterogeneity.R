#!/usr/bin/env Rscript
# Prospective population sensitivity analysis; no simulation defaults change.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: frozen_library scratch_output")
.libPaths(c(arguments[[1L]], .libPaths()))
suppressPackageStartupMessages(library(RoCE))
output <- arguments[[2L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
contrast_direction <- c(.8, -.3, .15, -.08)

evaluate_design <- function(order) {
  grid <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1, grid$X)
  weights <- grid$weights
  bound <- RoCE:::.bounded_covariate_spec()$bound
  rows <- list()
  fit_projection <- function(response, mass, family) {
    # Quasi families solve the canonical population score with fractional
    # responses; they do not change the binary finite-sample outcome model.
    fit <- glm.fit(design, response, weights = mass, family = family,
      control = glm.control(epsilon = 1e-12, maxit = 100))
    gradient <- max(abs(crossprod(design, mass * (fit$fitted.values - response))))
    stopifnot(isTRUE(fit$converged), all(is.finite(fit$coefficients)), gradient < 1e-9)
    list(prediction = fit$fitted.values, gradient = gradient,
      predictor_bound = abs(fit$coefficients[1]) + bound * sum(abs(fit$coefficients[-1])))
  }
  for (scenario in c("C1", "C2", "C3")) {
    strengths <- RoCE:::.face_misspecification_strengths(scenario, .75)
    outcome_basis <- RoCE:::.bounded_active_features(grid$X, strengths$outcome)
    assignment_basis <- RoCE:::.bounded_active_features(grid$X, strengths$propensity)
    base_predictor <- drop(outcome_basis %*% c(.25, .20, .10, .05))
    for (heterogeneity in c(0, .5, 1, 2)) {
      contrast <- heterogeneity * drop(outcome_basis %*% contrast_direction)
      true_means <- cbind(plogis(-.5 + base_predictor - contrast / 2),
                          plogis(.5 + base_predictor + contrast / 2))
      projection_fits <- lapply(1:2, function(arm) fit_projection(true_means[, arm], weights, quasibinomial()))
      projected_means <- do.call(cbind, lapply(projection_fits, `[[`, "prediction"))
      true_tate <- sum(weights * (true_means[, 2] - true_means[, 1]))
      predicted_contrast <- projected_means[, 2] - projected_means[, 1]
      common_variance <- sum(weights * (predicted_contrast - sum(weights * predicted_contrast))^2)
      private_variances <- numeric(2)
      candidate_errors <- numeric(2)
      maximum_gradient <- max(vapply(projection_fits, `[[`, numeric(1), "gradient"))
      predictor_bound <- max(vapply(projection_fits, `[[`, numeric(1), "predictor_bound"))
      for (site in 1:2) {
        shift <- (site / 2) * c(.25, -.20, .15, -.10)
        assignment_contrast <- c(.35, -.25, .125, -.0625) / 2
        slopes <- cbind(shift - assignment_contrast, shift + assignment_contrast)
        exponential <- exp(assignment_basis %*% slopes)
        joint <- exponential / sum(weights * rowSums(exponential))
        corrections <- second_moments <- numeric(2)
        for (arm in 1:2) {
          if (scenario == "C3") {
            derivative <- true_means[, arm] * (1 - true_means[, arm])
            tilt <- fit_projection(1 / joint[, arm], weights * joint[, arm] * derivative, quasipoisson())
            weight <- tilt$prediction
            maximum_gradient <- max(maximum_gradient, tilt$gradient)
            predictor_bound <- max(predictor_bound, tilt$predictor_bound)
          } else weight <- 1 / joint[, arm]
          residual_mean <- true_means[, arm] - projected_means[, arm]
          corrections[arm] <- sum(weights * joint[, arm] * weight * residual_mean)
          second_moments[arm] <- sum(weights * joint[, arm] * weight^2 *
            (true_means[, arm] * (1 - true_means[, arm]) + residual_mean^2))
        }
        candidate_errors[site] <- sum(weights * predicted_contrast) + diff(corrections) - true_tate
        private_variances[site] <- sum(second_moments) - diff(corrections)^2
      }
      stopifnot(max(abs(candidate_errors)) < 1e-9, all(private_variances > 0), predictor_bound < 12)
      rows[[length(rows) + 1L]] <- data.frame(config = scenario, effect_heterogeneity = heterogeneity,
        target_tate = true_tate,
        conditional_tate_sd = sqrt(sum(weights * ((true_means[, 2] - true_means[, 1]) - true_tate)^2)),
        common_target_variance_at_n1000 = common_variance / 1000,
        source_correlation = common_variance / sqrt(prod(common_variance + private_variances)),
        mean_target_variance_fraction = mean(common_variance / (common_variance + private_variances)),
        maximum_candidate_bias = max(abs(candidate_errors)), maximum_score_gradient = maximum_gradient,
        maximum_fitted_predictor_bound = predictor_bound)
    }
  }
  do.call(rbind, rows)
}

coarse <- evaluate_design(24L)
fine <- evaluate_design(32L)
metrics <- c("target_tate", "conditional_tate_sd", "common_target_variance_at_n1000",
             "source_correlation", "mean_target_variance_fraction")
difference <- max(abs(as.matrix(coarse[metrics]) - as.matrix(fine[metrics])))
write.csv(fine, file.path(output, "population_sensitivity.csv"), row.names = FALSE)
writeLines(c(paste("Maximum 24/32 quadrature difference:", format(difference, digits = 12)),
  paste("Contrast direction:", paste(contrast_direction, collapse = ", ")),
  "Prospective population calculation only. Existing bounded_joint_v3 data and defaults are unchanged.",
  "Population score roots and candidate identification are checked; no finite-sample coverage claim."),
  file.path(output, "checks.txt"))
stopifnot(difference < 1e-8)
print(fine)
writeLines("POPULATION_HETEROGENEITY_PROBE_PASSED", file.path(output, "COMPLETE"))
