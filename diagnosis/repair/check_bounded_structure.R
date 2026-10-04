#!/usr/bin/env Rscript
# Population quadrature audit of overlap and sparse nuisance targets. All
# optimization here has five coefficients; no high-dimensional Monte Carlo.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: check_bounded_structure.R LIBRARY NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
quadrature <- RoCE:::.bounded_quadrature(16L)
X <- quadrature$X
weights <- quadrature$weights
Z <- cbind(1, X)
bound <- RoCE:::.bounded_covariate_spec()$bound

fit_population <- function(kind, multiplier, response = NULL, linear = NULL) {
  mass <- weights * multiplier
  beta <- numeric(ncol(Z))
  if (kind == "logistic") beta[1L] <- qlogis(sum(mass * response) / sum(mass)) else
    beta[1L] <- log(sum(mass) / linear[1L])
  evaluate <- function(beta) {
    eta <- drop(Z %*% beta)
    if (kind == "logistic") {
      probability <- plogis(eta)
      list(value = sum(mass * (pmax(eta, 0) + log1p(exp(-abs(eta))) - response * eta)),
        gradient = drop(crossprod(Z, mass * (probability - response))),
        hessian = crossprod(Z, Z * (mass * probability * (1 - probability))))
    } else {
      tilted <- mass * exp(-eta)
      list(value = sum(linear * beta) + sum(tilted),
        gradient = linear - drop(crossprod(Z, tilted)), hessian = crossprod(Z, Z * tilted))
    }
  }
  state <- evaluate(beta)
  for (iteration in seq_len(100L)) {
    if (max(abs(state$gradient)) < 1e-10) break
    direction <- solve(state$hessian, state$gradient)
    step <- 1
    repeat {
      proposed <- beta - step * direction
      next_state <- evaluate(proposed)
      if (is.finite(next_state$value) && next_state$value <= state$value + 1e-13) break
      step <- step / 2
      if (step < 1e-10) stop("Population line search failed.")
    }
    beta <- proposed
    state <- next_state
  }
  stopifnot(max(abs(state$gradient)) < 1e-9)
  list(coefficients = beta, gradient = max(abs(state$gradient)),
       predictor_bound = abs(beta[1L]) + bound * sum(abs(beta[-1L])))
}

overlap <- population <- list()
add_population <- function(fit, config, site, arm, delta, xi, model) {
  population[[length(population) + 1L]] <<- data.frame(config, site, arm,
    covariate_shift = delta, source_treatment_scale = xi, model,
    max_gradient = fit$gradient, predictor_bound = fit$predictor_bound)
  fit$coefficients
}
for (config in c("C1", "C2", "C3")) {
  strengths <- RoCE:::.face_misspecification_strengths(config, .75)
  assignment <- RoCE:::.bounded_active_features(X, strengths$propensity)
  outcome <- RoCE:::.bounded_active_features(X, strengths$outcome)
  propensity <- plogis(drop(assignment %*% c(.35, -.25, .125, -.0625)))
  outcome_predictor <- drop(outcome %*% c(.25, .20, .10, .05))
  for (delta in c(0, .5, 1, 2)) for (xi in c(0, 1)) {
    for (arm in 0:1) {
      probability <- if (arm == 1L) propensity else 1 - propensity
      mean_y <- plogis(-.5 + arm + outcome_predictor)
      initial_or <- fit_population("logistic", probability, mean_y)
      initial_ps <- fit_population("tilt", probability,
        linear = drop(crossprod(Z, weights * (1 - probability))))
      initial_mean <- plogis(drop(Z %*% initial_or$coefficients))
      h <- initial_mean * (1 - initial_mean)
      calibrated_ps <- fit_population("tilt", probability * h,
        linear = drop(crossprod(Z, weights * (1 - probability) * h)))
      calibrated_or <- fit_population("logistic", probability * exp(-drop(Z %*% initial_ps$coefficients)), mean_y)
      for (name in c("initial_or", "initial_ps", "calibrated_or", "calibrated_ps")) {
        add_population(get(name), config, "t", arm, delta, xi, name)
      }
      for (j in 1:4) {
        common <- delta * (j / 4) * c(.25, -.20, .15, -.10)
        contrast <- xi * c(.35, -.25, .125, -.0625) / 2
        slopes <- cbind(common - contrast, common + contrast)
        eta <- assignment %*% slopes
        exponential <- exp(eta)
        normalizer <- sum(weights * rowSums(exponential))
        joint <- exponential / normalizer
        density <- rowSums(joint)
        joint_arm <- joint[, arm + 1L]
        q <- 1 / joint_arm
        arm_probability <- sum(weights * joint_arm)
        mean_x <- colSums(X * (weights * density))
        variance_x <- colSums(X^2 * (weights * density)) - mean_x^2
        gram <- crossprod(Z, Z * (weights * joint_arm))
        eta_range <- RoCE:::.bounded_predictor_range(slopes[, arm + 1L], strengths$propensity)
        ps_range <- plogis(RoCE:::.bounded_predictor_range(slopes[, 2L] - slopes[, 1L], strengths$propensity))
        overlap[[length(overlap) + 1L]] <- data.frame(config, site = paste0("s", j), arm,
          covariate_shift = delta, source_treatment_scale = xi,
          density_normalization_error = abs(sum(weights * density) - 1),
          source_arm_probability = arm_probability, propensity_min = ps_range[1L], propensity_max = ps_range[2L],
          max_abs_smd = max(abs(mean_x) / sqrt((1 + variance_x) / 2)),
          hellinger_squared = 1 - sum(weights * sqrt(density)),
          transport_ess_fraction = 1 / sum(weights / density),
          oracle_arm_effective_n = 1000 / sum(weights * q),
          joint_weight_min = normalizer * exp(-eta_range[2L]),
          joint_weight_max = normalizer * exp(-eta_range[1L]),
          population_gram_min_eigenvalue = min(arm_probability, min(eigen(gram, symmetric = TRUE, only.values = TRUE)$values)))
        initial_weight <- fit_population("tilt", joint_arm,
          linear = drop(crossprod(Z, weights)))
        final_weight <- fit_population("tilt", joint_arm * h,
          linear = drop(crossprod(Z, weights * h)))
        add_population(initial_weight, config, paste0("s", j), arm, delta, xi, "initial_weight")
        add_population(final_weight, config, paste0("s", j), arm, delta, xi, "calibrated_weight")
        for (rho in c(0, 2.5)) {
          source_mean <- plogis(-.5 + arm + outcome_predictor + rho)
          final_or <- fit_population("logistic", joint_arm * exp(-drop(Z %*% initial_weight$coefficients)), source_mean)
          add_population(final_or, config, paste0("s", j), arm, delta, xi, paste0("calibrated_or_intercept_shift", rho))
        }
      }
    }
    cat(config, "covariate_shift", delta, "treatment_scale", xi, "completed\n")
  }
}
overlap <- do.call(rbind, overlap)
population <- do.call(rbind, population)
write.csv(overlap, file.path(output, "overlap.csv"), row.names = FALSE)
write.csv(population, file.path(output, "population_fits.csv"), row.names = FALSE)
stopifnot(max(overlap$density_normalization_error) < 1e-10,
          min(overlap$population_gram_min_eigenvalue) > 0,
          min(overlap$propensity_min) > 0, max(overlap$propensity_max) < 1)
writeLines(c("Population quadrature and five-coefficient optimization only; not empirical coverage.",
  "Inactive coordinates are independent, centered and have unit variance in every site/arm.",
  "This establishes their zero population score structurally; no 200-dimensional fit was run.",
  "Reported predictor bounds are for fitted population projections on this quadrature grid.",
  sprintf("Largest population predictor bound: %.8f", max(population$predictor_bound)),
  "Numerical population checks complement, rather than prove, asymptotic inference conditions."),
  file.path(output, "scope.txt"))
writeLines("PASSED: structural numerical checks", file.path(output, "status.txt"))
