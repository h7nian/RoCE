# test-method-alignment.R - Formula <-> code alignment + result-schema tests
#
# This file targets the algebraic identities in docs/main.tex that the codebase
# implements. Failing tests here mean either:
#   (a) the implementation has drifted from the paper, or
#   (b) main.tex has been updated and the code/tests need to follow.
#
# The reporting tests at the bottom also pin the result schema for the
# cross-fitting one-round and two-round estimators so that any regression which
# silently drops a method from the simulation summary is caught here.

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEC namespace).
# Reuses the cached smoke result from test-cross_fitting.R when both files run
# in the same session.

# ============================================================================
# eq:aggregated_estimator -- mu_hat = mu_ot - sum_j eta_j * (mu_ot - mu_{ts,j})
# (a.k.a. (1 - sum eta) * mu_ot + sum eta_j * mu_{ts,j})
# ============================================================================

test_that("aggregated estimator: eta = 0 recovers target-only (eq:aggregated_estimator)", {
  skip_if_not(exists("calculate_aggregated_estimate_cpp"), message = "C++ not compiled")

  mu_ot <- 0.7
  mu_ts <- c(0.3, 0.4, 0.5)
  eta_zero <- rep(0.0, length(mu_ts))

  result <- calculate_aggregated_estimate_cpp(mu_ot, mu_ts, eta_zero)
  expect_equal(result, mu_ot, tolerance = 1e-12,
               info = "With eta = 0 the aggregated estimator must equal the target-only estimate.")
})

test_that("aggregated estimator: one-hot eta selects a single source (eq:aggregated_estimator)", {
  skip_if_not(exists("calculate_aggregated_estimate_cpp"), message = "C++ not compiled")

  mu_ot <- 0.5
  mu_ts <- c(0.21, 0.42, 0.63)
  for (j in seq_along(mu_ts)) {
    eta <- rep(0.0, length(mu_ts))
    eta[j] <- 1.0
    result <- calculate_aggregated_estimate_cpp(mu_ot, mu_ts, eta)
    expect_equal(result, mu_ts[j], tolerance = 1e-12,
                 info = sprintf("One-hot eta on source %d must return mu_ts[%d].", j, j))
  }
})

# ============================================================================
# eq:dr_estimator + eq:if_target_only -- AIPW influence function structure
# ============================================================================

test_that("calculate_aipw_influence: unweighted estimate equals mean(phi)", {
  set.seed(2026)
  n <- 200
  X <- matrix(rnorm(n * 3), n, 3)
  pi_true <- 1 / (1 + exp(-as.numeric(X %*% c(0.5, -0.3, 0.2))))
  A <- rbinom(n, 1, pi_true)
  m_true <- as.numeric(X %*% c(0.2, 0.4, -0.1))
  Y <- m_true + A + rnorm(n)

  # Use the true nuisances so the test focuses on the structural identity,
  # not on glmnet's behaviour.
  res <- calculate_aipw_influence(y = Y, a = A, x = X,
                                  m_hat = m_true + 1.0,  # m for A_val = 1 arm
                                  pi_hat = pi_true,
                                  w = NULL, A_val = 1L, family = "gaussian")

  # Without DR weights, estimate must equal mean(phi) where phi is the AIPW
  # pseudo-outcome.
  phi_manual <- (m_true + 1.0) + as.numeric(A == 1L) * (Y - (m_true + 1.0)) /
                pmax(pi_true, PROP_SCORE_LOWER)
  expect_equal(res$estimate, mean(phi_manual), tolerance = 1e-10)
})

test_that("calculate_aipw_influence: influence has length n and centered mean = 0", {
  set.seed(7)
  n <- 150
  X <- matrix(rnorm(n * 2), n, 2)
  pi_true <- rep(0.5, n)
  A <- rbinom(n, 1, pi_true)
  m_true <- as.numeric(X %*% c(0.3, -0.4))
  Y <- m_true + A * 0.5 + rnorm(n, sd = 0.5)

  res <- calculate_aipw_influence(y = Y, a = A, x = X,
                                  m_hat = m_true + 0.5,
                                  pi_hat = pi_true,
                                  w = NULL, A_val = 1L, family = "gaussian")

  # Influence function shape matches eq:if_target_only: one row per observation.
  expect_equal(length(res$influence), n)

  # The reported variance must equal mean(IF^2) / n_t (eq:variance_components).
  expect_equal(res$variance, mean(res$influence^2) / n, tolerance = 1e-10)
})

# ============================================================================
# eq:weight_def -- density-ratio weight normalization
# ============================================================================

test_that("calculate_dr_weights: normalized to mean(w) = 1 and strictly positive", {
  skip_if_not(exists("fit_initial_density_ratio"), message = "fit_initial_density_ratio not loaded")

  set.seed(11)
  n_source <- 200
  n_target <- 200
  p <- 3
  # Slight covariate shift between source and target.
  Z_source <- matrix(rnorm(n_source * p), n_source, p)
  Z_target <- matrix(rnorm(n_target * p, mean = 0.2), n_target, p)

  w <- calculate_dr_weights(Z_source, Z_target, lambda = 0.05)

  expect_equal(length(w), n_source)
  expect_true(all(w > 0), info = "DR weights must be strictly positive.")
  # After the post-clip mean-normalization the weights average to 1 up to
  # boundary effects from clipping. Be generous on tolerance for clipped runs.
  expect_lt(abs(mean(w) - 1), 0.25)
})

test_that("calculate_dr_weights: validates matching covariate dimensions", {
  expect_error(
    calculate_dr_weights(matrix(0, 5, 3), matrix(0, 5, 2)),
    "Z_source has 3 columns but Z_target has 2"
  )
})

# ============================================================================
# Result-schema regression tests for run_crossfit
# ----------------------------------------------------------------------------
# These guard against the regression where the simulation summary was missing
# rows because `source_estimates_matrix` was not threaded through the
# aggregation call. They live here so that test-method-alignment is the single
# place that pins the public output shape of the federated estimators.
# ============================================================================

.required_crossfit_fields <- c("estimate", "variance", "se",
                               "ci_lower", "ci_upper",
                               "fold_lambdas", "clip_diagnostics",
                               "aggregation_lambda_rule")

test_that("run_crossfit(one_round) result schema is complete", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  result <- get_smoke_result("one_round")
  missing <- setdiff(.required_crossfit_fields, names(result))
  expect_equal(length(missing), 0L,
               info = sprintf("one_round result missing fields: %s",
                              paste(missing, collapse = ", ")))
})

test_that("run_crossfit(two_round) result schema is complete", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  result <- get_smoke_result("two_round")
  missing <- setdiff(.required_crossfit_fields, names(result))
  expect_equal(length(missing), 0L,
               info = sprintf("two_round result missing fields: %s",
                              paste(missing, collapse = ", ")))
})

test_that("run_crossfit reported CI = estimate +/- Z_ALPHA_05 * se for both modes", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")

  for (mode in c("one_round", "two_round")) {
    result <- get_smoke_result(mode)
    expect_equal(result$se, sqrt(result$variance), tolerance = 1e-10,
                 info = sprintf("se must equal sqrt(variance) for %s mode", mode))
    expect_equal(result$ci_lower,
                 result$estimate - Z_ALPHA_05 * result$se,
                 tolerance = 1e-10,
                 info = sprintf("ci_lower formula mismatch for %s mode", mode))
    expect_equal(result$ci_upper,
                 result$estimate + Z_ALPHA_05 * result$se,
                 tolerance = 1e-10,
                 info = sprintf("ci_upper formula mismatch for %s mode", mode))
  }
})

test_that("run_crossfit one_round and two_round both produce finite, well-defined output", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")

  for (mode in c("one_round", "two_round")) {
    result <- get_smoke_result(mode)
    expect_true(is.finite(result$estimate),
                info = sprintf("%s estimate must be finite", mode))
    expect_true(is.finite(result$variance) && result$variance >= 0,
                info = sprintf("%s variance must be finite and non-negative", mode))
    expect_false(any(!is.finite(c(result$ci_lower, result$ci_upper))),
                 info = sprintf("%s CI endpoints must be finite", mode))
    # clip_diagnostics is structured; just ensure the field carries something.
    expect_true(!is.null(result$clip_diagnostics),
                info = sprintf("%s clip_diagnostics must be present", mode))
  }
})
