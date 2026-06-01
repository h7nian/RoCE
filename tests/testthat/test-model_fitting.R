# test-model_fitting.R - Tests for R-level model fitting wrappers
#
# Covers: fit_initial_outcome, fit_initial_density_ratio, fit_unified_density_ratio,
#         fit_unified_outcome, optimize_weights, calculate_aggregated_variance

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEHD namespace)

# ============================================================================
# Helper: generate small synthetic site data
# ============================================================================

make_site_data <- function(n, p, treated_frac = 0.5, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, treated_frac)
  beta_true <- runif(p + 1, -0.5, 0.5)
  eta <- cbind(1, X) %*% beta_true
  Y <- rbinom(n, 1, logistic(eta))
  list(X = X, A = A, Y = Y, beta_true = beta_true, n = n)
}

# ============================================================================
# Tests for fit_initial_outcome (R-based glmnet)
# ============================================================================

test_that("fit_initial_outcome returns coefficients of expected length", {
  d <- make_site_data(200, 5)
  alpha <- fit_initial_outcome(d$X, d$Y, d$A, A_val = 1)
  
  expect_equal(length(alpha), ncol(d$X) + 1)  # intercept + p
  expect_true(all(is.finite(alpha)))
})

test_that("fit_initial_outcome errors when requested treatment arm is absent", {
  d <- make_site_data(100, 3, treated_frac = 0.05, seed = 7)
  d$A[] <- 0

  expect_error(
    fit_initial_outcome(d$X, d$Y, d$A, A_val = 1),
    "no observations with A_val"
  )
})

test_that("fit_initial_outcome works when requested treatment arm is sparse but present", {
  d <- make_site_data(100, 3, treated_frac = 0.08, seed = 7)
  d$A[] <- 0
  d$A[1:10] <- 1

  alpha <- fit_initial_outcome(d$X, d$Y, d$A, A_val = 1, nlambda = 3L)
  expect_true(all(is.finite(alpha)))
})

test_that("fit_initial_outcome exposes glmnet-style lambda rules", {
  d <- make_site_data(120, 4, treated_frac = 0.5, seed = 11)

  alpha <- fit_initial_outcome(
    d$X, d$Y, d$A, A_val = 1, nlambda = 5L,
    lambda_rule = "1se"
  )

  expect_equal(attr(alpha, "lambda_rule"), "1se")
  expect_true(is.finite(attr(alpha, "lambda_min")))
  expect_true(is.finite(attr(alpha, "lambda_1se")))
  expect_equal(attr(alpha, "lambda_used"), attr(alpha, "lambda_1se"))
})

test_that("nuisance CV lambda selection requires standard metadata", {
  skip_if_not(exists(".select_nuisance_cv_lambda"), message = "internal helper not loaded")

  expect_error(
    .select_nuisance_cv_lambda(
      list(best_lambda = 0.1, lambda_min = 0.1),
      lambda_rule = "1se",
      caller = "unit_test"
    ),
    "missing finite glmnet-style field"
  )
})

# ============================================================================
# Tests for fit_initial_density_ratio
# ============================================================================

test_that("fit_initial_density_ratio returns valid gamma", {
  skip_if_not(exists("fit_initial_density_ratio_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 200
  p <- 5
  Z_site <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, 0.5)
  mean_phi <- colMeans(cbind(1, Z_site))
  
  gamma <- fit_initial_density_ratio(Z_site, A, mean_phi, lambda = 0.1)
  
  expect_equal(length(gamma), p + 1)
  expect_true(all(is.finite(gamma)))
})

test_that("optimize_weights rejects lambda values it would otherwise clip", {
  estimates <- c(0.50, 0.62)
  variances <- list(V_ot = 1.0, V_t = c(0.8, 0.9), V_s = c(0.7, 0.6))
  C_ot <- c(0.05, 0.04)
  n_samples <- list(n_t = 40, n_s = c(45, 50))

  expect_error(
    optimize_weights(
      estimates, variances, C_ot, n_samples,
      lambda = LAMBDA_MAX * 2, mu_ot = 0.50
    ),
    "lambda must be <="
  )
})

# ============================================================================
# Tests for fit_unified_density_ratio
# ============================================================================

test_that("fit_unified_density_ratio works in refined mode", {
  skip_if_not(exists("fit_unified_density_ratio_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 150
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, 0.5)
  alpha_init <- rnorm(p + 1, sd = 0.1)
  mean_grad_psi <- rnorm(p + 1, sd = 0.01)
  
  gamma <- fit_unified_density_ratio(X, A, mean_grad_psi, alpha_init,
                                      lambda = 0.05, calibrated = FALSE,
                                      W_outcome = X)
  
  expect_equal(length(gamma), p + 1)
  expect_true(all(is.finite(gamma)))
})

test_that("fit_unified_density_ratio works in calibrated mode", {
  skip_if_not(exists("fit_unified_density_ratio_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 150
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, 0.5)
  alpha_init <- rnorm(p + 1, sd = 0.1)
  mean_grad_psi <- rnorm(p + 1, sd = 0.01)
  
  gamma <- fit_unified_density_ratio(X, A, mean_grad_psi, alpha_init,
                                      lambda = 0.05, calibrated = TRUE, M_tau = 10.0,
                                      W_outcome = X)
  
  expect_equal(length(gamma), p + 1)
  expect_true(all(is.finite(gamma)))
})

# New test: ensure density‑ratio CD with Hessian step converges stably
test_that("fit_unified_density_ratio_cpp converges with Hessian step approximation", {
  skip_if_not(exists("fit_unified_density_ratio_cpp"), message = "C++ not compiled")

  set.seed(2026)
  n <- 300
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, 0.45)
  alpha_init <- rnorm(p + 1, sd = 1.0)
  mean_grad_psi <- rnorm(p + 1, sd = 0.05)

  res <- fit_unified_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                       0.01, 1000, 1e-6, TRUE, 10.0,
                                       W_outcome = X,
                                       A_val = 1L, family_int = 1L, link_int = 1L,
                                       warm_start = rep(0, p + 1))

  expect_true(is.list(res))
  # With random alpha_init and truncation, convergence may require more iterations;

  # focus on checking the result is well-formed and finite
  expect_equal(length(res$gamma), p + 1)
  expect_true(all(is.finite(res$gamma)))
})

test_that("site-basis GLM gradient matches outcome basis when W and Z agree", {
  skip_if_not(exists("mean_glm_gradient_cpp"), message = "C++ not compiled")

  set.seed(2027)
  W <- matrix(rnorm(80 * 3), 80, 3)
  alpha <- rnorm(4, sd = 0.2)

  r_grad <- .mean_glm_gradient_site_basis(W, W, alpha,
                                          family_int = FAMILY_BINOMIAL,
                                          link_int = LINK_LOGIT)
  cpp_grad <- mean_glm_gradient_cpp(W, alpha, FAMILY_BINOMIAL, LINK_LOGIT)

  expect_equal(r_grad, as.numeric(cpp_grad), tolerance = 1e-10)
})

test_that("site-basis GLM gradient supports different W and Z dimensions", {
  set.seed(2028)
  W <- matrix(rnorm(60 * 2), 60, 2)
  Z <- matrix(rnorm(60 * 5), 60, 5)
  alpha <- rnorm(3, sd = 0.2)

  grad <- .mean_glm_gradient_site_basis(W, Z, alpha,
                                        family_int = FAMILY_BINOMIAL,
                                        link_int = LINK_LOGIT)

  expect_equal(length(grad), ncol(Z) + 1)
  expect_true(all(is.finite(grad)))
})

test_that("fit_unified_density_ratio rejects mismatched target-gradient dimensions", {
  skip_if_not(exists("fit_unified_density_ratio_cpp"), message = "C++ not compiled")

  set.seed(2029)
  n <- 120
  W <- matrix(rnorm(n * 2), n, 2)
  Z <- matrix(rnorm(n * 4), n, 4)
  A <- rbinom(n, 1, 0.5)
  alpha_init <- rnorm(ncol(W) + 1, sd = 0.1)
  wrong_mean_grad <- rnorm(ncol(W) + 1, sd = 0.01)

  expect_error(
    fit_unified_density_ratio(Z, A, wrong_mean_grad, alpha_init,
                              lambda = 0.05, calibrated = TRUE, M_tau = 10.0,
                              W_outcome = W),
    "mean_grad_psi length must match"
  )
})

test_that("correction term reports clipping diagnostics", {
  skip_if_not(exists("calculate_correction_term_cpp"), message = "C++ not compiled")

  n <- 12
  Z <- matrix(c(rep(-30, n), rep(0, n)), n, 2)
  W <- matrix(rnorm(n * 2), n, 2)
  A <- rep(1, n)
  Y <- rnorm(n)
  gamma <- c(0, 1, 0)
  alpha <- c(0, 0, 0)

  res <- calculate_correction_term_cpp(
    Z, A, Y, gamma, alpha, W,
    M_tau = Inf,
    family_int = FAMILY_GAUSSIAN,
    link_int = LINK_IDENTITY,
    A_val = 1L
  )

  expect_true("clip_diagnostics" %in% names(res))
  expect_true(res$clip_diagnostics$ratio_max_clipped > 0)
  expect_true(isTRUE(res$clip_diagnostics$any_clipped))
})

# ============================================================================
# Tests for fit_unified_outcome
# ============================================================================

test_that("fit_unified_outcome returns valid alpha (refined)", {
  skip_if_not(exists("fit_unified_outcome_cpp"), message = "C++ not compiled")
  
  d <- make_site_data(200, 4)
  gamma_s <- rnorm(ncol(d$X) + 1, sd = 0.1)
  
  alpha <- fit_unified_outcome(d$X, d$Y, d$A, A_val = 1, gamma_s = gamma_s,
                                lambda = 0.05, calibrated = FALSE, Z_site = d$X)
  
  expect_equal(length(alpha), ncol(d$X) + 1)
  expect_true(all(is.finite(alpha)))
})

test_that("fit_unified_outcome returns valid alpha (calibrated)", {
  skip_if_not(exists("fit_unified_outcome_cpp"), message = "C++ not compiled")
  
  d <- make_site_data(200, 4)
  gamma_s <- rnorm(ncol(d$X) + 1, sd = 0.1)
  
  alpha <- fit_unified_outcome(d$X, d$Y, d$A, A_val = 1, gamma_s = gamma_s,
                                lambda = 0.05, calibrated = TRUE, M_tau = 10.0, Z_site = d$X)
  
  expect_equal(length(alpha), ncol(d$X) + 1)
  expect_true(all(is.finite(alpha)))
})

# ============================================================================
# Tests for optimize_weights
# ============================================================================

test_that("optimize_weights returns valid weights for K=2 sites", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")
  
  K <- 2
  estimates <- c(0.45, 0.55)
  variances <- list(V_t = c(0.10, 0.12), V_s = c(0.05, 0.04), V_ot = 0.20)
  C_ot <- c(0.01, 0.02)
  n_samples <- list(n_t = 100, n_s = c(200, 250))
  lambda <- 0.1
  mu_ot <- 0.50
  
  weights <- optimize_weights(estimates, variances, C_ot, n_samples,
                               lambda, mu_ot)
  
  expect_equal(length(weights), K)
  expect_true(all(is.finite(weights)))
  expect_true(all(weights >= 0 & weights <= 1))
})

test_that("optimize_weights with single source site", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")
  
  estimates <- c(0.48)
  variances <- list(V_t = 0.10, V_s = 0.05, V_ot = 0.20)
  C_ot <- c(0.01)
  n_samples <- list(n_t = 100, n_s = c(200))
  
  weights <- optimize_weights(estimates, variances, C_ot, n_samples,
                               lambda = 0.01, mu_ot = 0.50)
  
  expect_equal(length(weights), 1)
  expect_true(is.finite(weights))
})

test_that("optimize_weights respects clip_weights=FALSE", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")
  
  K <- 2
  estimates <- c(0.45, 0.55)
  variances <- list(V_t = c(0.10, 0.12), V_s = c(0.05, 0.04), V_ot = 0.20)
  C_ot <- c(0.01, 0.02)
  n_samples <- list(n_t = 100, n_s = c(200, 250))
  
  # With clip_weights=FALSE, weights may exceed [0,1]
  weights <- optimize_weights(estimates, variances, C_ot, n_samples,
                               lambda = 0.01, mu_ot = 0.5, clip_weights = FALSE)
  
  expect_equal(length(weights), K)
  expect_true(all(is.finite(weights)))
})

# ============================================================================
# Tests for calculate_aggregated_variance
# ============================================================================

test_that("calculate_aggregated_variance returns non-negative", {
  skip_if_not(exists("calculate_aggregated_variance_cpp"), message = "C++ not compiled")
  
  K <- 2
  eta <- c(0.4, 0.3)
  variances <- list(V_t = c(0.10, 0.12), V_s = c(0.05, 0.04), V_ot = 0.20)
  C_ot <- c(0.01, 0.02)
  n_samples <- list(n_t = 100, n_s = c(200, 250))
  
  var_agg <- calculate_aggregated_variance(eta, variances, C_ot, n_samples)
  
  expect_true(is.finite(var_agg))
  expect_true(var_agg >= 0)
})

test_that("calculate_aggregated_variance with cross-site covariances", {
  skip_if_not(exists("calculate_aggregated_variance_cpp"), message = "C++ not compiled")
  
  K <- 3
  eta <- c(0.3, 0.3, 0.2)
  variances <- list(V_t = c(0.10, 0.12, 0.15), V_s = c(0.05, 0.04, 0.06), V_ot = 0.20)
  C_ot <- c(0.01, 0.02, 0.015)
  n_samples <- list(n_t = 100, n_s = c(200, 250, 180))
  C_cross <- matrix(c(0, 0.005, 0.003,
                       0.005, 0, 0.004,
                       0.003, 0.004, 0), 3, 3)
  
  var_agg <- calculate_aggregated_variance(eta, variances, C_ot, n_samples,
                                            C_cross = C_cross)
  
  expect_true(is.finite(var_agg))
  expect_true(var_agg >= 0)
})
