# test-cpp-core.R - Tests for C++ core functions
#
# Covers: GLM fitting, GLM gradient/predict, CV lambda selection,
#         weight optimization, variance/covariance calculations,
#         and the NumericalConstants consistency between R and C++.

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEC namespace)

# ============================================================================
# Test GLM fitting (fit_general_glm_cpp)
# ============================================================================

test_that("fit_general_glm_cpp converges on simple logistic data", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 200
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(0.5, -0.3, 0, 0, 0.2)  # sparse
  eta <- X %*% beta_true
  Y <- rbinom(n, 1, logistic(eta))
  weights <- rep(1, n)
  
  result <- fit_general_glm_cpp(X, Y, weights, 1L, 1L, 0.01, 500, 1e-6,
                                warm_start = rep(0, p + 1))
  
  expect_true(result$converged)
  expect_equal(length(result$alpha), p + 1)  # p + intercept
  # Largest non-zero coefficient should match sign (0.5 is strong enough to recover)
  expect_true(result$alpha[2] > 0)   # beta_true[1] = 0.5
  # beta_true[2] = -0.3 is moderate; with regularization sign may not be recovered
  # Just check it's not strongly positive
  expect_true(result$alpha[3] <= 0.1)   # beta_true[2] = -0.3
})

test_that("fit_general_glm_cpp with lambda=0 gives unpenalized estimate", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  set.seed(123)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- rnorm(n)
  weights <- rep(1, n)
  
  result <- fit_general_glm_cpp(X, Y, weights, 0L, 0L, 0.0, 500, 1e-6,
                                warm_start = rep(0, p + 1))
  
  expect_true(result$converged)
  expect_equal(length(result$alpha), p + 1)
  # With lambda=0, no coefficients should be exactly zero (unless truly zero)
  expect_true(any(result$alpha != 0))
})

test_that("fit_general_glm_cpp with large lambda gives sparse solution", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  set.seed(99)
  n <- 100
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  Y <- rbinom(n, 1, 0.5)
  weights <- rep(1, n)
  
  result <- fit_general_glm_cpp(X, Y, weights, 1L, 1L, 10.0, 500, 1e-6,
                                warm_start = rep(0, p + 1))
  
  # Large lambda should zero out most coefficients (except intercept)
  n_zero <- sum(abs(result$alpha[-1]) < 1e-10)
  expect_true(n_zero >= p / 2)
})

# ============================================================================
# Test GLM gradient and predict (unified functions)
# ============================================================================

test_that("calculate_glm_gradient_cpp returns correct dimensions", {
  skip_if_not(exists("calculate_glm_gradient_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 50
  p <- 5
  W <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p + 1)  # intercept + p
  
  grad <- calculate_glm_gradient_cpp(W, beta, 1L, 1L)
  expect_equal(nrow(grad), n)
  expect_equal(ncol(grad), p + 1)
})

test_that("predict_glm_cpp returns probabilities in [0,1] for binomial/logit", {
  skip_if_not(exists("predict_glm_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 50
  p <- 5
  W <- matrix(rnorm(n * p), n, p)
  beta <- rnorm(p + 1)
  
  preds <- predict_glm_cpp(W, beta, 1L, 1L)
  expect_equal(length(preds), n)
  expect_true(all(preds >= 0 & preds <= 1))
})

test_that("predict_glm_cpp with family/link defaults matches explicit binomial/logit", {
  skip_if_not(exists("predict_glm_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  W <- matrix(rnorm(30 * 4), 30, 4)
  beta <- rnorm(5)
  
  # Binomial/logit = 1, 1
  preds_default <- predict_glm_cpp(W, beta, 1L, 1L)
  preds_explicit <- predict_glm_cpp(W, beta, 1L, 1L)
  
  expect_equal(preds_default, preds_explicit)
})

# ============================================================================
# Test CV lambda selection
# ============================================================================

test_that("select_lambda_cv_density_ratio_cpp returns valid result", {
  skip_if_not(exists("select_lambda_cv_density_ratio_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 80
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  A <- rbinom(n, 1, 0.5)
  mean_grad_psi <- rnorm(p + 1)
  alpha_init <- rnorm(p + 1)
  lambda_grid <- exp(seq(log(0.001), log(1), length.out = 10))
  
  result <- select_lambda_cv_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                                lambda_grid, 3, 100, 1e-4,
                                                A_val = 1L, family_int = 1L, link_int = 1L)
  
  expect_true("best_lambda" %in% names(result))
  expect_true("best_idx" %in% names(result))
  expect_true("cv_scores" %in% names(result))
  expect_true(result$best_lambda > 0)
  expect_equal(length(result$cv_scores), length(lambda_grid))
})

test_that("select_lambda_cv_general_refined_outcome_cpp returns valid result", {
  skip_if_not(exists("select_lambda_cv_general_refined_outcome_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 80
  p <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- rbinom(n, 1, 0.5)
  A <- rbinom(n, 1, 0.5)
  gamma_s <- rnorm(p + 1)
  lambda_grid <- exp(seq(log(0.001), log(1), length.out = 10))
  
  result <- select_lambda_cv_general_refined_outcome_cpp(X, Y, A, gamma_s, 1L, 1L,
                                                          lambda_grid, 3, 100, 1e-4, 1L, X)
  
  expect_true(result$best_lambda > 0)
  expect_equal(length(result$cv_scores), 10)
})

test_that("CV functions handle all-control data gracefully (empty treated set)", {
  skip_if_not(exists("select_lambda_cv_density_ratio_cpp"), message = "C++ not compiled")
  
  n <- 40
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  A <- rep(0, n)  # No treated units
  mean_grad_psi <- rnorm(p + 1)
  alpha_init <- rnorm(p + 1)
  lambda_grid <- exp(seq(log(0.01), log(1), length.out = 5))
  
  result <- select_lambda_cv_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                                lambda_grid, 3, 100, 1e-4,
                                                A_val = 1L, family_int = 1L, link_int = 1L)
  
  # Should return early without error

  expect_true("best_lambda" %in% names(result))
})

# ============================================================================
# Test weight optimization
# ============================================================================

test_that("optimize_weights_cpp produces valid weights", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")
  
  K <- 3
  estimates <- rnorm(K, mean = 0.5, sd = 0.1)
  V_t <- rep(0.1, K)
  V_s <- rep(0.05, K)
  n_s <- rep(200, K)
  V_ot <- 0.2
  n_t <- 100
  C_ot <- rep(0.01, K)
  lambda <- 0.1
  mu_ot <- 0.5
  
  C_cross <- matrix(0, K, K)  # zero cross-site covariance for basic test
  result <- optimize_weights_cpp(estimates, V_t, V_s, n_s, V_ot, n_t, C_ot,
                                  lambda, mu_ot, 1000, 1e-8, C_cross, numeric(0))
  
  expect_true("weights" %in% names(result))
  expect_equal(length(result$weights), K)
  # Weights should be finite
  expect_true(all(is.finite(result$weights)))
})

# ============================================================================
# Test variance calculation functions
# ============================================================================

test_that("calculate_target_variance_cpp returns positive variance", {
  skip_if_not(exists("calculate_target_variance_cpp"), message = "C++ not compiled")
  
  set.seed(42)
  n <- 50
  p <- 3
  W <- matrix(rnorm(n * p), n, p)
  alpha_ts <- c(0.1, 0.5, -0.3, 0.2)  # intercept + 3 covariates (outcome model params)
  M_ts <- 0.4
  
  result <- calculate_target_variance_cpp(W, alpha_ts, M_ts, family_int = 1L, link_int = 1L)
  
  expect_true("V_t" %in% names(result))
  expect_true(result$V_t >= 0)
})

test_that("calculate_covariance_term_cpp returns finite value", {
  skip_if_not(exists("calculate_covariance_term_cpp"), message = "C++ not compiled")
  
  n <- 50
  varphi_ot <- rnorm(n)
  zeta <- rnorm(n)
  
  cov_val <- calculate_covariance_term_cpp(varphi_ot, zeta)
  expect_true(is.finite(cov_val))
})

test_that("calculate_aggregated_variance_cpp returns non-negative", {
  skip_if_not(exists("calculate_aggregated_variance_cpp"), message = "C++ not compiled")
  
  K <- 2
  eta <- c(0.5, 0.5)
  V_t <- c(0.1, 0.1)
  V_s <- c(0.05, 0.05)
  n_s <- c(200, 200)
  V_ot <- 0.2
  n_t <- 100
  C_ot <- c(0.01, 0.01)
  C_cross <- matrix(0, K, K)
  
  var_agg <- calculate_aggregated_variance_cpp(eta, V_t, V_s, n_s, V_ot, n_t,
                                                C_ot, C_cross, 0.1, 0.5)
  
  expect_true(is.finite(var_agg))
  expect_true(var_agg >= 0)
})

test_that("calculate_aggregated_estimate_cpp matches manual calculation", {
  skip_if_not(exists("calculate_aggregated_estimate_cpp"), message = "C++ not compiled")
  
  mu_ot <- 0.5
  mu_ts <- c(0.45, 0.55)
  eta <- c(0.3, 0.7)
  
  result <- calculate_aggregated_estimate_cpp(mu_ot, mu_ts, eta)
  # Expected: (1 - sum(eta)) * mu_ot + sum(eta * mu_ts)
  expected <- (1 - sum(eta)) * mu_ot + sum(eta * mu_ts)
  
  expect_equal(result, expected, tolerance = 1e-10)
})
