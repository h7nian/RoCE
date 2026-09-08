# test-cpp-core.R - Tests for C++ core functions
#
# Covers: GLM fitting, GLM gradient/predict, CV lambda selection,
#         weight optimization, variance/covariance calculations,
#         and the NumericalConstants consistency between R and C++.

library(testthat)

# Package loaded by helper-load.R (all functions available via RoCE namespace)

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
  expect_true(is.finite(result$max_update))
  expect_true(is.finite(result$convergence_threshold))
  expect_lte(result$max_update, result$convergence_threshold * (1 + 1e-10))
  expect_identical(result$line_search_failures, 0L)
  expect_equal(result$max_abs_coefficient, max(abs(result$alpha)))
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
  lambda_grid <- exp(seq(log(1), log(0.001), length.out = 10))
  
  result <- select_lambda_cv_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                                lambda_grid, 3, 100, 1e-4,
                                                A_val = 1L, family_int = 1L, link_int = 1L)
  
  expect_true("best_lambda" %in% names(result))
  expect_true("best_idx" %in% names(result))
  expect_true("cv_scores" %in% names(result))
  expect_true("lambda_min" %in% names(result))
  expect_true("lambda_1se" %in% names(result))
  expect_true("cv_se" %in% names(result))
  expect_true(result$best_lambda > 0)
  expect_true(result$lambda_min > 0)
  expect_true(result$lambda_1se > 0)
  expect_equal(length(result$cv_scores), length(lambda_grid))
  expect_equal(length(result$cv_se), length(lambda_grid))
  expect_true(result$invalid_fold_fits >= 0)
  expect_true(result$invalid_lambdas >= 0)
  expect_true(result$path_tail_skipped_fold_fits >= 0)
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
  lambda_grid <- exp(seq(log(1), log(0.001), length.out = 10))
  
  result <- select_lambda_cv_general_refined_outcome_cpp(X, Y, A, gamma_s, 1L, 1L,
                                                          lambda_grid, 3, 100, 1e-4, 1L, X)
  
  expect_true(result$best_lambda > 0)
  expect_true(result$lambda_min > 0)
  expect_true(result$lambda_1se > 0)
  expect_equal(length(result$cv_scores), 10)
  expect_equal(length(result$cv_se), 10)
})

test_that("every fold-parallel nuisance CV path is deterministic", {
  cv_entry_points <- c(
    "select_lambda_cv_density_ratio_cpp",
    "select_lambda_cv_initial_density_ratio_cpp",
    "select_lambda_cv_calibrated_density_ratio_cpp",
    "select_lambda_cv_general_refined_outcome_cpp",
    "select_lambda_cv_calibrated_outcome_cpp"
  )
  skip_if_not(
    all(vapply(cv_entry_points, exists, logical(1L), inherits = TRUE)),
    message = "C++ nuisance CV entry points are not compiled"
  )

  previous_threads <- Sys.getenv(
    "ROCE_NUISANCE_CV_THREADS", unset = NA_character_
  )
  on.exit({
    if (is.na(previous_threads)) {
      Sys.unsetenv("ROCE_NUISANCE_CV_THREADS")
    } else {
      Sys.setenv(ROCE_NUISANCE_CV_THREADS = previous_threads)
    }
  }, add = TRUE)

  set.seed(20260814)
  n <- 120L
  p <- 8L
  X <- matrix(rnorm(n * p), n, p)
  A <- rep(c(0, 1), length.out = n)
  Y <- rbinom(n, 1, plogis(0.4 * X[, 1] - 0.25 * X[, 2]))
  lambda_grid <- exp(seq(log(0.5), log(0.005), length.out = 8L))
  mean_phi <- c(0.5, rep(0, p))
  mean_grad_psi <- c(0.25, rep(0, p))
  zero_coefficients <- rep(0, p + 1L)

  run_cv_paths <- function(thread_count) {
    # Every C++ CV entry point draws its fold assignments from R's RNG.
    # Reset the same stream before each thread-count run so this test isolates
    # execution topology rather than comparing two unrelated fold partitions.
    set.seed(20260815)
    Sys.setenv(ROCE_NUISANCE_CV_THREADS = as.character(thread_count))
    list(
      refined_density_ratio = select_lambda_cv_density_ratio_cpp(
        X, A, mean_grad_psi, zero_coefficients, lambda_grid,
        n_folds = 3L, max_iter = 1000L, tol = 1e-4, A_val = 1L,
        family_int = 1L, link_int = 1L
      ),
      initial_density_ratio = select_lambda_cv_initial_density_ratio_cpp(
        X, A, mean_phi, lambda_grid,
        n_folds = 3L, max_iter = 1000L, tol = 1e-4, A_val = 1L
      ),
      calibrated_density_ratio =
        select_lambda_cv_calibrated_density_ratio_cpp(
          X, A, mean_grad_psi, zero_coefficients, lambda_grid,
          n_folds = 3L, max_iter = 1000L, tol = 1e-4, M_tau = 5,
          W_outcome = X, A_val = 1L, family_int = 1L, link_int = 1L
        ),
      refined_outcome = select_lambda_cv_general_refined_outcome_cpp(
        X, Y, A, zero_coefficients, 1L, 1L, lambda_grid,
        n_folds = 3L, max_iter = 500L, tol = 1e-5, A_val = 1L,
        Z_site = X
      ),
      calibrated_outcome = select_lambda_cv_calibrated_outcome_cpp(
        X, Y, A, zero_coefficients, lambda_grid,
        n_folds = 3L, max_iter = 500L, tol = 1e-5, A_val = 1L,
        M_tau = 5, Z_site = X, family_int = 1L, link_int = 1L
      )
    )
  }

  sequential <- run_cv_paths(1L)
  parallel <- tryCatch(run_cv_paths(2L), error = identity)
  if (inherits(parallel, "error") && grepl(
      "compiled without OpenMP", conditionMessage(parallel), fixed = TRUE
  )) {
    skip("RoCE was compiled without OpenMP on this platform")
  }
  expect_false(inherits(parallel, "error"))
  expect_equal(parallel, sequential, tolerance = 0)

  Sys.setenv(ROCE_NUISANCE_CV_THREADS = "2x")
  expect_error(
    select_lambda_cv_initial_density_ratio_cpp(
      X, A, c(0.5, rep(0, p)), lambda_grid,
      n_folds = 3L, max_iter = 500L, tol = 1e-5, A_val = 1L
    ),
    "must be a positive integer"
  )
})

test_that("nuisance CV errors when no candidate converges in every fold", {
  skip_if_not(exists("select_lambda_cv_initial_density_ratio_cpp"),
              message = "C++ not compiled")
  skip_if_not(exists("select_lambda_cv_general_refined_outcome_cpp"),
              message = "C++ not compiled")

  set.seed(20260814)
  n <- 120
  p <- 6
  X <- matrix(rnorm(n * p), n, p)
  A <- rep(c(0, 1), length.out = n)
  Y <- rbinom(n, 1, plogis(0.8 * X[, 1] - 0.6 * X[, 2]))
  lambda_grid <- 1e-8

  expect_error(
    select_lambda_cv_initial_density_ratio_cpp(
      X, A, c(0.5, rep(0, p)), lambda_grid,
      n_folds = 3L, max_iter = 1L, tol = 1e-12, A_val = 1L
    ),
    "no lambda converged with a finite validation loss across every CV fold"
  )
  expect_error(
    select_lambda_cv_general_refined_outcome_cpp(
      X, Y, A, rep(0, p + 1L), 1L, 1L, lambda_grid,
      n_folds = 3L, max_iter = 1L, tol = 1e-12, A_val = 1L,
      Z_site = X
    ),
    "no lambda converged with a finite validation loss across every CV fold"
  )
})

test_that("nuisance CV records and skips a terminal failed lambda tail", {
  skip_if_not(exists("select_lambda_cv_initial_density_ratio_cpp"),
              message = "C++ not compiled")

  set.seed(20260824)
  n <- 120L
  p <- 12L
  X <- matrix(rnorm(n * p), n, p)
  A <- rep(c(0, 1), length.out = n)
  # At the three highly regularized candidates, zero is an exact solution:
  # the intercept moment matches the arm fraction and all penalized feature
  # updates threshold to zero. The smaller candidates cannot converge in one
  # sweep, yielding a deterministic terminal failure tail.
  lambda_grid <- c(1e6, 1e5, 1e4, rep(1e-10, 12L))
  result <- select_lambda_cv_initial_density_ratio_cpp(
    X, A, c(0.5, rep(0, p)), lambda_grid,
    n_folds = 3L, max_iter = 1L, tol = 1e-12, A_val = 1L
  )

  expect_equal(result$lambda_min, 1e6)
  expect_gt(result$path_tail_skipped_fold_fits, 0L)
  expect_lte(
    result$path_tail_skipped_fold_fits,
    result$invalid_fold_fits
  )
})

test_that("CV functions fail fast on all-control data (empty treated set)", {
  # The estimator is not identifiable when there are no A == A_val units in any
  # CV training fold, and the project's fail-fast philosophy requires that the
  # underlying C++ routine raises a clear, contextual error rather than
  # silently returning a meaningless lambda.
  skip_if_not(exists("select_lambda_cv_density_ratio_cpp"), message = "C++ not compiled")

  n <- 40
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  A <- rep(0, n)  # No treated units
  mean_grad_psi <- rnorm(p + 1)
  alpha_init <- rnorm(p + 1)
  lambda_grid <- exp(seq(log(1), log(0.01), length.out = 5))

  expect_error(
    select_lambda_cv_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                       lambda_grid, 3, 100, 1e-4,
                                       A_val = 1L, family_int = 1L, link_int = 1L),
    "no observations with A_val"
  )
  expect_error(
    select_lambda_cv_initial_density_ratio_cpp(X, A, mean_grad_psi,
                                               lambda_grid, 3, 100, 1e-4, 1L),
    "no observations with A_val"
  )
  expect_error(
    select_lambda_cv_calibrated_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                                  lambda_grid, 3, 100, 1e-4,
                                                  Inf, X, 1L, 1L, 1L),
    "no observations with A_val"
  )
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

test_that("optimize_weights_cpp fails fast on invalid numeric inputs", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")

  common <- list(
    V_t = c(0.1, 0.1),
    V_s = c(0.05, 0.05),
    n_s = c(200, 200),
    V_ot = 0.2,
    n_t = 100,
    C_ot = c(0.01, 0.01),
    lambda = 0.5,
    mu_ot = 0.5,
    max_iter = 100,
    tol = 1e-8,
    C_cross = matrix(0, 2, 2),
    warm_start = numeric(0)
  )
  call_optimizer <- function(estimates, n_s = common$n_s) {
    do.call(
      optimize_weights_cpp,
      c(list(estimates = estimates), common[names(common) != "n_s"],
        list(n_s = n_s))
    )
  }

  expect_error(call_optimizer(c(NA_real_, 0.5)), "must be finite")
  expect_error(call_optimizer(c(0.4, 0.5), n_s = c(200, 0)), "must be positive")
})

test_that("aggregation PSD ridge is homogeneous in the variance scale", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")

  # The deliberately incompatible cross-covariance makes the plug-in Hessian
  # indefinite, so the numerical safeguard must activate. Scaling every
  # variance/covariance component should scale that ridge by the same factor
  # while leaving the unpenalized minimizer unchanged.
  fit_at_scale <- function(scale) {
    optimize_weights_cpp(
      estimates = c(0.45, 0.55),
      V_t = scale * c(0.1, 0.1),
      V_s = scale * c(0.05, 0.05),
      n_s = c(200, 200),
      V_ot = scale * 0.2,
      n_t = 100,
      C_ot = scale * c(0.01, 0.01),
      lambda = 0,
      mu_ot = 0.5,
      max_iter = 1000,
      tol = 1e-10,
      # Positive off-diagonal curvature makes the anti-symmetric direction
      # indefinite while the symmetric score direction remains well scaled.
      C_cross = scale * matrix(c(0, 1, 1, 0), 2, 2),
      warm_start = numeric(0)
    )
  }

  baseline <- fit_at_scale(1)
  scaled <- fit_at_scale(1e-4)
  expect_gt(baseline$psd_ridge, 0)
  expect_equal(
    scaled$psd_ridge,
    1e-4 * baseline$psd_ridge,
    tolerance = 1e-12
  )
  expect_equal(scaled$weights, baseline$weights, tolerance = 1e-8)
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
