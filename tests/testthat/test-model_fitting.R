# test-model_fitting.R - Tests for R-level model fitting wrappers
#
# Covers: fit_initial_outcome, fit_initial_density_ratio, fit_unified_density_ratio,
#         fit_unified_outcome, optimize_weights, calculate_aggregated_variance

library(testthat)

# Package loaded by helper-load.R (all functions available via RoCE namespace)

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

  selected <- .select_nuisance_cv_lambda(
    list(
      lambda_min = 0.1, lambda_1se = 0.2,
      invalid_fold_fits = 3L, invalid_lambdas = 1L,
      path_tail_skipped_fold_fits = 2L
    ),
    lambda_rule = "min",
    caller = "unit_test"
  )
  expect_equal(attr(selected, "cv_invalid_fold_fits"), 3L)
  expect_equal(attr(selected, "cv_invalid_lambdas"), 1L)
  expect_equal(attr(selected, "cv_path_tail_skipped_fold_fits"), 2L)
})

test_that("lambda-grid size is validated before nuisance fitting", {
  expect_error(
    build_lambda_grid(lambda_max = 1, nlambda = 1),
    "integer >= 2"
  )
  expect_length(build_lambda_grid(lambda_max = 1, nlambda = 7), 7L)
})

test_that("glmnet prediction helper validates its lambda rule", {
  expect_error(
    fit_glmnet_cv(
      matrix(rnorm(80), nrow = 20L),
      rep(c(0, 1), each = 10L),
      matrix(rnorm(20), nrow = 5L),
      lambda_rule = "unsupported"
    ),
    "lambda_rule must be either 'min' or '1se'"
  )
})

test_that("glmnet prediction helper reports its constant-outcome fallback", {
  prediction <- fit_glmnet_cv(
    x_train = matrix(seq_len(60), nrow = 20L, ncol = 3L),
    y_train = rep(1, 20L),
    x_predict = matrix(0, nrow = 4L, ncol = 3L),
    family = "binomial",
    on_degenerate_response = "constant"
  )
  expect_equal(as.numeric(prediction), rep(1, 4L))
  expect_identical(attr(prediction, "outcome_degenerate"), 1L)
})

test_that("glmnet prediction helper applies the requested lambda rule", {
  set.seed(9137)
  x_train <- matrix(rnorm(600), nrow = 120L, ncol = 5L)
  y_train <- 0.8 * x_train[, 1L] - 0.4 * x_train[, 2L] + rnorm(120L)
  x_predict <- matrix(rnorm(50), nrow = 10L, ncol = 5L)
  nlambda <- 7L
  nfolds <- get_cv_fold_count(nrow(x_train), min_per_fold = 10L)

  # Reset to the same RNG state before each fit so cv.glmnet constructs the
  # same folds. This makes the test about lambda selection, not fold noise.
  set.seed(20260814)
  reference <- glmnet::cv.glmnet(
    x = x_train, y = y_train, family = "gaussian",
    alpha = 1, nfolds = nfolds, nlambda = nlambda
  )
  expected_one_se <- as.numeric(stats::predict(
    reference, newx = x_predict, s = "lambda.1se", type = "response"
  ))
  expected_min <- as.numeric(stats::predict(
    reference, newx = x_predict, s = "lambda.min", type = "response"
  ))

  set.seed(20260814)
  observed_one_se <- fit_glmnet_cv(
    x_train, y_train, x_predict,
    family = "gaussian", nlambda = nlambda, lambda_rule = "1se"
  )
  set.seed(20260814)
  observed_min <- fit_glmnet_cv(
    x_train, y_train, x_predict,
    family = "gaussian", nlambda = nlambda, lambda_rule = "min"
  )

  expect_equal(observed_one_se, expected_one_se, tolerance = 1e-12)
  expect_equal(observed_min, expected_min, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(
    observed_one_se, observed_min, tolerance = 1e-12
  )))
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
  expect_true(is.logical(attr(gamma, "converged")))
  expect_true(is.finite(attr(gamma, "iterations")))
  expect_true(is.finite(attr(gamma, "max_update")))
  expect_true(is.finite(attr(gamma, "convergence_threshold")))
  expect_true(is.finite(attr(gamma, "update_to_threshold_ratio")))
  expect_true(is.finite(attr(gamma, "max_abs_coefficient")))
  expect_identical(attr(gamma, "line_search_failures"), 0L)
  if (isTRUE(attr(gamma, "converged"))) {
    expect_lte(attr(gamma, "update_to_threshold_ratio"), 1 + 1e-10)
  }
  expect_gte(attr(gamma, "cv_seconds"), 0)
  expect_gte(attr(gamma, "final_fit_seconds"), 0)
})

test_that("density-ratio support floor controls source-arm null directions", {
  skip_if_not(exists("fit_initial_density_ratio_cpp"), message = "C++ not compiled")

  set.seed(481)
  n <- 120L
  A <- rep(c(0, 1), each = n / 2L)
  Z_site <- cbind(
    source_arm_constant = c(rep(0, n / 2L), rnorm(n / 2L)),
    varying = rnorm(n)
  )
  mean_phi <- c(1, 0.4, 0)

  support <- .density_ratio_support_penalty_floor(
    Z_site, A, mean_phi, A_val = 0L, caller = "unit_test"
  )
  expect_gt(support$floor, 0.4)
  expect_equal(support$constant_feature_count, 1L)

  gamma <- fit_initial_density_ratio(
    Z_site, A, mean_phi, lambda = 0.01, A_val = 0L,
    max_iter = 2000L
  )
  expect_true(isTRUE(attr(gamma, "support_penalty_floor_applied")))
  expect_equal(attr(gamma, "lambda_selected_on_path"), 0.01)
  expect_equal(attr(gamma, "lambda_selected_before_support_floor"), 0.01)
  expect_false(isTRUE(attr(gamma, "support_path_floor_applied")))
  expect_equal(attr(gamma, "support_constant_feature_count"), 1L)
  expect_gt(attr(gamma, "lambda_used"), 0.4)
  expect_true(all(is.finite(gamma)))
  expect_lt(max(abs(gamma)), 0.9 * PARAM_MAX)
})

test_that("density-ratio support floor is inactive without constant features", {
  set.seed(482)
  Z_site <- matrix(rnorm(160), 80L, 2L)
  A <- rep(c(0, 1), each = 40L)
  support <- .density_ratio_support_penalty_floor(
    Z_site, A, c(1, 0, 0), A_val = 0L, caller = "unit_test"
  )

  expect_equal(support$floor, 0)
  expect_equal(support$constant_feature_count, 0L)
})

test_that("density-ratio CV grid records support-floor constraints", {
  constrained <- .constrain_density_ratio_lambda_grid(
    c(1, 0.1, 0.01), list(floor = 0.2)
  )
  untouched <- .constrain_density_ratio_lambda_grid(
    c(1, 0.5, 0.3), list(floor = 0.2)
  )

  expect_equal(as.numeric(constrained), c(1, 0.2, 0.2))
  expect_true(isTRUE(attr(constrained, "support_path_floor_applied")))
  expect_equal(as.numeric(untouched), c(1, 0.5, 0.3))
  expect_false(isTRUE(attr(untouched, "support_path_floor_applied")))

  applied <- .apply_density_ratio_support_floor(
    lambda = 0.2,
    Z_site = matrix(c(0, 0, 1, 2), ncol = 1L),
    A = c(0, 0, 1, 1),
    linear_moment = c(1, 0.2),
    A_val = 0L,
    caller = "unit_test",
    support = list(floor = 0.2, constant_feature_count = 1L),
    path_floor_applied = TRUE
  )
  expect_equal(attr(applied, "lambda_selected_on_path"), 0.2)
  expect_true(is.na(attr(applied, "lambda_selected_before_support_floor")))
  expect_true(isTRUE(attr(applied, "support_path_floor_applied")))
})

test_that("density-ratio warm starts reject unsafe previous fits", {
  make_fit <- function(values, converged = TRUE, failures = 0L) {
    attr(values, "converged") <- converged
    attr(values, "line_search_failures") <- failures
    values
  }

  safe <- make_fit(c(0.2, -0.3))
  expect_equal(.safe_density_ratio_warm_start(safe), as.numeric(safe))
  expect_null(.safe_density_ratio_warm_start(
    make_fit(c(0.2, -0.3), converged = FALSE)
  ))
  expect_null(.safe_density_ratio_warm_start(
    make_fit(c(0.2, -0.3), failures = 1L)
  ))
  expect_null(.safe_density_ratio_warm_start(
    make_fit(c(0, 0.95 * PARAM_MAX))
  ))
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

test_that("optimize_weights rejects non-finite inputs and invalid sample sizes", {
  estimates <- c(0.50, 0.62)
  variances <- list(V_ot = 1.0, V_t = c(0.8, 0.9), V_s = c(0.7, 0.6))
  C_ot <- c(0.05, 0.04)
  n_samples <- list(n_t = 40, n_s = c(45, 50))

  expect_error(
    optimize_weights(
      c(0.50, NA_real_), variances, C_ot, n_samples,
      lambda = 0.5, mu_ot = 0.50
    ),
    "must be finite"
  )
  expect_error(
    optimize_weights(
      estimates, variances, C_ot, list(n_t = 40, n_s = c(45, 0)),
      lambda = 0.5, mu_ot = 0.50
    ),
    "must be positive"
  )
  expect_error(
    optimize_weights(
      estimates, variances, C_ot, n_samples,
      lambda = 0.5, mu_ot = 0.50,
      C_cross = matrix(c(0, NA, NA, 0), 2, 2)
    ),
    "C_cross must contain only finite"
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
  expect_gte(attr(gamma, "cv_seconds"), 0)
  expect_gte(attr(gamma, "final_fit_seconds"), 0)
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

# Ensure the shared density-ratio solver handles strongly heterogeneous scales.
test_that("fit_unified_density_ratio_cpp converges with backtracking", {
  skip_if_not(exists("fit_unified_density_ratio_cpp"), message = "C++ not compiled")

  set.seed(2026)
  n <- 300
  p <- 6
  X <- sweep(matrix(rnorm(n * p), n, p), 2, c(0.01, 0.1, 1, 10, 100, 1000), `*`)
  X <- scale(X)
  A <- rbinom(n, 1, 0.45)
  alpha_init <- rnorm(p + 1, sd = 1.0)
  eta <- drop(cbind(1, X) %*% alpha_init)
  psi_prime <- stats::plogis(eta) * (1 - stats::plogis(eta))
  mean_grad_psi <- colMeans(cbind(1, X) * psi_prime)

  res <- fit_unified_density_ratio_cpp(X, A, mean_grad_psi, alpha_init,
                                       0.01, 1000, 1e-6, TRUE, 10.0,
                                       W_outcome = X,
                                       A_val = 1L, family_int = 1L, link_int = 1L,
                                       warm_start = rep(0, p + 1))

  expect_true(is.list(res))
  expect_equal(length(res$gamma), p + 1)
  expect_true(all(is.finite(res$gamma)))
  expect_true(res$converged)
  expect_lte(res$max_update, res$convergence_threshold * (1 + 1e-10))
  expect_identical(res$line_search_failures, 0L)
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

test_that("correction diagnostics distinguish configured logit truncation", {
  skip_if_not(exists("calculate_correction_term_cpp"), message = "C++ not compiled")

  n <- 12
  Z <- matrix(c(rep(-30, n), rep(0, n)), n, 2)
  W <- matrix(rnorm(n * 2), n, 2)
  A <- rep(c(1, 0), length.out = n)
  res <- calculate_correction_term_cpp(
    Z, A, rnorm(n), c(0, 1, 0), c(0, 0, 0), W,
    M_tau = 5,
    family_int = FAMILY_GAUSSIAN,
    link_int = LINK_IDENTITY,
    A_val = 1L
  )

  expect_equal(res$clip_diagnostics$n_obs, sum(A == 1))
  expect_equal(res$clip_diagnostics$logit_truncated, sum(A == 1))
  expect_equal(res$clip_diagnostics$logit_truncation_fraction, 1)
  expect_true(isTRUE(res$clip_diagnostics$any_truncated))
  expect_false(isTRUE(res$clip_diagnostics$any_safety_clipped))
})

# ============================================================================
# Tests for fit_unified_outcome
# ============================================================================

expect_valid_outcome_fit_diagnostics <- function(alpha) {
  expect_true(is.logical(attr(alpha, "converged")))
  expect_true(is.finite(attr(alpha, "iterations")))
  expect_true(is.finite(attr(alpha, "max_update")))
  expect_true(is.finite(attr(alpha, "convergence_threshold")))
  expect_true(is.finite(attr(alpha, "update_to_threshold_ratio")))
  expect_true(is.finite(attr(alpha, "max_abs_coefficient")))
  expect_identical(attr(alpha, "line_search_failures"), 0L)
  if (isTRUE(attr(alpha, "converged"))) {
    expect_lte(attr(alpha, "update_to_threshold_ratio"), 1 + 1e-10)
  }
}

test_that("fit_unified_outcome returns valid alpha (refined)", {
  skip_if_not(exists("fit_unified_outcome_cpp"), message = "C++ not compiled")
  
  d <- make_site_data(200, 4)
  gamma_s <- rnorm(ncol(d$X) + 1, sd = 0.1)
  
  alpha <- fit_unified_outcome(d$X, d$Y, d$A, A_val = 1, gamma_s = gamma_s,
                                lambda = 0.05, calibrated = FALSE, Z_site = d$X)
  
  expect_equal(length(alpha), ncol(d$X) + 1)
  expect_true(all(is.finite(alpha)))
  expect_valid_outcome_fit_diagnostics(alpha)
})

test_that("fit_unified_outcome returns valid alpha (calibrated)", {
  skip_if_not(exists("fit_unified_outcome_cpp"), message = "C++ not compiled")
  
  d <- make_site_data(200, 4)
  gamma_s <- rnorm(ncol(d$X) + 1, sd = 0.1)
  
  alpha <- fit_unified_outcome(d$X, d$Y, d$A, A_val = 1, gamma_s = gamma_s,
                                lambda = 0.05, calibrated = TRUE, M_tau = 10.0, Z_site = d$X)
  
  expect_equal(length(alpha), ncol(d$X) + 1)
  expect_true(all(is.finite(alpha)))
  expect_valid_outcome_fit_diagnostics(alpha)
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
  expect_true(is.finite(attr(weights, "optimizer_iterations")))
  expect_gte(attr(weights, "optimizer_iterations"), 1L)
  expect_true(is.finite(attr(weights, "psd_ridge")))
  expect_gte(attr(weights, "psd_ridge"), 0)
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

test_that("optimize_weights uses the N_all-scaled main.tex objective", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")

  estimates <- 0.35
  mu_ot <- 0.50
  V_ot <- 0.20
  V_t <- 0.10
  V_s <- 0.05
  C_ot <- 0.01
  n_t <- 100
  n_s <- 200
  lambda <- 0.50

  var_disc <- V_ot / n_t + V_t / n_t + V_s / n_s -
    2 * C_ot / n_t
  t_stat <- abs(mu_ot - estimates) / sqrt(var_disc)
  penalty <- max(lambda * t_stat - 1, 0)
  eta_unpenalized <- ((V_ot - C_ot) / n_t) / var_disc
  expected <- max(
    eta_unpenalized - penalty / ((n_t + n_s) * 2 * var_disc),
    0
  )

  result <- optimize_weights(
    estimates = estimates,
    variances = list(V_t = V_t, V_s = V_s, V_ot = V_ot),
    C_ot = C_ot,
    n_samples = list(n_t = n_t, n_s = n_s),
    lambda = lambda,
    mu_ot = mu_ot,
    C_cross = matrix(0, 1, 1),
    clip_weights = FALSE
  )

  expect_gt(penalty, 0)
  expect_equal(as.numeric(result), expected, tolerance = 1e-10)
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

test_that("unconstrained weight optimization does not silently zero large solutions", {
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ not compiled")

  result <- optimize_weights(
    estimates = 0.5,
    variances = list(V_ot = 1, V_t = 0.9801, V_s = 1e-6),
    C_ot = 0.99,
    n_samples = list(n_t = 100, n_s = 100),
    lambda = 0,
    mu_ot = 0.5,
    C_cross = matrix(0, 1, 1),
    clip_weights = FALSE
  )

  expect_gt(as.numeric(result), 10)
  expect_equal(
    as.numeric(result),
    ((1 - 0.99) / 100) /
      (1 / 100 + 0.9801 / 100 + 1e-6 / 100 - 2 * 0.99 / 100),
    tolerance = 1e-8
  )
})

test_that("source nuisance iteration budget is validated before fitting", {
  expect_error(
    process_source_site(
      s = "s1",
      source_folds = NULL,
      target_folds = NULL,
      k1 = 1L,
      n_folds = 2L,
      A_val = 1L,
      M_tau = M_TAU_DEFAULT,
      data_split = NULL,
      get_fold_inputs = function(...) NULL,
      nuisance_max_iter = 0L
    ),
    "nuisance_max_iter must be one positive integer",
    fixed = TRUE
  )
})

test_that("run_crossfit forwards and audits the source nuisance grid size", {
  run_crossfit_body <- paste(deparse(body(run_crossfit)), collapse = "\n")
  expect_match(
    run_crossfit_body,
    "nuisance_nlambda = nlambda_init",
    fixed = TRUE
  )
  expect_match(
    run_crossfit_body,
    "source nuisance grid mismatch",
    fixed = TRUE
  )

  expect_error(
    process_source_site(
      s = "s1",
      source_folds = NULL,
      target_folds = NULL,
      k1 = 1L,
      n_folds = 2L,
      A_val = 1L,
      M_tau = M_TAU_DEFAULT,
      data_split = NULL,
      get_fold_inputs = function(...) NULL,
      nuisance_nlambda = 1L
    ),
    "nuisance_nlambda must be one integer >= 2",
    fixed = TRUE
  )
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

test_that("prepare_cross_matrix rejects malformed non-null covariance inputs", {
  expect_equal(dim(prepare_cross_matrix(NULL, 2L)), c(0L, 0L))
  expect_error(
    prepare_cross_matrix(c(0, 1, 1, 0), 2L),
    "NULL or a numeric matrix"
  )
  expect_error(
    prepare_cross_matrix(matrix(0, 2L, 3L), 2L),
    "dimension 2 x 2"
  )
  expect_error(
    prepare_cross_matrix(matrix(c(0, NA, NA, 0), 2L, 2L), 2L),
    "only finite values"
  )
})
