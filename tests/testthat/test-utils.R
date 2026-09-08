# test-utils.R - Tests for utility functions

library(testthat)

# Package loaded by helper-load.R (all functions available via RoCE namespace)

test_that("scale_center works on vectors", {
  x <- c(1, 2, 3, 4, 5)
  scaled <- scale_center(x)
  
  expect_true(abs(mean(scaled)) < 1e-10)
  expect_true(abs(sd(scaled) - 1) < 1e-10)
})

test_that("scale_center works on matrices", {
  X <- matrix(rnorm(100), nrow = 20, ncol = 5)
  scaled <- scale_center(X)
  
  col_means <- colMeans(scaled)
  col_sds <- apply(scaled, 2, sd)
  
  expect_true(all(abs(col_means) < 1e-10))
  expect_true(all(abs(col_sds - 1) < 1e-10))
})

test_that("scale_center handles zero variance (constant values)", {
  # Vector with zero variance
  x_const <- rep(5, 10)
  scaled_const <- scale_center(x_const)
  
  # Should only center (all zeros), not throw error
  expect_true(all(scaled_const == 0))
  
  # Matrix with one constant column
  X_mixed <- cbind(1:5, rep(3, 5))
  scaled_mixed <- scale_center(X_mixed)
  
  # First column should be scaled, second column should be centered only
  expect_true(abs(mean(scaled_mixed[, 1])) < 1e-10)
  expect_true(all(scaled_mixed[, 2] == 0))
})

test_that("clip_to_range respects bounds", {
  x <- c(-5, -2, 0, 2, 5)
  clipped <- clip_to_range(x, lower = -2.5, upper = 2.5)
  
  expect_true(all(clipped >= -2.5))
  expect_true(all(clipped <= 2.5))
  expect_equal(clipped[3], 0)  # Unchanged
})

test_that("logistic function is bounded", {
  x <- c(-100, -1, 0, 1, 100)
  result <- logistic(x)
  
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
  expect_equal(result[3], 0.5, tolerance = 1e-10)  # logistic(0) = 0.5
  # For moderate inputs, result is strictly between 0 and 1
  expect_true(logistic(5) > 0 && logistic(5) < 1)
  expect_true(logistic(-5) > 0 && logistic(-5) < 1)
})

test_that("normalize_to_unit produces unit norm", {
  x <- c(3, 4)  # Should become (0.6, 0.8)
  result <- normalize_to_unit(x)
  
  expect_equal(sqrt(sum(result^2)), 1, tolerance = 1e-10)
  expect_equal(result, c(0.6, 0.8), tolerance = 1e-10)
})

test_that("normalize_to_unit handles zero vector", {
  x <- c(0, 0, 0)
  result <- normalize_to_unit(x)
  
  expect_equal(result, x)  # Zero vector unchanged
})

test_that("empirical_expectation computes correct mean", {
  x <- c(1, 2, 3, 4, 5)
  result <- empirical_expectation(x)
  
  expect_equal(result, 3)
})

test_that("empirical_expectation handles weights", {
  x <- c(1, 2, 3)
  weights <- c(1, 2, 1)  # More weight on x=2
  result <- empirical_expectation(x, weights)
  
  expected <- sum(x * weights) / sum(weights)  # (1 + 4 + 3) / 4 = 2
  expect_equal(result, expected)
})

test_that("safe_log handles edge cases", {
  x <- c(0, 0.5, 1, 2)
  result <- safe_log(x, epsilon = 1e-10)
  
  # First element should not be -Inf
  expect_true(is.finite(result[1]))
  expect_equal(result[2], log(0.5), tolerance = 1e-10)
})

test_that("safe_var handles single values and NAs", {
  x <- c(5)
  result <- safe_var(x)
  expect_equal(result, 1e-8)  # Returns min_var for single values
  
  x_with_na <- c(1, 2, NA, 4)
  result <- safe_var(x_with_na)
  expect_true(is.finite(result))
})

test_that("build_lambda_grid returns a glmnet-style descending path", {
  grid <- build_lambda_grid(lambda_max = 2, lambda_min_ratio = 1e-3, nlambda = 6)

  expect_equal(grid[1], 2, tolerance = 1e-12)
  expect_equal(grid[length(grid)], 2e-3, tolerance = 1e-12)
  expect_true(all(diff(grid) < 0))
})

test_that("lambda_max helpers fail fast on empty treatment arm", {
  Z <- matrix(rnorm(20), nrow = 10, ncol = 2)
  W <- matrix(rnorm(20), nrow = 10, ncol = 2)
  A <- rep(0L, 10)
  Y <- rnorm(10)
  mean_grad_psi <- c(1, 0, 0)
  alpha_init <- c(0, 0, 0)
  gamma_s <- c(0, 0, 0)

  expect_error(
    compute_lambda_max_refined_dr(Z, A, mean_grad_psi, alpha_init, A_val = 1L),
    "no observations with A_val"
  )
  expect_error(
    compute_lambda_max_outcome(W, Y, A, gamma_s, A_val = 1L, Z_site = Z),
    "no observations with A_val"
  )
})

test_that("validate_algorithm_inputs catches invalid data", {
  # Create minimal valid data
  valid_data <- list(
    t = list(
      W_outcome = matrix(rnorm(20), nrow = 10, ncol = 2),
      Z_site = matrix(rnorm(20), nrow = 10, ncol = 2),
      A = sample(c(0, 1), 10, replace = TRUE),
      Y = sample(c(0, 1), 10, replace = TRUE),
      n = 10
    ),
    s1 = list(
      W_outcome = matrix(rnorm(20), nrow = 10, ncol = 2),
      Z_site = matrix(rnorm(20), nrow = 10, ncol = 2),
      A = sample(c(0, 1), 10, replace = TRUE),
      Y = sample(c(0, 1), 10, replace = TRUE),
      n = 10
    )
  )
  
  # Should not throw error for valid data
  expect_silent(validate_algorithm_inputs(valid_data, "cv"))

  # Test with non-binary Y (should give warning or error)
  invalid_data <- valid_data
  invalid_data$t$Y <- rnorm(10)  # Continuous, not binary

  expect_error(validate_algorithm_inputs(invalid_data, "cv"))
})

test_that("validate_algorithm_inputs accepts lambda_selection = NULL", {
  # Comparison methods (estimate_tilted_aipw, estimate_federated_dr,
  # estimate_pooled_dr) do not use a lambda parameter; they pass
  # lambda_selection = NULL so the validator skips that check.
  valid_data <- list(
    t = list(
      W_outcome = matrix(rnorm(20), nrow = 10, ncol = 2),
      Z_site = matrix(rnorm(20), nrow = 10, ncol = 2),
      A = rep(c(0, 1), 5),
      Y = rep(c(0, 1), 5),
      n = 10
    ),
    s1 = list(
      W_outcome = matrix(rnorm(20), nrow = 10, ncol = 2),
      Z_site = matrix(rnorm(20), nrow = 10, ncol = 2),
      A = rep(c(0, 1), 5),
      Y = rep(c(0, 1), 5),
      n = 10
    )
  )

  expect_silent(validate_algorithm_inputs(valid_data))
  expect_silent(validate_algorithm_inputs(valid_data, lambda_selection = NULL))
  expect_silent(validate_algorithm_inputs(valid_data, lambda_selection = 2))
  expect_silent(validate_algorithm_inputs(valid_data, family = "binomial"))
  expect_error(validate_algorithm_inputs(valid_data, lambda_selection = LAMBDA_MAX * 2),
               "finite numeric scalar")
  expect_error(validate_algorithm_inputs(valid_data, lambda_selection = c("cv", "bad")),
               "must be 'cv'")
})

test_that("comparison estimators validate inputs at entry", {
  # Regression test: every exported comparison estimator must reject malformed
  # data_split BEFORE attempting any nuisance fitting. The shared entry
  # validator (validate_algorithm_inputs) is responsible for these checks.
  bad_inputs <- list(
    "NULL"          = NULL,
    "empty list"    = list(),
    "no target 't'" = list(s1 = list())
  )

  comparison_estimators <- list(
    estimate_tilted_aipw           = estimate_tilted_aipw,
    estimate_federated_dr          = estimate_federated_dr,
    estimate_pooled_dr             = estimate_pooled_dr,
    estimate_sample_size_weighted  = estimate_sample_size_weighted,
    estimate_inverse_variance_weighted = estimate_inverse_variance_weighted
  )

  for (est_name in names(comparison_estimators)) {
    est_fn <- comparison_estimators[[est_name]]
    for (label in names(bad_inputs)) {
      bad <- bad_inputs[[label]]
      expect_error(est_fn(bad),
                   info = paste(est_name, "must reject:", label))
    }
  }
})

test_that("RHC imputers fail fast when no non-NA values remain", {
  # Regression test: previously these helpers silently returned all-NA or
  # unchanged data, deferring failure to a far less informative call site.
  na_vec <- rep(NA_real_, 5L)
  expect_error(.rhc_impute_continuous(na_vec, var_name = "test_cont"),
               "no non-NA values")
  expect_error(.rhc_impute_mode(rep(NA_character_, 5L), var_name = "test_mode"),
               "no non-NA values")

  # And succeed on the happy path so the strict check does not over-fire.
  cont <- c(1.0, 2.0, NA, 4.0, 5.0)
  imputed <- .rhc_impute_continuous(cont, var_name = "happy")
  expect_equal(sum(is.na(imputed)), 0L)
  expect_equal(imputed[3L], stats::median(cont, na.rm = TRUE))

  cat <- c("a", "b", NA, "a", "")
  imputed_mode <- .rhc_impute_mode(cat, var_name = "happy_mode")
  expect_equal(sum(is.na(imputed_mode)), 0L)
  expect_true(all(imputed_mode != ""))
  expect_equal(imputed_mode[3L], "a")  # mode of {"a","b","a"}
})

test_that("RHC continuous imputer rejects non-numeric tokens explicitly", {
  expect_error(
    .rhc_impute_continuous(c("1.2", "not_numeric", NA), var_name = "bad_cont"),
    "non-numeric non-missing"
  )
})

test_that("solve_with_ridge fails instead of switching solvers", {
  singular_mat <- matrix(c(1, 1, 1, 1), nrow = 2L)

  expect_error(
    solve_with_ridge(singular_mat, ridge = 0),
    "refusing to switch"
  )
})

test_that("calculate_weighted_site_aipw warns when too few treated units", {
  # Regression test: previously this function silently returned NA / Inf with
  # a zero influence function when the treated-arm cell was below
  # MIN_TREATED_FOR_MODEL. It must now emit a warning carrying the counts.
  set.seed(1)
  n <- 50
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  Y <- rbinom(n, 1, 0.5)
  A <- rep(0L, n)              # zero treated units
  A[1L] <- 1L                  # one treated unit -- below MIN_TREATED_FOR_MODEL

  expect_warning(
    res <- calculate_weighted_site_aipw(y = Y, a = A, X = X, A_val = 1L),
    "calculate_weighted_site_aipw: only 1 unit"
  )
  expect_true(is.na(res$estimate))
  expect_identical(typeof(res$estimate), "double")  # NA_real_, not NA (logical)
  expect_true(is.infinite(res$variance))
  expect_equal(res$psi, rep(0, n))
})
