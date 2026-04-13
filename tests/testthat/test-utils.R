# test-utils.R - Tests for utility functions

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEC namespace)

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

