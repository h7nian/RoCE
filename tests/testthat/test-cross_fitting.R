# test-cross_fitting.R - Integration tests for cross-fitting algorithms
#
# Covers: run_crossfit with communication_mode = "two_round" / "one_round"
# Uses small synthetic datasets to verify the full pipeline executes
# and produces outputs with the expected structure and properties.

library(testthat)

# Package loaded by helper-load.R (all functions available via RoCE namespace)

# make_small_data_split + get_smoke_data_split + get_smoke_result are shared
# helpers loaded via tests/testthat/helper-data.R before this file is sourced.

# ============================================================================
# Tests for run_crossfit(..., communication_mode = "two_round")
# ============================================================================

test_that("run_crossfit(two_round) returns expected structure", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  result <- get_smoke_result("two_round")
  
  # Check required output fields
  expect_true("estimate" %in% names(result))
  expect_true("variance" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("fold_lambdas" %in% names(result))
  expect_true("fold_lambda_info" %in% names(result))
  expect_true("clip_diagnostics" %in% names(result))
  expect_equal(result$aggregation_lambda_rule, "min")
  
  # Estimate should be finite
  expect_true(is.finite(result$estimate))
  
  # Variance should be non-negative
  expect_true(result$variance >= 0)
  
  # SE should be sqrt of variance
  expect_equal(result$se, sqrt(result$variance), tolerance = 1e-10)
  
  # CI should contain the estimate
  expect_true(result$ci_lower <= result$estimate)
  expect_true(result$ci_upper >= result$estimate)
})

test_that("run_crossfit(two_round) rejects too few folds", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")

  # get_smoke_data_split() generates continuous Y; pass family="gaussian" so
  # validate_algorithm_inputs does not raise its binary-Y check before the
  # n_folds validation we are exercising here.
  data_split <- get_smoke_data_split()

  expect_error(
    run_crossfit(data_split, n_folds = 2, communication_mode = "two_round",
                 family = "gaussian", verbose = FALSE),
    "too small for two-level cross-fitting"
  )
})

# ============================================================================
# Tests for run_crossfit(..., communication_mode = "one_round")
# ============================================================================

test_that("run_crossfit(one_round) returns expected structure", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  result <- get_smoke_result("one_round")
  
  # Check required output fields
  expect_true("estimate" %in% names(result))
  expect_true("variance" %in% names(result))
  expect_true("se" %in% names(result))
  expect_true("ci_lower" %in% names(result))
  expect_true("ci_upper" %in% names(result))
  expect_true("fold_lambdas" %in% names(result))
  expect_true("fold_lambda_info" %in% names(result))
  expect_true("clip_diagnostics" %in% names(result))
  expect_equal(result$aggregation_lambda_rule, "min")
  
  # Basic sanity checks
  expect_true(is.finite(result$estimate))
  expect_true(result$variance >= 0)
  expect_true(result$ci_lower <= result$estimate)
  expect_true(result$ci_upper >= result$estimate)
})

# ============================================================================
# Tests comparing one-round and two-round estimates
# ============================================================================

test_that("one-round and two-round produce finite but potentially different estimates", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  
  res_2r <- get_smoke_result("two_round")
  res_1r <- get_smoke_result("one_round")
  
  # Both should give valid estimates
  expect_true(is.finite(res_2r$estimate))
  expect_true(is.finite(res_1r$estimate))
  
  # Both should be finite
  expect_true(is.finite(res_2r$estimate))
  expect_true(is.finite(res_1r$estimate))
  
  # Both should have non-negative variance
  expect_true(res_2r$variance >= 0)
  expect_true(res_1r$variance >= 0)
})
