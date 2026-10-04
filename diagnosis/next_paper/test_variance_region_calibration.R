#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
source(file.path(arguments[[1L]], "source_confidence_sets.R"))
source(file.path(arguments[[1L]], "variance_region_calibration.R"))
library(testthat)

test_that("conditional envelope is below the exact loading probabilities", {
  for (common_noise in c(0, 1, 2, 4)) {
    for (critical in c(.5, 2)) {
      lower <- loading_region_probabilities(critical, common_noise, .2, .95)
      loading <- seq(.2, .95, length.out = 101)
      exact <- normal_interval_probability(critical, loading * common_noise, sqrt(1 - loading^2))
      expect_lte(lower, min(exact) + 1e-14)
    }
  }
  expect_equal(normal_interval_probability(1, c(0, 1, 2), c(0, 0, 0)), c(1, 1, 0))
})

test_that("point variance regions recover the known-variance calculation", {
  private <- c(.2, .5, 1.2, 2)
  for (shared in c(0, .6)) {
    expect_equal(variance_region_vote_cdf(2, 3, 2, rep(shared, 2), cbind(private, private)),
      shared_target_vote_cdf(2, 3, 2, shared, private), tolerance = 1e-8)
    region <- make_source_region_calibration(3, 2, rep(shared, 2), cbind(private, private),
      1, variance_alpha = 0)
    known <- make_source_calibration(4, 3, 2, method = "shared_target_gaussian",
      shared_variance = shared, source_private_variances = private)
    expect_equal(region$source_critical, known$source_critical, tolerance = 2e-7)
  }
})

test_that("the envelope protects every fixed valid subset inside the region", {
  shared_bounds <- c(.3, .8)
  private_bounds <- cbind(c(.1, .3, .7, 1), c(.3, .8, 1.5, 3))
  lower <- variance_region_vote_cdf(2, 3, 2, shared_bounds, private_bounds)
  for (shared in c(.3, .55, .8)) {
    for (fraction in c(0, .5, 1)) {
      private <- (1 - fraction) * private_bounds[, 1] + fraction * private_bounds[, 2]
      exact <- vapply(combn(4, 3, simplify = FALSE), function(index)
        shared_target_vote_cdf(2, 3, 2, shared, private[index]), numeric(1))
      expect_lte(lower, min(exact) + 1e-8)
    }
  }
})

test_that("calibration respects its probability budget and scale equivariance", {
  private_bounds <- cbind(c(.1, .3, .7, 1), c(.3, .8, 1.5, 3))
  calibration <- make_source_region_calibration(3, 2, c(.3, .8), private_bounds, 2)
  scaled <- make_source_region_calibration(3, 2, 4 * c(.3, .8), 4 * private_bounds, 8)
  expect_equal(calibration$anchor_alpha + calibration$source_alpha + calibration$variance_alpha, .05)
  expect_equal(calibration$source_critical, scaled$source_critical, tolerance = 1e-7)
  reference <- make_source_calibration(4, 3, 2, alpha = .045)
  expect_lte(calibration$source_critical, reference$source_critical)
  result <- source_region_interval(0, c(-.3, .2, .4, 1), calibration)
  transformed <- source_region_interval(3, 3 + 2 * c(-.3, .2, .4, 1), scaled)
  expect_equal(transformed$intervals, 3 + 2 * result$intervals, tolerance = 1e-7)
})

test_that("variance bounds have exact normal-sample marginal tail probabilities", {
  variances <- c(.3, 1, 4)
  degrees_freedom <- c(9, 99, 999)
  bounds <- gaussian_variance_bounds(variances, degrees_freedom, .006)
  expect_equal(pchisq(degrees_freedom * variances / bounds[, 1], degrees_freedom,
                     lower.tail = FALSE), rep(.001, 3), tolerance = 1e-12)
  expect_equal(pchisq(degrees_freedom * variances / bounds[, 2], degrees_freedom),
               rep(.001, 3), tolerance = 1e-12)
  expect_equal(gaussian_variance_bounds(4 * variances, 99, .006),
               4 * gaussian_variance_bounds(variances, 99, .006))
  expect_equal(unname(gaussian_variance_bounds(0, 99, .006)), matrix(c(0, 0), 1))
  expect_error(gaussian_variance_bounds(variances, c(9, 99), .006), "Invalid")
})

test_that("wide and degenerate regions remain defined with explicit invalid-input checks", {
  wide <- make_source_region_calibration(3, 2, c(0, 10), cbind(rep(0, 4), rep(10, 4)), 20)
  expect_identical(wide$calibration_bound, "marginal_bound")
  expect_true(is.finite(wide$source_critical))
  expect_equal(variance_region_vote_cdf(2, 3, 2, c(0, 0), cbind(rep(0, 4), rep(1, 4))),
               equicorrelated_vote_cdf(2, 3, 2, 0))
  expect_error(make_source_region_calibration(3, 2, c(1, .5), cbind(1:4, 2:5), 1), "shared_bounds")
  expect_error(make_source_region_calibration(3, 2, c(0, 1), cbind(1:4, 0:3), 1), "private_bounds")
  expect_error(make_source_region_calibration(3, 2, c(0, 1), cbind(1:4, 2:5), 1,
                                             variance_alpha = .05), "budget")
})

test_that("the recorded heterogeneous pipeline case converges without relaxing tolerance", {
  # The third draw in cell 2 previously failed at a kink during root finding.
  set.seed(91002)
  components <- c(.5, 1, exp(seq(log(.25), log(4), length.out = 16))) / 100
  for (iteration in 1:3) {
    rnorm(1, sd = sqrt(components[1]))
    rnorm(17, sd = sqrt(components[-1]))
    estimates <- components * rchisq(18, 99) / 99
  }
  bounds <- gaussian_variance_bounds(estimates, 99, .005)
  calibration <- make_source_region_calibration(12, 9, bounds[1, ],
    bounds[-c(1, 2), ], sum(bounds[1:2, "upper"]))
  expect_null(calibration$integration_fallback_reason)
  expect_identical(calibration$calibration_bound, "conditional_variance_region")
  expect_equal(variance_region_vote_cdf(calibration$source_critical, 12, 9,
    bounds[1, ], bounds[-c(1, 2), ]), 1 - calibration$source_alpha, tolerance = 2e-6)
})

cat("VARIANCE_REGION_TESTS_PASSED\n")
