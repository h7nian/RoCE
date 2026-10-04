#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Usage: dispersion_source_file")
source(arguments[[1L]])
library(testthat)

test_that("vectorized noncentrality inversion agrees with probability inversion", {
  for (degrees in c(3, 15, 63)) for (probability in c(.01, .025)) {
    parameters <- c(0, 1, 5, 30)
    values <- qchisq(probability, degrees, ncp = parameters)
    upper <- noncentrality_upper_bound(values, degrees, probability)
    expect_true(all(upper + 1e-7 >= parameters))
    expect_lt(max(abs(upper - parameters)), 1e-6)
    actual <- pchisq(values, degrees, ncp = upper)
    expect_true(all(actual <= probability + 1e-12))
    for (index in 2:4) {
      reference <- uniroot(function(value) pchisq(values[index], degrees, ncp = value) - probability,
                            c(0, 200), tol = 1e-11)$root
      expect_equal(upper[index], reference, tolerance = 1e-7)
    }
  }
  expect_equal(noncentrality_upper_bound(0, 3, .025), 0)
  expect_error(noncentrality_upper_bound(-1, 3, .025), "nonnegative")
  expect_error(noncentrality_upper_bound(1, 0, .025), "integer df")
})

test_that("exact and inexpensive bounds match compound symmetry", {
  covariance <- matrix(.2, 6, 6) + diag(.7, 6)
  exact <- make_dispersion_calibration(covariance, 3, 1L, "subset_exact")
  upper <- make_dispersion_calibration(covariance, 3, 1L, "uniform_upper")
  expect_equal(exact$variance, .2 + .7 / 6, tolerance = 1e-14)
  expect_equal(exact$weights, rep(1/6, 6), tolerance = 1e-14)
  expect_equal(exact$bias_factor, .7 * (1/3 - 1/6), tolerance = 1e-14)
  expect_equal(upper$bias_factor, exact$bias_factor, tolerance = 1e-14)
  expect_equal(exact$subsets_examined, choose(5, 2))
  expect_equal(make_dispersion_calibration(covariance, 6)$bias_factor, 0)
})

test_that("a known valid anchor restricts the admissible subsets", {
  covariance <- diag(c(1, 3, 4, 5))
  restricted <- make_dispersion_calibration(covariance, 2, 1L)
  unrestricted <- make_dispersion_calibration(covariance, 2)
  upper <- make_dispersion_calibration(covariance, 2, 1L, "uniform_upper")
  expect_lt(restricted$bias_factor, unrestricted$bias_factor)
  expect_lte(restricted$bias_factor, upper$bias_factor + 1e-12)
  expect_lte(upper$bias_factor, covariance[1, 1] - restricted$variance + 1e-12)
  expect_error(make_dispersion_calibration(covariance, 1, c(1, 2)), "valid-candidate")
  expect_error(make_dispersion_calibration(diag(20), 10, max_subsets = 10), "enumeration")
})

test_that("two-arm noise variance retains covariance and error allocations", {
  covariance <- diag(6)
  covariance[1:3, 4:6] <- diag(.3, 3)
  covariance[4:6, 1:3] <- diag(.3, 3)
  calibration <- make_arm_dispersion_calibration(covariance, 2, 2)
  expect_equal(calibration$variance, 2 * (1 - .3) / 3, tolerance = 1e-14)
  values <- matrix(c(.6, .7, .5, .4, .3, .5), 1)
  interval <- arm_dispersion_intervals(values, calibration)
  radius <- qnorm(.975) * sqrt(calibration$variance)
  expect_equal(interval$pooled_estimate, .2, tolerance = 1e-14)
  expect_equal(c(interval$lower, interval$upper), .2 + c(-1, 1) * radius, tolerance = 1e-14)
  expect_equal(interval$bias_allowance, 0)
  no_sources <- make_arm_dispersion_calibration(covariance, 0, 0)
  anchor <- arm_dispersion_intervals(values, no_sources)
  expect_equal(c(anchor$lower, anchor$upper), .2 + c(-1, 1) * qnorm(.975) * sqrt(1.4), tolerance = 1e-14)
  partial <- make_arm_dispersion_calibration(covariance, 1, 2)
  bounded <- arm_dispersion_intervals(values, partial, anchor_fraction = .5)
  anchor_radius <- qnorm(.9875) * sqrt(1.4)
  expect_gte(bounded$lower, .2 - anchor_radius)
  expect_lte(bounded$upper, .2 + anchor_radius)
  expect_error(checked_dispersion_covariance(matrix(1, 2, 2)), "positive-definite")
  expect_error(arm_dispersion_intervals(values, partial, alpha = 0), "allocations")
})

cat("DISPERSION_INTERVAL_TESTS_PASSED\n")
