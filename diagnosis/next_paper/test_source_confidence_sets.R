#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
source(arguments[[1L]])
library(testthat)

test_that("the sweep preserves disjoint intervals and closed endpoint votes", {
  expect_equal(unname(quorum_intervals(c(0, 1, 4), c(2, 3, 5), 2)), matrix(c(1, 2), 1))
  expect_equal(unname(quorum_intervals(c(0, 1), c(1, 2), 2)), matrix(c(1, 1), 1))
  expect_equal(unname(quorum_intervals(c(0, 0, 3, 3), c(1, 1, 4, 4), 2)),
               matrix(c(0, 3, 1, 4), 2))
  expect_equal(nrow(quorum_intervals(c(0, 3), c(1, 4), 2)), 0L)
  expect_error(quorum_intervals(c(2, 0), c(1, 4), 1), "ordered")
})

test_that("the endpoint sweep agrees with direct votes across randomized intervals", {
  set.seed(941)
  for (iteration in 1:30) {
    lower <- sample(-5:5, 12, replace = TRUE)
    upper <- lower + sample(0:5, 12, replace = TRUE)
    required <- sample(1:12, 1)
    intervals <- quorum_intervals(lower, upper, required)
    endpoints <- sort(unique(c(lower, upper)))
    points <- sort(unique(c(endpoints, (head(endpoints, -1) + tail(endpoints, -1)) / 2)))
    direct <- vapply(points, function(x) sum(lower <= x & x <= upper) >= required, logical(1))
    swept <- vapply(points, function(x) any(intervals[, 1] <= x & x <= intervals[, 2]), logical(1))
    expect_identical(swept, direct)
  }
})

test_that("the marginal bound has the stated probability budget", {
  calibration <- make_source_calibration(64, 48, 33)
  expect_equal(48 * calibration$marginal_tail / (48 - 33 + 1), calibration$source_alpha)
  expect_equal(calibration$source_alpha + calibration$anchor_alpha, .05)
  expect_error(make_source_calibration(4, 5, 3), "valid_minimum")
  expect_error(make_source_calibration(4, 3, 4), "votes_required")
  expect_error(make_source_calibration(4, 3, 2, source_correlation = .5), "not used")
})

test_that("Gaussian critical values match exact beta and dependence limits", {
  for (correlation in c(0, .5, 1)) {
    calibration <- make_source_calibration(8, 6, 5, method = "equicorrelated_gaussian",
                                          source_correlation = correlation)
    expect_equal(equicorrelated_vote_cdf(calibration$source_critical, 6, 5, correlation), .975,
                 tolerance = 1e-7)
    bound <- make_source_calibration(8, 6, 5)
    expect_lte(calibration$source_critical, bound$source_critical + 1e-8)
  }
  independent <- make_source_calibration(8, 6, 5, method = "equicorrelated_gaussian", source_correlation = 0)
  expect_equal(independent$source_critical, qnorm((1 + qbeta(.975, 5, 2)) / 2))
  single <- make_source_calibration(4, 1, 1, method = "equicorrelated_gaussian", source_correlation = .4)
  expect_equal(single$source_critical, qnorm(.9875))
})

test_that("inference is equivariant and retains its empty-set fallback", {
  calibration <- make_source_calibration(4, 3, 3)
  estimate <- c(-.2, .1, .4, 5)
  standard_error <- rep(.2, 4)
  fit <- source_search_interval(0, .5, estimate, standard_error, calibration)
  shifted <- source_search_interval(7, 1, 7 + 2 * estimate, 2 * standard_error, calibration)
  expect_equal(shifted$intervals, 7 + 2 * fit$intervals)
  permuted <- source_search_interval(0, .5, rev(estimate), rev(standard_error), calibration)
  expect_equal(permuted$intervals, fit$intervals)
  fallback <- source_search_interval(0, .1, rep(100, 4), rep(.1, 4), calibration)
  expect_true(fallback$used_anchor_fallback)
  expect_equal(c(fallback$lower, fallback$upper), c(-1, 1) * qnorm(.9875) * .1)
  expect_error(source_search_interval(0, .1, estimate, c(0, .1, .1, .1), calibration), "positive")
})

test_that("no valid-source assumption gives exactly the target-only interval", {
  calibration <- make_source_calibration(6, 0, 0)
  fit <- source_search_interval(.2, .03, calibration = calibration)
  expect_equal(c(fit$lower, fit$upper), .2 + c(-1, 1) * qnorm(.975) * .03)
  expect_identical(calibration$method, "target_only")
  expect_equal(fit$source_components, 0L)
  expect_error(make_source_calibration(6, 0, 1), "votes_required")
})

test_that("Poisson-binomial tails agree with exhaustive outcomes", {
  probabilities <- c(0, .2, .7, 1, .4)
  outcomes <- as.matrix(expand.grid(rep(list(0:1), length(probabilities))))
  probability <- apply(outcomes, 1L, function(row) prod(ifelse(row == 1, probabilities, 1 - probabilities)))
  for (threshold in 0:length(probabilities)) {
    expect_equal(poisson_binomial_tail(probabilities, threshold),
                 sum(probability[rowSums(outcomes) >= threshold]), tolerance = 1e-13)
  }
  expect_error(poisson_binomial_tail(c(.2, 1.1), 1), "probabilities")
})

test_that("shared-factor calibration has the required exact reductions", {
  private <- c(.1, .2, .5, 1, 2, 4)
  for (critical in c(1, 2, 3)) {
    expect_equal(shared_target_vote_cdf(critical, 4, 3, 0, private),
                 equicorrelated_vote_cdf(critical, 4, 3, 0), tolerance = 1e-12)
    expect_equal(shared_target_vote_cdf(critical, 4, 3, .4, rep(.6, 6)),
                 equicorrelated_vote_cdf(critical, 4, 3, .4), tolerance = 1e-12)
  }
  calibration <- make_source_calibration(6, 4, 3, method = "shared_target_gaussian",
                                         shared_variance = 0, source_private_variances = private)
  reference <- make_source_calibration(6, 4, 3, method = "equicorrelated_gaussian", source_correlation = 0)
  expect_equal(calibration$source_critical, reference$source_critical, tolerance = 2e-7)
})

test_that("the conditional adversary lower-bounds every fixed valid subset", {
  private <- c(.2, .5, 1.2, 2)
  lower_bound <- shared_target_vote_cdf(2, 3, 2, .6, private)
  subsets <- combn(seq_along(private), 3, simplify = FALSE)
  exact <- vapply(subsets, function(index) shared_target_vote_cdf(2, 3, 2, .6, private[index]), numeric(1))
  expect_lte(lower_bound, min(exact) + 1e-8)
  expect_gt(lower_bound, 0)
  expect_lt(lower_bound, 1)
})

test_that("factor calibration is scale invariant and enforces its SE specification", {
  private <- c(.2, .5, 1.2, 2)
  calibration <- make_source_calibration(4, 3, 2, method = "shared_target_gaussian",
                                         shared_variance = .6, source_private_variances = private)
  scaled <- make_source_calibration(4, 3, 2, method = "shared_target_gaussian",
                                    shared_variance = 2.4, source_private_variances = 4 * private)
  expect_equal(calibration$source_critical, scaled$source_critical, tolerance = 2e-7)
  expect_lte(calibration$source_critical, make_source_calibration(4, 3, 2)$source_critical)
  expect_equal(shared_target_vote_cdf(calibration$factor_critical, 3, 2, .6, private), .975,
               tolerance = 2e-6)
  expect_error(source_search_interval(0, 1, rep(0, 4), rep(1, 4), calibration), "standard errors")
  fitted <- source_search_interval(0, 1, c(-.3, .2, .4, 1), sqrt(.6 + private), calibration)
  expect_true(fitted$lower < fitted$upper)
})

cat("SOURCE_CONFIDENCE_SET_TESTS_PASSED\n")
