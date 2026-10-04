#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
source(arguments[[1L]])
library(testthat)

target <- cbind(target_anchor = c(1, 0, 2, -1, 3), s1 = c(.1, .3, .5, .8, .6),
                s2 = c(-.1, .3, .8, .7, .9))
residuals <- list(s1 = c(-1, .5, 0, .2), s2 = c(.2, -.4, .3, .7, -.1, .9))
indices <- list(list(target = c(1L, 3L, 5L), source = list(c(1L, 4L), c(1L, 2L))),
                list(target = c(2L, 4L), source = list(c(2L, 3L), 3:6)))
folds <- lapply(indices, function(index) {
  target_values <- target[index$target, , drop = FALSE]
  source_values <- Map(function(values, positions) values[positions], residuals, index$source)
  list(target_idx = index$target, source_idx = index$source,
    fold_target_estimate = mean(target_values[, 1L]),
    varphi_ot = target_values[, 1L] - mean(target_values[, 1L]),
    mu_pred_ts = colMeans(target_values[, -1L, drop = FALSE]),
    zeta_components = lapply(2:3, function(column) target_values[, column] - mean(target_values[, column])),
    delta_ts = vapply(source_values, mean, numeric(1L)),
    V_s = vapply(source_values, function(values) mean((values - mean(values))^2), numeric(1L)),
    xi_components = lapply(source_values, function(values) values - mean(values)))
})
data_split <- list(t = list(n = 5L), s1 = list(n = 4L), s2 = list(n = 6L))
fitted <- list(source_estimates = colMeans(target[, -1L]) + vapply(residuals, mean, numeric(1L)),
  target_only = list(estimate = mean(target[, 1L]), variance = sum((target[, 1L] - mean(target[, 1L]))^2) / 25),
  n_folds = 2L, intermediates = list(fold_info = folds))

test_that("unequal folds reconstruct each site's mean and covariance exactly", {
  summary <- extract_candidate_score_summary(fitted, data_split, keep_scores = TRUE)
  expected_covariance <- stats::cov(target) * 4 / 25
  diag(expected_covariance)[-1L] <- diag(expected_covariance)[-1L] +
    vapply(residuals, function(values) stats::var(values) * (length(values) - 1) / length(values)^2, numeric(1))
  expect_equal(summary$covariance, expected_covariance, tolerance = 1e-14)
  expect_equal(summary$score_records$target, target, tolerance = 1e-14)
  expect_equal(summary$score_records$source, residuals, tolerance = 1e-14)
  expect_equal(summary$estimates[-1L], fitted$source_estimates)
  expect_equal(summary$anchor_variance_identity_error, 0, tolerance = 1e-14)
  expect_equal(summary$source_message_identity_error, 0, tolerance = 1e-14)
})

test_that("summary pooling retains between-fold mean variation", {
  pooled <- pool_source_fold_moments(matrix(c(0, 2), 2, 1), matrix(0, 2, 1), matrix(c(2, 3), 2, 1))
  values <- c(0, 0, 2, 2, 2)
  expect_equal(unname(pooled$means), mean(values))
  expect_equal(unname(pooled$variances), sum((values - mean(values))^2) / length(values)^2)
  expect_gt(pooled$variances[1], 0)
})

test_that("rank one alone is not mistaken for equal target components", {
  base <- c(-2, -1, 0, 1, 2)
  equal <- common_target_diagnostics(cbind(s1 = base, s2 = base), c("s1", "s2"))
  unequal <- common_target_diagnostics(cbind(s1 = base, s2 = 2 * base), c("s1", "s2"))
  expect_equal(equal$relative_covariance_error, 0)
  expect_equal(unequal$leading_eigenvalue_fraction, 1, tolerance = 1e-14)
  expect_gt(unequal$relative_covariance_error, .1)
  single <- common_target_diagnostics(cbind(s1 = base, s2 = 2 * base), "s1")
  expect_true(is.na(single$relative_covariance_error))
})

test_that("bad fold provenance and inconsistent saved means are rejected", {
  duplicate <- fitted
  duplicate$intermediates$fold_info[[2L]]$target_idx <- c(1L, 4L)
  expect_error(extract_candidate_score_summary(duplicate, data_split), "exactly once")
  shifted <- fitted
  shifted$source_estimates[1] <- shifted$source_estimates[1] + .01
  expect_error(extract_candidate_score_summary(shifted, data_split), "disagree")
  expect_error(extract_candidate_score_summary(fitted, data_split, comparison_sources = "s3"), "Comparison")
})

cat("CANDIDATE_SCORE_SUMMARY_TESTS_PASSED\n")
