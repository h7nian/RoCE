#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Usage: research_source_directory")
for (name in c("source_confidence_sets.R", "candidate_score_summary.R",
               "arm_candidate_score_summary.R", "arm_pair_confidence_sets.R")) {
  source(file.path(arguments[[1L]], name))
}
library(testthat)

make_arm_fixture <- function(target, residuals, indices) {
  folds <- lapply(indices, function(index) {
    target_values <- target[index$target, , drop = FALSE]
    source_values <- Map(function(values, rows) values[rows], residuals, index$source)
    list(target_idx = index$target, source_idx = index$source,
      fold_target_estimate = mean(target_values[, 1L]),
      varphi_ot = target_values[, 1L] - mean(target_values[, 1L]),
      mu_pred_ts = colMeans(target_values[, -1L, drop = FALSE]),
      zeta_components = lapply(2:ncol(target_values), function(column)
        target_values[, column] - mean(target_values[, column])),
      delta_ts = vapply(source_values, mean, numeric(1L)),
      V_s = vapply(source_values, function(values) mean((values - mean(values))^2), numeric(1L)),
      xi_components = lapply(source_values, function(values) values - mean(values)))
  })
  list(source_estimates = colMeans(target[, -1L, drop = FALSE]) + vapply(residuals, mean, numeric(1L)),
    target_only = list(estimate = mean(target[, 1L]), variance = as.numeric(score_mean_covariance(target[, 1L]))),
    n_folds = length(folds), intermediates = list(fold_info = folds))
}

test_that("joint arm extraction retains cross-arm and between-fold covariance", {
  target1 <- cbind(target_anchor = c(1, 0, 2, -1, 3), s1 = c(.1, .3, .5, .8, .6),
                   s2 = c(-.1, .3, .8, .7, .9))
  target0 <- target1 / 2 + c(.1, -.3, .2, -.4, .8)
  residual1 <- list(s1 = c(-1, .5, 0, .2), s2 = c(.2, -.4, .3, .7, -.1, .9))
  residual0 <- list(s1 = c(.4, -.7, .3, 0), s2 = c(.4, .8, -.5, .2, .4, -.2))
  indices <- list(list(target = c(1L, 3L, 5L), source = list(c(1L, 4L), c(1L, 2L))),
                  list(target = c(2L, 4L), source = list(c(2L, 3L), 3:6)))
  mu1 <- make_arm_fixture(target1, residual1, indices)
  mu0 <- make_arm_fixture(target0, residual0, indices)
  fitted <- make_arm_fixture(target1 - target0, Map(`-`, residual1, residual0), indices)
  fitted$arm_results <- list(mu1 = mu1, mu0 = mu0)
  data <- list(t = list(n = 5L), s1 = list(n = 4L), s2 = list(n = 6L))
  observed <- extract_joint_arm_score_summary(fitted, data, keep_scores = TRUE)
  expected <- stats::cov(cbind(target1, target0)) * 4 / 25
  for (index in 1:2) {
    values <- cbind(residual1[[index]], residual0[[index]])
    positions <- c(index + 1L, index + 4L)
    expected[positions, positions] <- expected[positions, positions] +
      stats::cov(values) * (nrow(values) - 1) / nrow(values)^2
  }
  expect_equal(unname(observed$covariance), unname(expected), tolerance = 1e-14)
  expect_equal(observed$covariance_identity_error, 0, tolerance = 1e-14)
  expect_equal(observed$source_message_identity_error, 0, tolerance = 1e-14)
  pairs <- arm_pair_moments(observed$means, observed$covariance)
  for (index in seq_len(nrow(pairs))) {
    contrast <- numeric(6L)
    contrast[match(pairs$treated[index], rownames(observed$means))] <- 1
    contrast[3 + match(pairs$control[index], rownames(observed$means))] <- -1
    expect_equal(pairs$variance[index], drop(t(contrast) %*% expected %*% contrast), tolerance = 1e-14)
  }
  bad <- fitted
  bad$arm_results$mu0$intermediates$fold_info[[1L]]$source_idx[[1L]] <- c(4L, 1L)
  expect_error(extract_joint_arm_score_summary(bad, data), "observation order")
  expect_null(extract_joint_arm_score_summary(fitted, data)$score_records)
})

test_that("cross-arm pooling includes different fold means", {
  mu1 <- c(0, 0, 2, 2, 2)
  mu0 <- c(3, 3, -1, -1, -1)
  observed <- pool_source_fold_cross_covariance(matrix(c(0, 2), 2), matrix(c(3, -1), 2),
    matrix(0, 2), matrix(c(2, 3), 2))
  expect_equal(as.numeric(observed), stats::cov(mu1, mu0) * 4 / 25, tolerance = 1e-14)
  expect_lt(observed[1L], 0)
})

test_that("pair confidence sets handle degeneracy and explicit validity counts", {
  means <- cbind(mu1 = c(.6, .6, .6), mu0 = c(.4, .4, .4))
  rownames(means) <- c("target_anchor", "s1", "s2")
  covariance <- matrix(0, 6, 6)
  labels <- c(paste0("mu1:", rownames(means)), paste0("mu0:", rownames(means)))
  dimnames(covariance) <- list(labels, labels)
  result <- arm_pair_search_interval(means, covariance, 1, 2)
  expect_equal(c(result$lower, result$upper), c(.2, .2), tolerance = 1e-14)
  expect_equal(result$valid_pairs, 6)
  expect_equal(result$votes_required, 3)
  expect_equal(result$pair_tail, .025 * 4 / 6)
  diag(covariance) <- .01
  anchor <- arm_pair_search_interval(means, covariance, 0, 0)
  expect_true(anchor$anchor_only)
  expect_equal(c(anchor$lower, anchor$upper), .2 + c(-1, 1) * qnorm(.975) * sqrt(.02))
  expect_error(arm_pair_search_interval(means, covariance, 3, 1), "valid_mu1")
  expect_error(arm_pair_search_interval(means, covariance, 1, 1, votes_required = 5), "votes_required")
  bad <- covariance
  bad[1, 2] <- bad[2, 1] <- 1
  expect_error(arm_pair_moments(means, bad), "positive semidefinite")
  expect_error(arm_pair_moments(means, covariance[6:1, 6:1]), "labels")
})

cat("ARM_PAIR_CONFIDENCE_TESTS_PASSED\n")
