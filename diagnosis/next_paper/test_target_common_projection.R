#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
source(file.path(arguments[[1L]], "candidate_score_summary.R"))
source(file.path(arguments[[1L]], "target_common_projection.R"))
library(testthat)
set.seed(94103)
design <- matrix(rnorm(1200), 300, 4)
target_data <- list(n = 300L, W_outcome = design, Z_site = design,
  A = rep(0:1, 150), Y = rbinom(300, 1, plogis(.3 + design[, 1] - .4 * design[, 2])))
labels <- rep(1:5, length.out = 300)
coefficients <- c(0, .2, -.1, 0, 0)

test_that("held-out outcomes do not affect their projection model", {
  original <- fit_target_projection_fold(target_data, labels, 1, coefficients, 1, nlambda = 12)
  changed <- target_data
  changed$Y[labels == 1] <- 1 - changed$Y[labels == 1]
  refitted <- fit_target_projection_fold(changed, labels, 1, coefficients, 1, nlambda = 12)
  expect_equal(original$prediction, refitted$prediction, tolerance = 1e-14)
  expect_equal(original$coefficients, refitted$coefficients, tolerance = 1e-14)
  expect_false(any(original$training_indices %in% original$evaluation_indices))
})

test_that("constant propensity weights reproduce ordinary weighted-likelihood fitting", {
  ipw <- fit_target_projection_fold(target_data, labels, 1, rep(0, 5), 0, nlambda = 12)
  unweighted <- fit_target_projection_fold(target_data, labels, 1, rep(0, 5), 0,
                                           weighting = "unweighted", nlambda = 12)
  expect_equal(ipw$prediction, unweighted$prediction, tolerance = 1e-12)
  expect_equal(ipw$weight_range, c(2, 2))
  expect_equal(ipw$propensity_clipped_fraction, 0)
})

cat("TARGET_COMMON_PROJECTION_TESTS_PASSED\n")
