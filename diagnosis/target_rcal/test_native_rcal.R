.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
suppressPackageStartupMessages(library(testthat))
source("diagnosis/target_rcal/native_rcal_helpers.R")

test_that("the compiled calibrated propensity objective matches RCAL", {
  set.seed(196)
  x <- scale(matrix(rnorm(1000), 200, 5))
  indicator <- rbinom(200, 1, plogis(.2 + .5 * x[, 1L] - .3 * x[, 2L]))
  for (lambda in c(.02, .08, .2)) {
    reference <- RCAL::glm.regu(y = indicator, x = x, loss = "cal",
                                rhos = rep(lambda, ncol(x)),
                                n.iter = 1000L, eps = 1e-10)
    candidate <- fit_rcal_calibration(x, indicator, lambda)
    expect_gt(reference$conv, 0)
    expect_true(candidate$converged)
    expect_lt(candidate$kkt_error, 1e-5)
    # Absolute tolerance: all.equal's relative scaling exaggerates differences
    # near zero. RCAL and the compiled solver use different stopping criteria.
    expect_lt(max(abs(candidate$coefficients - c(reference$inter, reference$bet))), 2e-5)
  }
})

test_that("weighted glmnet and RCAL minimize the same fixed-penalty outcome loss", {
  set.seed(196)
  x <- scale(matrix(rnorm(1000), 200, 5))
  y <- rbinom(200, 1, plogis(.2 + .5 * x[, 1L]))
  weights <- exp(.3 * x[, 2L])
  lambda <- .04
  reference <- RCAL::glm.regu(y = y, x = x, iw = weights, loss = "ml",
                              rhos = rep(lambda, ncol(x)), n.iter = 1000L,
                              eps = 1e-10)
  candidate <- glmnet::glmnet(x = x, y = y, weights = weights, family = "binomial",
                              lambda = lambda, standardize = FALSE, thresh = 1e-12)
  expect_gt(reference$conv, 0)
  expect_lt(max(abs(as.numeric(coef(candidate)) - c(reference$inter, reference$bet))), 2e-5)
})
