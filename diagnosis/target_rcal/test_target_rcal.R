#!/usr/bin/env Rscript
# Focused contracts for the diagnostic, separate from the frozen package suite.
.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source("diagnosis/target_rcal/target_rcal_helpers.R")
source("diagnosis/target_rcal/native_rcal_helpers.R")
suppressPackageStartupMessages(library(testthat))

test_that("AIPW reconstruction and oracle remainder are algebraically consistent", {
  eval_data <- list(n = 4L, A = c(0, 1, 0, 1), Y = c(0, 1, 1, 0))
  outcome <- c(.2, .4, .6, .8)
  probability <- c(.3, .4, .5, .6)
  treated <- target_prediction_fit(eval_data, 1L, outcome, probability)
  flipped <- eval_data
  flipped$A <- 1 - eval_data$A
  control <- target_prediction_fit(flipped, 0L, outcome, probability)
  expect_equal(treated$estimate, control$estimate, tolerance = 1e-14)
  expect_equal(treated$varphi_ot, control$varphi_ot, tolerance = 1e-14)
  expect_equal(mean(treated$varphi_ot), 0, tolerance = 1e-14)
  expect_equal(treated$variance, treated$V_ot / eval_data$n)
  expect_error(target_prediction_fit(eval_data, c(0L, 1L), outcome, probability),
               "invalid arm")
  expect_error(target_prediction_fit(eval_data, 1L, outcome, rep(0, 4)),
               "invalid arm")
})

for (rcal_backend in c("reference", "native")) {
  test_that(paste("RCAL", rcal_backend, "preserves arm labels and evaluation exclusion"), {
    set.seed(901)
    x <- matrix(rnorm(720), 240, 3)
    a <- rbinom(240, 1, plogis(.1 + .2 * x[, 1L]))
    y <- rbinom(240, 1, plogis(.2 + .4 * x[, 2L] + .3 * a))
    site <- function(rows) list(n = length(rows), Z_site = x[rows, , drop = FALSE],
                                W_outcome = x[rows, , drop = FALSE],
                                A = a[rows], Y = y[rows])
    train <- site(1:200)
    evaluation <- site(201:240)
    before <- .Random.seed
    first <- fit_target_rcal_predictions(train, evaluation, 1L, 501L, nlambda = 8L,
                                         rcal_backend = rcal_backend)
    expect_identical(.Random.seed, before)
    changed <- evaluation
    changed$Y <- 1 - changed$Y
    changed$A <- 1 - changed$A
    second <- fit_target_rcal_predictions(train, changed, 1L, 501L, nlambda = 8L,
                                          rcal_backend = rcal_backend)
    expect_identical(first, second)
    flipped_train <- train
    flipped_train$A <- 1 - train$A
    third <- fit_target_rcal_predictions(flipped_train, changed, 0L, 501L, nlambda = 8L,
                                         rcal_backend = rcal_backend)
    expect_identical(first$arm_probability, third$arm_probability)
    expect_identical(first$outcome, third$outcome)
    expect_true(all(is.finite(first$outcome)))
    expect_lt(abs(first$diagnostics$calibration_intercept), 1e-3)
    expect_identical(train, site(1:200))
  })
}

test_that("the legacy RCAL switch is demonstrably ignored in cross-fit mode", {
  set.seed(102)
  x <- matrix(rnorm(480), 160, 3)
  site <- list(n = 160L, W_outcome = x, Z_site = x,
               A = rep(0:1, 80), Y = rbinom(160, 1, .5))
  ordinary <- RoCE:::with_seed(302, RoCE:::fit_site_aipw(
    site, use_rcal = FALSE, use_crossfit = TRUE, n_folds = 4L
  ))
  requested_rcal <- RoCE:::with_seed(302, RoCE:::fit_site_aipw(
    site, use_rcal = TRUE, use_crossfit = TRUE, n_folds = 4L
  ))
  expect_identical(ordinary, requested_rcal)
})
