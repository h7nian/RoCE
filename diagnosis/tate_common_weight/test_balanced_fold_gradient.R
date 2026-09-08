library(testthat)
source(if (file.exists("balanced_fold_gradient.R")) "balanced_fold_gradient.R" else
         "diagnosis/tate_common_weight/balanced_fold_gradient.R")

test_that("shared treatment-total variance replaces independent fold composition", {
  fold <- rep(1:3, each = 10)
  A <- rep(c(rep(1, 3), rep(0, 7)), 3)
  coefficients <- c(1, 2, 4)
  score <- coefficients[fold]*A
  gradient <- (score-mean(score))/length(score)
  result <- .balanced_arm_fold_gradient(gradient, A, fold)
  expected <- mean(coefficients)^2*.3*.7/30
  expect_equal(result$variance, expected, tolerance = 1e-12)
  expect_equal(result$treatment_slope, mean(coefficients), tolerance = 1e-12)
  expect_equal(result$within_cell, rep(0, 30), tolerance = 1e-12)
  expect_gt(sum(gradient^2), expected)
  expect_equal(.balanced_arm_fold_gradient((A-.3)/30, A, fold)$gradient,
               (A-.3)/30, tolerance = 1e-12)
})

test_that("the map preserves contrasts and joint covariance", {
  fold <- rep(1:3, each = 10); A <- rep(c(rep(1, 3), rep(0, 7)), 3)
  residual <- sin(seq_along(A))
  cell <- interaction(fold, A)
  residual <- residual-ave(residual, cell)
  u <- (c(1, 2, 4)[fold]*A+residual)/30
  v <- (c(2, -.5, 1)[fold]*A-.3*residual)/30
  first <- .balanced_arm_fold_gradient(u, A, fold)
  second <- .balanced_arm_fold_gradient(v, A, fold)
  contrast <- .balanced_arm_fold_gradient(u-v, A, fold)
  expect_equal(contrast$gradient, first$gradient-second$gradient, tolerance = 1e-12)
  expected_covariance <- sum(residual*(-.3*residual))/30^2+
    mean(c(1, 2, 4))*mean(c(2, -.5, 1))*.3*.7/30
  expect_equal(sum(first$gradient*second$gradient), expected_covariance, tolerance = 1e-12)
  expect_equal(contrast$variance, first$variance+second$variance-
                 2*sum(first$gradient*second$gradient), tolerance = 1e-12)
  order <- rev(seq_along(A))
  expect_equal(.balanced_arm_fold_gradient(u[order], A[order], fold[order])$gradient,
               first$gradient[order], tolerance = 1e-12)
})

test_that("the specific fold design is enforced", {
  expect_error(.balanced_arm_fold_gradient(rep(0, 6), c(1, 1, 0, 0, 0, 0),
                                           c(1, 1, 1, 2, 2, 2)), "nonempty")
  expect_error(.balanced_arm_fold_gradient(c(0, NA), 0:1, c(1, 1)), "invalid")
  expect_error(.balanced_arm_fold_gradient(c(0, 0), 0:1, c(2, 2)), "consecutive")
})
