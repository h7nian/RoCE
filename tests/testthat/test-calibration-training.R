test_that("compact calibration preserves fold-specific initial predictors", {
  designs <- list(matrix(1:6, 3, 2), matrix(7:14, 4, 2))
  coefficients <- list(c(.2, .3, -.1), c(-.5, -.2, .6))
  block <- .calibration_plugin(designs, coefficients, "block")
  compact <- .calibration_plugin(designs, coefficients, "compact")
  expected <- c(drop(cbind(1, designs[[1]]) %*% coefficients[[1]]),
                 drop(cbind(1, designs[[2]]) %*% coefficients[[2]]))
  expect_equal(drop(cbind(1, block$design) %*% block$coefficients), expected, tolerance = 1e-14)
  expect_equal(drop(cbind(1, compact$design) %*% compact$coefficients), expected, tolerance = 1e-14)
  expect_lt(length(compact$design), length(block$design))
})

test_that("fold-summed fitting agrees with the pooled weighted least-squares objective", {
  set.seed(49)
  x <- matrix(rnorm(720), 240, 3)
  a <- rep(0:1, 120)
  y <- .5 + x[, 1] - .2 * x[, 2] + rep(c(-.5, .5), each = 120) + rnorm(240, sd = .2)
  indices <- list(k2_2 = 1:120, k2_3 = 121:240)
  blocks <- lapply(indices, function(rows) list(
    W_outcome = x[rows, , drop = FALSE], Z_site = x[rows, , drop = FALSE],
    A = a[rows], Y = y[rows], n = length(rows), original_idx = rows
  ))
  initial_outcome <- list(k2_2 = c(.1, .1, .2, 0), k2_3 = c(-.1, .2, .1, 0))
  initial_weight <- list(k2_2 = c(-.7, .1, -.1, 0), k2_3 = c(-.5, -.1, .2, 0))
  expected_tilt <- c(-.6, .1, -.2, .05)
  design <- cbind(1, x)
  linear_moment <- colMeans(design * (a * exp(-drop(design %*% expected_tilt))))
  fits <- lapply(c("block", "compact"), function(layout) {
    .fit_fold_summed_calibration(
      blocks, initial_outcome, initial_weight, linear_moment,
      site = "s1", training_folds = 2:3, A_val = 1L,
      family_int = 0L, link_int = 0L, M_tau = 5,
      nlambda = 5L, max_iter = 1000L, lambda_rule = "min", layout = layout,
      weight_lambda = 0, outcome_lambda = 0, tol = 1e-10
    )
  })
  initial_predictor <- unlist(Map(function(rows, coefficients) {
    drop(design[rows, , drop = FALSE] %*% coefficients)
  }, indices, initial_weight), use.names = FALSE)
  weights <- a * exp(-initial_predictor)
  expected_outcome <- drop(solve(crossprod(design, design * weights), crossprod(design, y * weights)))
  for (fit in fits) {
    expect_equal(as.numeric(fit$weight), expected_tilt, tolerance = 1e-7)
    expect_equal(as.numeric(fit$outcome), expected_outcome, tolerance = 1e-7)
  }
  expect_equal(as.numeric(fits[[1]]$weight), as.numeric(fits[[2]]$weight), tolerance = 1e-9)
  expect_equal(as.numeric(fits[[1]]$outcome), as.numeric(fits[[2]]$outcome), tolerance = 1e-9)
})

test_that("binary calibration layouts preserve the weighted derivative moments and CV choice", {
  set.seed(53)
  x <- matrix(rnorm(720), 240, 3)
  design <- cbind(1, x)
  a <- rep(0:1, 120)
  y <- rbinom(240, 1, plogis(.2 + .3 * x[, 1]))
  indices <- list(k2_2 = 1:120, k2_3 = 121:240)
  blocks <- lapply(indices, function(rows) list(
    W_outcome = x[rows, , drop = FALSE], Z_site = x[rows, , drop = FALSE],
    A = a[rows], Y = y[rows], n = length(rows)
  ))
  initial_outcome <- list(k2_2 = c(.3, .2, .1, 0), k2_3 = c(-.3, .4, -.1, 0))
  initial_weight <- list(k2_2 = c(-.7, .1, -.1, 0), k2_3 = c(-.5, -.1, .2, 0))
  initial_mean <- unlist(Map(function(rows, coefficients) {
    plogis(drop(design[rows, , drop = FALSE] %*% coefficients))
  }, indices, initial_outcome), use.names = FALSE)
  moment <- colMeans(design * (a * exp(-drop(design %*% c(-.6, .1, -.2, .05))) *
                                initial_mean * (1 - initial_mean)))
  fits <- lapply(c("block", "compact"), function(layout) {
    .fit_fold_summed_calibration(
      blocks, initial_outcome, initial_weight, moment,
      site = "s1", training_folds = 2:3, A_val = 1L,
      family_int = 1L, link_int = 1L, M_tau = 5,
      nlambda = 8L, max_iter = 1000L, lambda_rule = "min", layout = layout,
      tol = 1e-10
    )
  })
  for (model in c("weight", "outcome")) {
    expect_equal(attr(fits[[1]][[model]], "lambda_used"),
                 attr(fits[[2]][[model]], "lambda_used"), tolerance = 1e-12)
    expect_equal(as.numeric(fits[[1]][[model]]), as.numeric(fits[[2]][[model]]), tolerance = 1e-8)
  }
})
