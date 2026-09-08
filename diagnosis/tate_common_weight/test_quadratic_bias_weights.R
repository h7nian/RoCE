library(testthat)
source(if (file.exists("quadratic_bias_weights.R")) "quadratic_bias_weights.R" else
         "diagnosis/tate_common_weight/quadratic_bias_weights.R")

weight_fixture <- function(n = 1000, discrepancy = .1) {
  # n_target * Var = 0.5 - 0.5*eta + 0.5*eta^2.
  list(V_ot = .5, V_t = .25, V_s = .25, C_ot = .25,
       C_cross = matrix(0, 1, 1), n_t = n, n_s = n,
       avg_target_est = .2, avg_source_est = .2+discrepancy)
}

test_that("one-source solution matches the exact penalized quadratic", {
  for (delta in c(0, -.1, .1, .5)) {
    fit <- .quadratic_bias_weights(weight_fixture(discrepancy = delta))
    expect_equal(fit$weights, .25/(.5+1000^.75*delta^2), tolerance = 1e-12)
    expect_lte(fit$relative_normal_equation_error, 1e-12)
  }
})

test_that("full outcome rescaling preserves learned weights", {
  m <- weight_fixture(); reference <- .quadratic_bias_weights(m)$weights
  for (factor in c(.01, 100)) {
    scaled <- m
    for (field in c("V_ot", "V_t", "V_s", "C_ot", "C_cross")) scaled[[field]] <- m[[field]]*factor^2
    for (field in c("avg_target_est", "avg_source_est")) scaled[[field]] <- m[[field]]*factor
    expect_equal(.quadratic_bias_weights(scaled)$weights, reference, tolerance = 1e-12)
  }
})

test_that("the deterministic rate examples satisfy the two necessary limits", {
  n <- c(1e4, 1e8, 1e12, 1e16)
  null_error <- vapply(n, function(size)
    abs(.quadratic_bias_weights(weight_fixture(size, 1/sqrt(size)))$weights-.5), numeric(1))
  biased_root_n <- vapply(n, function(size)
    sqrt(size)*.3*.quadratic_bias_weights(weight_fixture(size, .3))$weights, numeric(1))
  expect_true(all(diff(null_error) < 0))
  expect_true(all(diff(biased_root_n) < 0))
  expect_lt(tail(null_error, 1), .001)
  expect_lt(tail(biased_root_n, 1), .001)
})

test_that("source permutation and covariance coupling are retained", {
  m <- list(V_ot = .5, V_t = c(.3, .4), V_s = c(.2, .3), C_ot = c(.2, .1),
    C_cross = matrix(c(0, .1, .1, 0), 2), n_t = 800, n_s = c(700, 900),
    avg_target_est = .2, avg_source_est = c(.21, .35))
  original <- .quadratic_bias_weights(m)
  reordered <- m
  for (field in c("V_t", "V_s", "C_ot", "n_s", "avg_source_est")) reordered[[field]] <- m[[field]][2:1]
  reordered$C_cross <- m$C_cross[2:1, 2:1]
  expect_equal(.quadratic_bias_weights(reordered)$weights, rev(original$weights), tolerance = 1e-12)
  H <- original$variance_quadratic+diag(original$penalty)
  expect_equal(original$weights, drop(solve(H, -original$linear)), tolerance = 1e-12)
})

test_that("unsupported rate and covariance inputs fail rather than get ridged", {
  expect_error(.quadratic_bias_weights(weight_fixture(), .5), "strictly")
  expect_error(.quadratic_bias_weights(weight_fixture(), 1), "strictly")
  m <- weight_fixture(); m$C_ot <- 2
  expect_error(.quadratic_bias_weights(m), "curvature")
  m <- weight_fixture(); m$n_s <- 2.5
  expect_error(.quadratic_bias_weights(m), "sample size")
})
