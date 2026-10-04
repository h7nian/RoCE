test_that("initial weight CV extends a wholly failed grid without changing folds", {
  design <- matrix(c(1, rep(0, 39)), ncol = 1L)
  treatment <- rep(1L, 40)
  moment <- c(1, .02)
  grid <- build_lambda_grid(compute_lambda_max_initial_dr(design, treatment, moment),
    lambda_min_ratio = 1e-4, nlambda = 12L)
  previous <- .set_nuisance_solver("proximal_newton")
  on.exit(.restore_nuisance_solver(previous), add = TRUE)
  set.seed(919)
  expect_error(select_lambda_cv_initial_density_ratio_cpp(design, treatment, moment,
    grid, 5L, 10000L, 1e-10, 1L, 2, use_kkt_certificate = TRUE), "no lambda converged")
  original_rng <- .Random.seed
  set.seed(919)
  fitted <- fit_initial_density_ratio(design, treatment, moment, nlambda = 12L,
    M_tau = 2, tol = 1e-10, nuisance_solver = "proximal_newton", nuisance_cv_certificate = TRUE)
  expect_identical(.Random.seed, original_rng)
  retry <- attr(fitted, "cv_grid_retry")
  expect_type(retry, "list")
  expect_gt(retry$added_lambdas, 0L)
  expect_equal(retry$original_grid, as.numeric(grid), tolerance = 0)
  expect_true(all(as.numeric(grid) %in% retry$retried_grid))
  expect_equal(min(retry$retried_grid), min(grid), tolerance = 0)
  expect_true(isTRUE(attr(fitted, "converged")))
  set.seed(919)
  manual <- select_lambda_cv_initial_density_ratio_cpp(design, treatment, moment,
    retry$retried_grid, 5L, 10000L, 1e-10, 1L, 2, use_kkt_certificate = TRUE)
  expect_equal(attr(fitted, "lambda_used"), manual$lambda_min, tolerance = 0)
  for (solver in c("proximal_newton", "coordinate_descent")) {
    set.seed(919)
    plain <- fit_initial_density_ratio(design, treatment, moment, nlambda = 12L,
      M_tau = 2, tol = 1e-10, nuisance_solver = solver, nuisance_cv_certificate = FALSE)
    expect_equal(as.numeric(plain), as.numeric(fitted), tolerance = 1e-8)
    expect_equal(attr(plain, "lambda_used"), attr(fitted, "lambda_used"), tolerance = 0)
  }
})

test_that("successful and fixed-penalty initial fits do not extend their grid", {
  design <- matrix(seq(-1, 1, length.out = 60L), ncol = 1L)
  treatment <- rep(1L, 60L); moment <- c(1, .3)
  previous <- .set_nuisance_solver("proximal_newton")
  on.exit(.restore_nuisance_solver(previous), add = TRUE)
  grid <- build_lambda_grid(compute_lambda_max_initial_dr(design, treatment, moment),
    lambda_min_ratio = 1e-4, nlambda = 12L)
  set.seed(920)
  reference <- select_lambda_cv_initial_density_ratio_cpp(design, treatment, moment,
    grid, 5L, 10000L, 1e-10, 1L, 2, use_kkt_certificate = TRUE)
  original_rng <- .Random.seed
  set.seed(920)
  fitted <- fit_initial_density_ratio(design, treatment, moment, nlambda = 12L,
    M_tau = 2, tol = 1e-10, nuisance_solver = "proximal_newton", nuisance_cv_certificate = TRUE)
  expect_null(attr(fitted, "cv_grid_retry"))
  expect_equal(attr(fitted, "lambda_used"), reference$lambda_min, tolerance = 0)
  expect_identical(.Random.seed, original_rng)
  fixed <- fit_initial_density_ratio(design, treatment, moment, lambda = .1,
    M_tau = 2, tol = 1e-10, nuisance_solver = "proximal_newton")
  expect_null(attr(fixed, "cv_grid_retry"))
  expect_identical(.Random.seed, original_rng)
})

test_that("derivative-weighted density CV recovers the same rare-support failure", {
  design <- matrix(c(1, rep(0, 39)), ncol = 1L)
  treatment <- rep(1L, 40L); moment <- c(.25, .005)
  previous <- .set_nuisance_solver("proximal_newton")
  on.exit(.restore_nuisance_solver(previous), add = TRUE)
  grid <- build_lambda_grid(compute_lambda_max_refined_dr(design, treatment, moment,
    c(0, 0), W_outcome = design, calibrated = TRUE, M_tau = 2), 1e-4, 12L)
  set.seed(921)
  expect_error(select_lambda_cv_calibrated_density_ratio_cpp(design, treatment, moment,
    c(0, 0), grid, 5L, 10000L, 1e-10, 2, design, 1L, 1L, 1L,
    truncate_initial_outcome = FALSE, use_kkt_certificate = TRUE), "no lambda converged")
  original_rng <- .Random.seed
  for (calibrated in c(TRUE, FALSE)) for (solver in c("proximal_newton", "coordinate_descent")) {
    set.seed(921)
    fitted <- fit_unified_density_ratio(design, treatment, moment, c(0, 0),
      W_outcome = design, nlambda = 12L, M_tau = 2, calibrated = calibrated,
      truncate_initial_outcome = FALSE, tol = 1e-10,
      nuisance_solver = solver, nuisance_cv_certificate = TRUE)
    expect_identical(.Random.seed, original_rng)
    expect_type(attr(fitted, "cv_grid_retry"), "list")
    expect_true(isTRUE(attr(fitted, "converged")))
    expect_equal(as.numeric(fitted), c(0, 0), tolerance = 1e-8)
  }
})

test_that("density CV retry propagates other errors and cannot change a second failure", {
  x <- matrix(c(0, 1), ncol = 1L)
  calls <- 0L
  run_cv <- function(grid) { calls <<- calls + 1L; stop("unrelated failure") }
  expect_error(.density_ratio_cv_with_retry(run_cv, .01, x, c(1, 1), c(1, .5), 1),
    "unrelated failure")
  expect_identical(calls, 1L)
  calls <- 0L
  run_cv <- function(grid) { calls <<- calls + 1L
    stop("aggregate_cv_results: no lambda converged with a finite validation loss across every CV fold.") }
  expect_error(.density_ratio_cv_with_retry(run_cv, .01, x, c(1, 1), c(1, .5), 1), "no lambda converged")
  expect_identical(calls, 2L)
})
