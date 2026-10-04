test_that("both nuisance solvers reject an unbounded tilt stopped at the coefficient limit", {
  # Direction d=(-1,1) has Xd > 0 and target_moment'd + lambda*|d1| < 0.
  # Thus the objective has no finite minimizer. A bounded coefficient returned
  # by a stopped optimizer is not a valid stationary solution.
  x <- matrix(seq(1.1, 2, length.out = 40), ncol = 1L)
  previous <- set_nuisance_solver_cpp("coordinate_descent")
  on.exit(set_nuisance_solver_cpp(previous), add = TRUE)
  for (solver in c("coordinate_descent", "proximal_newton")) {
    set_nuisance_solver_cpp(solver)
    fit <- fit_initial_density_ratio_cpp(x, rep(1, 40), c(1, -1),
      lambda = .1, max_iter = 1000L, tol = 1e-7, A_val = 1L, M_tau = 5,
      warm_start = numeric())
    expect_false(fit$converged)
    expect_gt(fit$kkt_residual, fit$kkt_threshold)
  }
})

test_that("reported nuisance KKT residuals match the actual calibrated objective", {
  set.seed(1070)
  x <- matrix(rnorm(900), 300, 3)
  a <- rep(0:1, 150)
  moment <- c(1, .05, -.03, .01)
  for (solver in c("coordinate_descent", "proximal_newton")) {
    fit <- fit_initial_density_ratio(x, a, moment, lambda = .02,
      M_tau = 5, tol = 1e-8, nuisance_solver = solver)
    design <- cbind(1, x)
    eta <- drop(design %*% fit)
    gradient <- moment - colMeans(design * a * exp(-pmax(-5, pmin(5, eta))))
    residual <- max(abs(c(gradient[1L], ifelse(fit[-1L] != 0,
      gradient[-1L] + .02 * sign(fit[-1L]), pmax(abs(gradient[-1L]) - .02, 0)))))
    expect_lt(abs(attr(fit, "kkt_residual") - residual), 1e-12)
    expect_true(attr(fit, "converged"))
    expect_lte(residual, 1e-8)
  }
})

test_that("Newton retains a resolvable model decrease at strict nuisance tolerance", {
  set.seed(1081)
  data <- split_data_by_site(generate_simulation_data(
    n_total = 3000L, K = 2L, p = 4L, config = "C1", dgp_type = "face",
    outcome_type = "binary", n_target = 1000L, n_source_sizes = c(1000L, 1000L),
    warn_ignored = FALSE))
  folds <- build_crossfit_folds(data, 4L)
  source_train <- combine_folds(folds$source_folds$s1, 2:3)
  target_train <- combine_folds(folds$target_folds, 2:3)
  fits <- lapply(c("coordinate_descent", "proximal_newton"), function(solver) {
    .fit_nuisance_training_subset("fit_initial_density_ratio", list(
      Z_site = source_train$Z_site, A = source_train$A,
      mean_phi = c(1, colMeans(target_train$Z_site)), A_val = 1L,
      nlambda = 8L, M_tau = 5, tol = 1e-10, nuisance_solver = solver),
      "s1", 2:3)
  })
  for (fit in fits) expect_lte(attr(fit, "kkt_residual"), 1e-10)
  expect_equal(as.numeric(fits[[1L]]), as.numeric(fits[[2L]]), tolerance = 1e-8)
  expect_equal(attr(fits[[1L]], "lambda_used"), attr(fits[[2L]], "lambda_used"), tolerance = 1e-12)
})
