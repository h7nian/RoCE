test_that("CV certificate controls are scoped, validated and restored on failure", {
  prior <- options(RoCE.nuisance_cv_certificate = NULL)
  on.exit(options(prior), add = TRUE)
  expect_false(.nuisance_cv_certificate_enabled())
  inherited <- .set_nuisance_cv_certificate(TRUE)
  expect_true(.nuisance_cv_certificate_enabled())
  expect_null(.set_nuisance_cv_certificate(NULL))
  expect_error(fit_initial_density_ratio(matrix(0, 20, 1), rep(1, 20), c(1, 0),
    lambda = NA_real_, nuisance_cv_certificate = FALSE), "finite numeric")
  expect_true(.nuisance_cv_certificate_enabled())
  for (invalid in list(1, "TRUE", NA, c(TRUE, FALSE))) {
    expect_error(.set_nuisance_cv_certificate(invalid), "TRUE or FALSE")
    expect_true(.nuisance_cv_certificate_enabled())
  }
  .restore_nuisance_cv_certificate(inherited)
  expect_false(.nuisance_cv_certificate_enabled())
  expect_null(getOption("RoCE.nuisance_cv_certificate"))
})

test_that("public density fitters preserve lambda and coefficients with certified failures", {
  x <- matrix(seq(.2, 1, length.out = 60), ncol = 1L)
  a <- rep(1L, nrow(x))
  for (solver in c("proximal_newton", "coordinate_descent")) {
    for (type in c("initial", "refined", "calibrated")) {
      fit <- function(enabled) {
        set.seed(6011)
        if (type == "initial") {
          fit_initial_density_ratio(x, a, c(1, -.5), nlambda = 12L, M_tau = Inf,
            nuisance_solver = solver, nuisance_cv_certificate = enabled, tol = 1e-10)
        } else {
          fit_unified_density_ratio(x, a, c(.25, -.125), c(0, .3), W_outcome = x,
            nlambda = 12L, M_tau = 12, calibrated = type == "calibrated",
            nuisance_solver = solver, nuisance_cv_certificate = enabled, tol = 1e-10,
            truncate_initial_outcome = FALSE)
        }
      }
      reference <- fit(FALSE)
      accelerated <- fit(TRUE)
      expect_equal(as.numeric(accelerated), as.numeric(reference), tolerance = 0)
      for (field in c("lambda_used", "lambda_min", "lambda_1se", "cv_invalid_fold_fits",
                      "cv_invalid_lambdas", "cv_path_tail_skipped_fold_fits")) {
        expect_identical(attr(accelerated, field), attr(reference, field))
      }
      expect_gt(attr(accelerated, "cv_certified_nonconvergent_fold_fits"), 0L)
      expect_null(attr(reference, "cv_certified_nonconvergent_fold_fits"))
      expect_true(isTRUE(attr(accelerated, "converged")))
    }
  }
})

test_that("certificate cache provenance changes without changing fold seeds or caller RNG", {
  skip_if(Sys.getenv("ROCE_TEST_INSTALLED") != "1", "Persistent cache uses an installed package")
  directory <- tempfile("certificate_cache_")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  prior <- options(RoCE.nuisance_cv_certificate = FALSE)
  on.exit(options(prior), add = TRUE)
  arguments <- list(Z_site = matrix(seq(.2, 1, length.out = 60), ncol = 1L),
    A = rep(1L, 60), mean_phi = c(1, -.5), A_val = 1L, nlambda = 12L, M_tau = Inf)
  cache <- .new_nuisance_cache(TRUE, checkpoint_dir = directory)
  set.seed(9811)
  rng <- .Random.seed
  reference <- .fit_nuisance_training_subset("fit_initial_density_ratio", arguments, "s1", 2:3, cache)
  options(RoCE.nuisance_cv_certificate = TRUE)
  accelerated <- .fit_nuisance_training_subset("fit_initial_density_ratio", arguments, "s1", 2:3, cache)
  cached <- .fit_nuisance_training_subset("fit_initial_density_ratio", arguments, "s1", 2:3, cache)
  expect_identical(.Random.seed, rng)
  expect_identical(attr(accelerated, "cv_seed"), attr(reference, "cv_seed"))
  expect_equal(as.numeric(accelerated), as.numeric(reference), tolerance = 0)
  expect_gt(attr(accelerated, "cv_certified_nonconvergent_fold_fits"), 0L)
  expect_identical(attr(cached, "cv_certified_nonconvergent_fold_fits"),
                   attr(accelerated, "cv_certified_nonconvergent_fold_fits"))
  expect_length(list.files(directory, pattern = "[.]rds$", recursive = TRUE), 2L)
})

test_that("complete TATE fitting retains the certificate control across arms and aggregation", {
  skip_on_cran()
  set.seed(6731)
  data <- split_data_by_site(generate_bounded_data(n_target = 600L,
    n_source_sizes = 600L, p = 4L, config = "C2"))
  arguments <- list(data_split = data, n_folds = 4L, communication_mode = "one_round",
    crossfit_layers = 3L, nlambda_init = 12L, M_tau = 12, M_tau_inference = 12,
    target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
    calibration_layout = "compact", nuisance_solver = "proximal_newton",
    calibration_control = list(recipe = "score_derivative",
      target_propensity_initialization = "calibrated", target_radius = 12), verbose = FALSE)
  reference <- do.call(run_tate_crossfit, c(arguments, list(n_cores = 1L, nuisance_cv_certificate = FALSE)))
  accelerated <- do.call(run_tate_crossfit, c(arguments, list(n_cores = 2L,
    parallel_arms = TRUE, nuisance_cv_certificate = TRUE)))
  expect_identical(accelerated$nuisance_cv_certificate, TRUE)
  for (arm in c("mu1", "mu0")) {
    expect_true(accelerated$arm_results[[arm]]$nuisance_cv_certificate)
  }
  for (mode in c("common_tate", "separate_arms", "joint_tate")) {
    a <- reaggregate_tate_crossfit(data, reference, aggregation_mode = mode, M_tau_inference = 12)
    b <- reaggregate_tate_crossfit(data, accelerated, aggregation_mode = mode, M_tau_inference = 12)
    expect_equal(a$estimate, b$estimate, tolerance = 1e-12)
    expect_equal(a$se, b$se, tolerance = 1e-12)
    expect_equal(a$fold_weights, b$fold_weights, tolerance = 1e-12)
    expect_true(b$nuisance_cv_certificate)
  }
})
