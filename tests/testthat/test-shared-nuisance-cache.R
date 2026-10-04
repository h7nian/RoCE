test_that("shared nuisance cache verifies inputs across separate worker environments", {
  cache <- .new_nuisance_cache(TRUE, tempdir())
  path <- attr(cache, "shared_directory")
  on.exit(.release_nuisance_cache(cache), add = TRUE)
  set.seed(1171)
  x <- matrix(stats::rnorm(1000 * 4L), 1000, 4L)
  arguments <- list(W_outcome = x, Y = stats::rbinom(1000, 1L, .5), A = rep(1L, 1000),
                    A_val = 1L, nlambda = 8L)
  fit <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 3:10, cache)
  fresh_worker_cache <- function() {
    worker <- new.env(parent = emptyenv())
    attr(worker, "shared_directory") <- path
    worker
  }
  testthat::local_mocked_bindings(fit_initial_outcome = function(...) {
    stop("fitter was called")
  }, .package = "RoCE")
  state <- .Random.seed
  reused <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 3:10,
                                         fresh_worker_cache())
  expect_identical(as.numeric(reused), as.numeric(fit))
  expect_identical(attr(reused, "lambda_used"), attr(fit, "lambda_used"))
  expect_identical(.Random.seed, state)
  expect_identical(attr(reused, "cv_seed"), attr(fit, "cv_seed"))
  changed <- arguments
  changed$W_outcome[1L, 1L] <- changed$W_outcome[1L, 1L] + 1e-12
  expect_error(.fit_nuisance_training_subset("fit_initial_outcome", changed, "s1", 3:10,
                                            fresh_worker_cache()), "fitter was called", fixed = TRUE)
  file <- list.files(path, pattern = "[.]rds$", full.names = TRUE)
  stopifnot(length(file) == 1L)
  corrupted <- readRDS(file)
  corrupted$fit[1L] <- corrupted$fit[1L] + 1
  saveRDS(corrupted, file)
  expect_error(.fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 3:10,
                                            fresh_worker_cache()), "integrity check", fixed = TRUE)
})

test_that("only the cache owner removes a temporary shared cache", {
  expect_error(.new_nuisance_cache(FALSE, tempdir()), "requires use_lambda_cache=TRUE", fixed = TRUE)
  expect_error(.new_nuisance_cache(TRUE, tempfile()), "existing writable", fixed = TRUE)
  cache <- .new_nuisance_cache(TRUE, tempdir())
  path <- attr(cache, "shared_directory")
  on.exit(.release_nuisance_cache(cache), add = TRUE)
  if (.Platform$OS.type != "windows") {
    child <- parallel::mcparallel({ .release_nuisance_cache(cache); dir.exists(path) })
    expect_true(parallel::mccollect(child)[[1L]])
  }
  .release_nuisance_cache(cache)
  expect_false(dir.exists(path))
})

test_that("shared workers preserve three-level TATE results and clean up their cache", {
  skip_on_cran()
  set.seed(1173)
  data <- split_data_by_site(generate_simulation_data(n_total = 3000L, n_target = 1000L,
    n_source_sizes = c(1000L, 1000L), K = 2L, p = 4L, config = "C3",
    dgp_type = "face", outcome_type = "binary", warn_ignored = FALSE))
  folds <- build_crossfit_folds(data, 4L)
  arguments <- list(data_split = data, precomputed_folds = folds, n_folds = 4L,
    communication_mode = "one_round", nlambda_init = 8L, nuisance_tol = 1e-10,
    target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
    calibration_control = list(recipe = "score_derivative",
      target_propensity_initialization = "calibrated", target_radius = log(9)),
    nuisance_solver = "proximal_newton", calibration_layout = "compact", verbose = FALSE)
  parent <- tempfile("cache_parent_")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  reference <- do.call(run_tate_crossfit, c(arguments, list(n_cores = 1L)))
  shared <- do.call(run_tate_crossfit, c(arguments,
    list(n_cores = 2L, parallel_arms = TRUE, nuisance_cache_dir = parent)))
  expect_equal(reference$estimate, shared$estimate, tolerance = 1e-10)
  expect_equal(reference$se, shared$se, tolerance = 1e-10)
  expect_equal(reference$fold_weights, shared$fold_weights, tolerance = 1e-10)
  expect_identical(shared$nuisance_cache_dir, parent)
  expect_length(list.files(parent, all.files = TRUE, no.. = TRUE), 0L)
  for (mode in c("separate_arms", "joint_tate")) {
    a <- reaggregate_tate_crossfit(data, reference, M_tau_inference = reference$M_tau_inference,
                                  aggregation_mode = mode)
    b <- reaggregate_tate_crossfit(data, shared, M_tau_inference = shared$M_tau_inference,
                                  aggregation_mode = mode)
    expect_equal(a$estimate, b$estimate, tolerance = 1e-10)
    expect_equal(a$se, b$se, tolerance = 1e-10)
    expect_equal(a$fold_weights, b$fold_weights, tolerance = 1e-10)
    expect_identical(b$nuisance_cache_dir, parent)
  }
})
