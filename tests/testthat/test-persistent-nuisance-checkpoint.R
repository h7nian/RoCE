test_that("persistent nuisance fits survive worker exit and reject different inputs", {
  directory <- tempfile("persistent_nuisance_")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  set.seed(8011)
  arguments <- list(W_outcome = matrix(rnorm(800), 200, 4),
    Y = rbinom(200, 1, .5), A = rep(1L, 200), A_val = 1L, nlambda = 4L)
  cache <- .new_nuisance_cache(TRUE, checkpoint_dir = directory)
  path <- attr(cache, "shared_directory")
  fitted <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 2:4, cache)
  .release_nuisance_cache(cache)
  expect_true(dir.exists(path))
  testthat::local_mocked_bindings(fit_initial_outcome = function(...) stop("cache miss"), .package = "RoCE")
  restarted <- .new_nuisance_cache(TRUE, checkpoint_dir = directory)
  state <- .Random.seed
  recovered <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 2:4, restarted)
  expect_identical(as.numeric(recovered), as.numeric(fitted))
  expect_identical(attr(recovered, "lambda_used"), attr(fitted, "lambda_used"))
  expect_identical(.Random.seed, state)
  arguments$W_outcome[1, 1] <- arguments$W_outcome[1, 1] + 1e-10
  expect_error(.fit_nuisance_training_subset("fit_initial_outcome", arguments, "s1", 2:4, restarted), "cache miss")
  expect_error(.new_nuisance_cache(FALSE, checkpoint_dir = directory), "requires use_lambda_cache")
  expect_error(.new_nuisance_cache(TRUE, tempdir(), checkpoint_dir = directory), "only one")
})

test_that("persistent fits preserve results after interruption and a cutoff change", {
  directory <- tempfile("restart_simulation_")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  arguments <- list(sim_id = 8012L, p = 4L, K = 1L, config = "C3",
    n_target = 1000L, n_source_sizes = 1000L, n_folds = 4L, nlambda_init = 4L,
    methods = c("one_round_crossfit", "target_only"), estimate_ate = TRUE,
    nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
    calibration_layout = "compact", M_tau = 12, M_tau_inference = 12,
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "calibrated", target_radius = 12),
    include_quadratic_bias_rule = FALSE, additional_aggregation_modes = c("separate_arms", "joint_tate"),
    n_cores_internal = 1L, verbose = FALSE)
  reference <- do.call(run_single_simulation, arguments)
  interrupt_run <- function() {
    original <- .fit_nuisance_training_subset
    calls <- 0L
    testthat::local_mocked_bindings(.fit_nuisance_training_subset = function(...) {
      result <- original(...)
      calls <<- calls + 1L
      if (calls == 8L) stop("simulated interruption after checkpoint")
      result
    }, .package = "RoCE")
    do.call(run_single_simulation, c(arguments, list(checkpoint_dir = directory)))
  }
  expect_error(interrupt_run(), "simulated interruption")
  expect_gt(length(list.files(directory, pattern = "[.]rds$", recursive = TRUE)), 0L)
  resumed <- do.call(run_single_simulation, c(arguments, list(checkpoint_dir = directory)))
  columns <- c("method", "estimate", "se", "bias", "coverage", "truth")
  expect_equal(resumed[columns], reference[columns], tolerance = 1e-12)
  resumed_parallel <- do.call(run_single_simulation,
    c(arguments, list(checkpoint_dir = directory, parallel_treatment_arms = TRUE)))
  expect_equal(resumed_parallel[columns], reference[columns], tolerance = 1e-12)
  # Aggregation does not enter nuisance fitting. Reusing cutoff-1 nuisance
  # checkpoints for cutoff 2 must reproduce a fresh cutoff-2 computation.
  arguments$aggregation_lambda <- 0.5
  cutoff2_reference <- do.call(run_single_simulation, arguments)
  cutoff2_restored <- do.call(run_single_simulation,
    c(arguments, list(checkpoint_dir = directory, parallel_treatment_arms = TRUE)))
  expect_equal(cutoff2_restored[columns], cutoff2_reference[columns], tolerance = 1e-12)
  expect_true(all(cutoff2_restored$aggregation_lambda == 0.5))
  expect_true(all(cutoff2_restored$aggregation_cutoff == 2))
})
