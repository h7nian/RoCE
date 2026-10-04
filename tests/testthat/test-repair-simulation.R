test_that("simulation compares aggregation modes without moving the fixed target benchmark", {
  skip_on_cran()
  arguments <- list(sim_id = 1044L, n_total = 2000L, n_target = 1000L,
    n_source_sizes = 1000L, K = 1L, p = 4L, config = "C1", dgp_type = "face",
    outcome_type = "continuous", n_folds = 4L, nlambda_init = 6L,
    methods = c("one_round_crossfit", "target_only"), estimate_ate = TRUE,
    include_quadratic_bias_rule = FALSE, verbose = FALSE, n_cores_internal = 1L)
  ordinary <- suppressWarnings(do.call(run_single_simulation, arguments))
  arguments$target_nuisance_method <- "hou_calibrated"
  arguments$source_validation_method <- "calibrated"
  arguments$calibration_layout <- "compact"
  arguments$additional_aggregation_modes <- c("separate_arms", "joint_tate")
  calibrated <- suppressWarnings(do.call(run_single_simulation, arguments))
  expect_true(all(c("one_round_crossfit_ate_separate_arms", "one_round_crossfit_ate_joint_tate") %in%
                    calibrated$method))
  expect_true(all(is.finite(calibrated$estimate)))
  expect_true(all(is.finite(calibrated$se)))
  alternative_rows <- calibrated$method %in% c(
    "one_round_crossfit_ate_separate_arms", "one_round_crossfit_ate_joint_tate")
  expect_true(all(calibrated$estimand_scope[alternative_rows] == "tate"))
  expect_equal(calibrated$truth[alternative_rows],
               rep(calibrated$truth[calibrated$method == "target_only_ate"], 2L))
  expect_true("target_anchor_ate" %in% calibrated$method)
  expect_false("target_anchor_ate" %in% ordinary$method)
  expect_equal(calibrated$estimate[calibrated$method == "target_only_ate"],
               ordinary$estimate[ordinary$method == "target_only_ate"], tolerance = 1e-12)
  expect_equal(calibrated$se[calibrated$method == "target_only_ate"],
               ordinary$se[ordinary$method == "target_only_ate"], tolerance = 1e-12)
  joint <- calibrated[calibrated$method == "one_round_crossfit_ate_joint_tate", ]
  expect_identical(joint$aggregation_mode, "joint_tate")
  expect_true(all(c("target_anchor_weight_mu1", "target_anchor_weight_mu0") %in% names(joint)))
  expect_identical(joint$source_validation_method, "calibrated")
  expect_true(all(c("calibration_layout", "compiled_nuisance_solver", "target_nuisance_method",
                    "source_validation_method", "aggregation_mode") %in%
                   .simulation_diagnostic_group_columns(calibrated)))
})

test_that("simulation options are not silently ignored", {
  expect_error(run_single_simulation(1, additional_aggregation_modes = "joint_tate"),
               "requires estimate_ate=TRUE", fixed = TRUE)
  expect_error(run_simulation_study(n_cores = 0), "positive integer CPU budget", fixed = TRUE)
})

test_that("study resumes preserve method configuration and respect the CPU budget", {
  observed <- new.env(parent = emptyenv())
  observed$checkpoint <- NULL
  observed$calls <- 0L
  testthat::local_mocked_bindings(
    run_single_simulation = function(sim_id, n_cores_internal, ...) {
      observed$calls <- observed$calls + 1L
      observed$cores <- n_cores_internal
      data.frame(sim_id = sim_id, method = "fixture", estimate = 0)
    },
    load_checkpoint = function(...) observed$checkpoint,
    save_checkpoint = function(state, ...) {
      observed$checkpoint <- state
      TRUE
    },
    check_preempt_signal = function(...) FALSE,
    cleanup_checkpoint = function(...) invisible(TRUE),
    .package = "RoCE"
  )
  arguments <- list(n_sims = 1L, n_total_vec = 3000L, K_vec = 2L,
    p_vec = 4L, configs = "C1", n_cores = 1L,
    checkpoint_config = list(file = tempfile(), preempt_signal = tempfile(), saved_signal = tempfile()))
  invisible(capture.output(first <- suppressWarnings(do.call(run_simulation_study, arguments))))
  expect_equal(observed$cores, 1L)
  expect_identical(observed$calls, 1L)
  expect_true(is.list(observed$checkpoint$run_configuration))
  arguments$n_total_vec <- 3000  # numeric storage must not change the method identity
  arguments$n_cores <- 2L       # a completed checkpoint does not need a new cluster
  invisible(capture.output(resumed <- suppressWarnings(do.call(run_simulation_study, arguments))))
  expect_identical(resumed, first)
  expect_identical(observed$calls, 1L)
  arguments$nuisance_tol <- 1e-9
  expect_error(suppressWarnings(do.call(run_simulation_study, arguments)),
               "Checkpoint method/configuration does not match", fixed = TRUE)
  expect_identical(observed$calls, 1L)
})
