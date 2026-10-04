test_that("layer control selects complete outer reuse or nested calibration", {
  control <- .validate_calibration_control()
  default <- c("initial", "calibrated")
  expect_identical(.resolve_crossfit_validation(default, NULL, control, "hou_calibrated"), "initial")
  expect_identical(.resolve_crossfit_validation(default, 2, control, "hou_calibrated"), "outer_fit")
  expect_identical(.resolve_crossfit_validation(default, 3, control, "hou_calibrated"), "calibrated")
  expect_identical(.resolve_crossfit_validation("complete", 2, control, "lasso"), "outer_fit")
  expect_error(.resolve_crossfit_validation(default, 4, control, "lasso"), "crossfit_layers")
  expect_error(.resolve_crossfit_validation("initial", 2, control, "lasso"), "legacy")
  expect_error(.resolve_crossfit_validation("outer_fit", 3, control, "lasso"), "two layers")
  expect_error(.resolve_crossfit_validation("typo", 2, control, "lasso"), "arg")
  standard <- .validate_calibration_control(list(source_nuisance_method = "standard"))
  expect_error(.resolve_crossfit_validation(default, 2, standard, "lasso"), "final calibration")
})

test_that("ordinary target training-score retention preserves outer predictions", {
  set.seed(4811)
  data <- split_data_by_site(generate_bounded_data(n_target = 600, n_source_sizes = 600,
    K = 1, p = 4, config = "C2"))
  folds <- build_crossfit_folds(data, 4L)$target_folds
  cache <- new.env(parent = emptyenv())
  ordinary <- .get_target_only_fold_fit(folds, 1, 4, "binomial", 1, cache)
  retained <- .get_target_only_fold_fit(folds, 1, 4, "binomial", 1, cache, return_training_scores = TRUE)
  for (field in c("estimate", "V_ot", "varphi_ot", "m_pred", "prop_scores")) {
    expect_equal(as.numeric(ordinary[[field]]), as.numeric(retained[[field]]), tolerance = 1e-12)
  }
  expected_rows <- unlist(lapply(folds[2:4], `[[`, "original_idx"), use.names = FALSE)
  expect_identical(retained$training_scores$original_idx, expected_rows)
  expect_identical(retained$training_scores$training_folds, 2:4)
  reused <- .reuse_outer_target_scores(retained, folds, 1, 2, 1, "binomial", 12)
  expect_true(reused$reused_outer_fit)
  expect_equal(reused$outcome_degenerate, 0L)
  expect_error(.reuse_outer_target_scores(retained, folds, 1, 1, 1, "binomial", 12), "exclude")
})

test_that("two layers retain the same final calibrated models in both communication modes", {
  skip_on_cran()
  source(file.path(repo_root, "diagnosis", "repair", "layer_comparison_summary.R"), local = TRUE)
  set.seed(4812)
  data <- split_data_by_site(generate_bounded_data(n_target = 600, n_source_sizes = 600,
    K = 1, p = 4, config = "C2"))
  run <- function(layers, protocol, target) {
    control <- list(recipe = "score_derivative")
    if (target == "hou_calibrated") {
      control$target_propensity_initialization <- "calibrated"
      control$target_radius <- 12
    }
    run_tate_crossfit(data, n_folds = 4L, communication_mode = protocol, nlambda_init = 5L,
      target_nuisance_method = target, crossfit_layers = layers, calibration_control = control,
      nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
      M_tau = 12, M_tau_inference = 12, lambda_selection = .5, aggregation_mode = "joint_tate", verbose = FALSE)
  }
  for (protocol in c("one_round", "two_round")) for (target in c("hou_calibrated", "lasso")) {
    two <- run(2L, protocol, target)
    three <- run(3L, protocol, target)
    expect_equal(two$crossfit_levels, 2L)
    expect_equal(three$crossfit_levels, 3L)
    expect_identical(two$source_validation_method, "outer_fit")
    expect_equal(two$target_only, three$target_only, tolerance = 1e-12)
    expect_equal(two$source_estimates, three$source_estimates, tolerance = 1e-12)
    expect_identical(layer_outer_fit_fingerprint(two), layer_outer_fit_fingerprint(three))
    expect_true(all(is.finite(c(two$estimate, two$se, three$estimate, three$se))))
    for (arm in c("mu1", "mu0")) for (outer in 1:4) {
      first <- two$arm_results[[arm]]$fold_results[[outer]]
      second <- three$arm_results[[arm]]$fold_results[[outer]]
      expect_null(.source_inner_fits(first$source_results$s1))
      expect_false(is.null(.source_inner_fits(second$source_results$s1)))
      for (parameter in c("gamma_s", "alpha_ts")) {
        expect_equal(as.numeric(first$source_results$s1[[parameter]]),
                     as.numeric(second$source_results$s1[[parameter]]), tolerance = 1e-12)
      }
      for (key in names(first$target_only_inner)) {
        expect_true(first$target_only_inner[[key]]$reused_outer_fit)
        expect_setequal(first$target_only_inner[[key]]$training_folds, setdiff(1:4, outer))
      }
      info <- two$arm_results[[arm]]$intermediates$inner_fold_info[[outer]]
      old_info <- three$arm_results[[arm]]$intermediates$inner_fold_info[[outer]]
      expect_equal(info$n_t, old_info$n_t)
      expect_equal(info$n_s, old_info$n_s)
      expect_false(any(unlist(lapply(info$components, `[[`, "target_idx")) %in%
                         two$arm_results[[arm]]$intermediates$fold_info[[outer]]$target_idx))
    }
  }
})

test_that("two-layer eta is unaffected by changes confined to its outer evaluation fold", {
  skip_on_cran()
  set.seed(4813)
  data <- split_data_by_site(generate_bounded_data(n_target = 600, n_source_sizes = 600,
    K = 1, p = 4, config = "C3"))
  folds <- build_crossfit_folds(data, 4L)
  changed <- data
  for (site in names(data)) {
    views <- if (site == "t") folds$target_folds else folds$source_folds[[site]]
    indices <- views[[1L]]$original_idx
    changed[[site]]$Y[indices] <- 1 - changed[[site]]$Y[indices]
    changed[[site]]$W_outcome[indices, 1L] <- changed[[site]]$W_outcome[indices, 1L] + .5
    changed[[site]]$Z_site[indices, 1L] <- changed[[site]]$Z_site[indices, 1L] + .5
  }
  fit <- function(input) run_tate_crossfit(input, n_folds = 4L, communication_mode = "one_round",
    target_nuisance_method = "hou_calibrated", crossfit_layers = 2L, nlambda_init = 5L,
    nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "calibrated", target_radius = 12),
    M_tau = 12, M_tau_inference = 12, lambda_selection = .5, aggregation_mode = "joint_tate", verbose = FALSE)
  original <- fit(data)
  perturbed <- fit(changed)
  expect_equal(original$fold_weights[1L, ], perturbed$fold_weights[1L, ], tolerance = 1e-12)
  for (arm in c("mu1", "mu0")) {
    first <- original$arm_results[[arm]]$intermediates$inner_fold_info[[1L]]
    second <- perturbed$arm_results[[arm]]$intermediates$inner_fold_info[[1L]]
    for (field in c("avg_target_est", "avg_source_est", "V_ot", "V_t", "V_s", "C_ot", "C_cross")) {
      expect_equal(first[[field]], second[[field]], tolerance = 1e-12)
    }
  }
})

test_that("two-layer simulation preserves checkpoint and source-refit equivalence", {
  skip_on_cran()
  directory <- tempfile("two_layer_checkpoint_")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  arguments <- list(sim_id = 4831L, p = 4L, K = 2L, config = "C3", dgp_type = "bounded",
    n_target = 1000L, n_source_sizes = c(1000L, 1000L), n_folds = 3L, nlambda_init = 4L,
    methods = c("one_round_crossfit", "target_only"), estimate_ate = TRUE, return_fitted_tate = TRUE,
    nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    target_nuisance_method = "hou_calibrated", crossfit_layers = 2L,
    calibration_layout = "compact", M_tau = 12, M_tau_inference = 12, aggregation_lambda = .5,
    calibration_control = list(recipe = "score_derivative", target_propensity_initialization = "calibrated", target_radius = 12),
    include_quadratic_bias_rule = FALSE, additional_aggregation_modes = c("separate_arms", "joint_tate"),
    n_cores_internal = 1L, verbose = FALSE)
  first <- do.call(run_single_simulation, c(arguments, list(checkpoint_dir = directory)))
  restored <- do.call(run_single_simulation, c(arguments,
    list(checkpoint_dir = directory, parallel_treatment_arms = TRUE)))
  columns <- c("method", "estimate", "se", "truth", "roce_crossfit_levels", "source_validation_method")
  expect_equal(restored[columns], first[columns], tolerance = 1e-12)
  expect_true(all(first$roce_crossfit_levels == 2L))
  artifacts <- attr(first, "roce_simulation_artifacts")
  expect_named(artifacts$arm_truth, c("mu1", "mu0"))
  expect_equal(unname(diff(rev(artifacts$arm_truth))), artifacts$tate_truth, tolerance = 1e-12)
  source(file.path(repo_root, "diagnosis", "repair", "layer_comparison_summary.R"), local = TRUE)
  fit <- artifacts$direct_tate_results$one_round_crossfit
  expect_silent(summarize_layer_arms(fit, artifacts$data_split))
  changed <- artifacts$data_split
  treated <- which(changed$s1$A == 1L)
  changed$s1$Y[treated] <- 1 - changed$s1$Y[treated]
  reused <- .reuse_one_round_tate_across_rho(artifacts$data_split, changed, fit, "s1", .5)
  fresh <- run_tate_crossfit(changed, n_folds = 3L, communication_mode = "one_round",
    target_nuisance_method = "hou_calibrated", crossfit_layers = 2L, nlambda_init = 4L,
    nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
    calibration_control = arguments$calibration_control,
    M_tau = 12, M_tau_inference = 12, lambda_selection = .5, verbose = FALSE)
  for (field in c("estimate", "variance", "fold_weights", "source_estimates", "crossfit_levels")) {
    expect_equal(reused[[field]], fresh[[field]], tolerance = 1e-12)
  }
  standard <- arguments
  standard$calibration_control$source_nuisance_method <- "standard"
  standard_fit <- do.call(run_single_simulation, standard)
  expect_true(all(standard_fit$roce_crossfit_levels == 2L))
  expect_true(all(standard_fit$source_validation_method == "outer_fit"))
})
