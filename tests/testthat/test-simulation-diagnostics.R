make_diagnostic_fixture <- function(n = 200L, method = "target_only_ate") {
  z <- stats::qnorm(0.975)
  bias <- seq(-0.02, 0.02, length.out = n)
  se <- rep(stats::sd(bias), n)
  truth <- rep(0, n)
  estimate <- truth + bias
  data.frame(
    sim_id = seq_len(n),
    experiment = "test",
    outcome_family = "gaussian",
    config = "C1",
    p = 100L,
    K = 2L,
    rho = 0,
    cutoff = 2,
    n_site = 1000L,
    n_folds = 5L,
    estimand_scope = "tate",
    method = method,
    estimate = estimate,
    truth = truth,
    bias = bias,
    se = se,
    ci_lower = estimate - z * se,
    ci_upper = estimate + z * se,
    coverage = truth >= estimate - z * se & truth <= estimate + z * se,
    stringsAsFactors = FALSE
  )
}

test_that("TATE method classification includes diagnostic aggregation rules", {
  expect_true(all(.is_tate_method(c(
    "one_round_crossfit_ate",
    "one_round_crossfit_ate_armwise",
    "one_round_crossfit_ate_hard_threshold",
    "one_round_crossfit_ate_quadratic_bias"
  ))))
  expect_false(any(.is_tate_method(c(
    "one_round_crossfit", "target_only", "treated_mean"
  ))))
})

test_that("weight-bootstrap controls and seeds are deterministic", {
  expect_equal(RoCE:::.validate_weight_bootstrap_replicates(0), 0L)
  expect_equal(RoCE:::.validate_weight_bootstrap_replicates(20), 20L)
  expect_error(
    RoCE:::.validate_weight_bootstrap_replicates(1),
    "0 \\(disabled\\) or at least 2"
  )
  expect_error(
    RoCE:::.validate_weight_bootstrap_replicates(2.5),
    "non-negative integer"
  )
  soft_seed <- RoCE:::.tate_weight_bootstrap_seed(17L, "soft_penalty")
  expect_identical(
    soft_seed,
    RoCE:::.tate_weight_bootstrap_seed(17L, "soft_penalty")
  )
  expect_false(identical(
    soft_seed,
    RoCE:::.tate_weight_bootstrap_seed(17L, "hard_threshold")
  ))
  expect_error(
    RoCE:::.tate_weight_bootstrap_seed(0L), "positive integer"
  )
})

test_that("simulation weight-bootstrap adapter is opt-in and RNG-neutral", {
  capture <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    estimate_tate_weight_bootstrap = function(
        tate_result, B, seed, relearn_weights, screening_rule,
        save_draws) {
      capture$arguments <- list(
        tate_result = tate_result, B = B, seed = seed,
        relearn_weights = relearn_weights,
        screening_rule = screening_rule, save_draws = save_draws
      )
      list(
        variance_weight_relearn_bootstrap = 0.04,
        se_weight_relearn_bootstrap = 0.2,
        se_fixed_weight_bootstrap = 0.1,
        weight_uncertainty_ratio = 2,
        n_bootstrap = B,
        bootstrap_multiplier = "site_stratified_exponential",
        bootstrap_seed = seed,
        bootstrap_failures = 0L
      )
    },
    .package = "RoCE"
  )
  expect_identical(
    RoCE:::.run_tate_weight_bootstrap(list(), 0L, 7L),
    list()
  )
  set.seed(511L)
  expected_next <- stats::runif(1L)
  set.seed(511L)
  result <- RoCE:::.run_tate_weight_bootstrap(
    tate_result = list(marker = "fit"),
    n_weight_bootstrap = 20L,
    sim_id = 7L,
    screening_rule = "soft_penalty"
  )
  actual_next <- stats::runif(1L)
  expect_equal(actual_next, expected_next, tolerance = 0)
  expect_equal(capture$arguments$B, 20L)
  expect_equal(capture$arguments$seed, result$weight_bootstrap_seed)
  expect_true(capture$arguments$relearn_weights)
  expect_false(capture$arguments$save_draws)
  expect_equal(result$se_weight_relearn_bootstrap, 0.2)
  expect_equal(result$weight_uncertainty_ratio, 2)
})

test_that("weight-relearning bootstrap diagnostics fail closed", {
  fixture <- make_diagnostic_fixture(
    n = 3L, method = "one_round_crossfit_ate"
  )
  fixture$n_weight_bootstrap <- 20L
  fixture$variance_weight_relearn_bootstrap <- 0.04
  fixture$se_weight_relearn_bootstrap <- 0.2
  fixture$se_fixed_weight_bootstrap <- 0.1
  fixture$weight_uncertainty_ratio <- 2
  fixture$weight_relearn_n_bootstrap <- 20L
  fixture$weight_bootstrap_multiplier <- "exponential"
  fixture$weight_bootstrap_seed <- 1001:1003
  fixture$weight_bootstrap_failures <- 0L
  fixture$weight_bootstrap_relearn_weights <- TRUE
  fixture$weight_bootstrap_screening_rule <- "soft_penalty"

  valid <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )
  expect_false(grepl(
    "invalid_weight_relearn_bootstrap_diagnostics",
    valid$diagnostic_status
  ))
  expect_equal(valid$weight_relearn_bootstrap_replications, 3L)
  expect_equal(valid$mean_se_weight_relearn_bootstrap, 0.2)
  expect_equal(valid$mean_se_fixed_weight_bootstrap, 0.1)
  expect_equal(valid$mean_weight_uncertainty_ratio, 2)
  expect_equal(valid$weight_bootstrap_failures_total, 0)
  expect_equal(valid$weight_bootstrap_variance_identity_error, 0)
  expect_equal(valid$weight_bootstrap_ratio_identity_error, 0)

  two_round <- fixture
  two_round$method <- "two_round_crossfit_ate"
  two_round[c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
    "weight_bootstrap_failures", "weight_bootstrap_multiplier",
    "weight_bootstrap_relearn_weights",
    "weight_bootstrap_screening_rule"
  )] <- NA
  combined <- RoCE:::diagnose_simulation_results(
    rbind(fixture, two_round), expected_replications = 3L
  )
  two_round_diagnostic <- combined[
    combined$method == "two_round_crossfit_ate", , drop = FALSE
  ]
  expect_equal(nrow(two_round_diagnostic), 1L)
  expect_false(grepl(
    "invalid_weight_relearn_bootstrap_diagnostics",
    two_round_diagnostic$diagnostic_status
  ))
  expect_equal(
    two_round_diagnostic$weight_relearn_bootstrap_replications, 0L
  )

  invalid <- fixture
  invalid$weight_uncertainty_ratio[[1L]] <- 3
  diagnosed <- RoCE:::diagnose_simulation_results(
    invalid, expected_replications = 3L
  )
  expect_match(
    diagnosed$diagnostic_status,
    "invalid_weight_relearn_bootstrap_diagnostics"
  )

  disabled <- fixture
  disabled$n_weight_bootstrap <- 0L
  diagnosed <- RoCE:::diagnose_simulation_results(
    disabled, expected_replications = 3L
  )
  expect_match(
    diagnosed$diagnostic_status,
    "invalid_weight_relearn_bootstrap_diagnostics"
  )

  partial <- make_diagnostic_fixture(n = 3L)
  partial$n_weight_bootstrap <- 0L
  partial$se_weight_relearn_bootstrap <- NA_real_
  diagnosed <- RoCE:::diagnose_simulation_results(
    partial, expected_replications = 3L
  )
  expect_match(
    diagnosed$diagnostic_status,
    "invalid_weight_relearn_bootstrap_diagnostics"
  )

  invalid_logical <- fixture
  invalid_logical$weight_bootstrap_relearn_weights <- "TRUE"
  diagnosed <- RoCE:::diagnose_simulation_results(
    invalid_logical, expected_replications = 3L
  )
  expect_match(
    diagnosed$diagnostic_status,
    "invalid_weight_relearn_bootstrap_diagnostics"
  )

  false_logical <- fixture
  false_logical$weight_bootstrap_relearn_weights <- FALSE
  diagnosed <- RoCE:::diagnose_simulation_results(
    false_logical, expected_replications = 3L
  )
  expect_match(
    diagnosed$diagnostic_status,
    "invalid_weight_relearn_bootstrap_diagnostics"
  )
})

test_that("simulation diagnostics check complete replicate-level settings", {
  fixture <- make_diagnostic_fixture()
  diagnostic <- RoCE:::diagnose_simulation_results(fixture)

  expect_equal(nrow(diagnostic), 1L)
  expect_equal(diagnostic$n_unique_replications, 200L)
  expect_equal(diagnostic$n_complete_replications, 200L)
  expect_equal(diagnostic$rmse, sqrt(mean(fixture$bias^2)))
  expect_equal(diagnostic$rmse, diagnostic$rmse_from_moments, tolerance = 1e-12)
  expect_equal(diagnostic$rmse_identity_error, 0, tolerance = 1e-12)
  expect_false(grepl("incomplete_replications", diagnostic$diagnostic_status))
  expect_false(grepl("rmse_identity_failed", diagnostic$diagnostic_status))
  expect_equal(
    diagnostic$truth_below_interval_fraction +
      diagnostic$truth_above_interval_fraction,
    1 - diagnostic$coverage
  )
  expect_true(is.finite(diagnostic$bias_skewness))
  expect_true(is.finite(diagnostic$bias_excess_kurtosis))
})

test_that("binary cell diagnostics expose sparse site-arm-outcome cells", {
  make_site <- function(a, y) list(A = a, Y = y, n = length(y))
  split <- list(
    t = make_site(c(0, 0, 1, 1), c(0, 1, 0, 1)),
    s1 = make_site(c(0, 0, 0, 1), c(0, 0, 1, 1))
  )

  diagnostic <- RoCE:::.binary_cell_diagnostics(split, "binomial")

  expect_equal(diagnostic$min_site_arm_outcome_cell_n, 0L)
  expect_equal(diagnostic$min_target_arm_outcome_cell_n, 1L)
  expect_equal(diagnostic$n_site_arm_outcome_cells_below_8, 8L)
  expect_true(all(is.na(unlist(
    RoCE:::.binary_cell_diagnostics(split, "gaussian")
  ))))
  split$s1$Y[[1L]] <- 2
  expect_error(
    RoCE:::.binary_cell_diagnostics(split, "binomial"),
    "aligned binary A and Y"
  )
})

test_that("setting diagnostics summarize design sparsity and DR clipping", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "federated_dr_ate")
  fixture$min_site_arm_outcome_cell_n <- c(4L, 7L, 9L)
  fixture$min_target_arm_outcome_cell_n <- c(12L, 11L, 10L)
  fixture$n_site_arm_outcome_cells_below_8 <- c(2L, 1L, 0L)
  fixture$dr_weight_n <- 2000L
  fixture$dr_weight_n_clipped <- c(20L, 0L, 40L)
  fixture$dr_weight_fraction_clipped <- c(0.01, 0, 0.02)
  fixture$dr_weight_max_site_fraction_clipped <- c(0.02, 0, 0.04)
  fixture$dr_weight_min_before_clipping <- c(0.08, 0.11, 0.04)
  fixture$dr_weight_max_before_clipping <- c(11, 8, 14)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$min_site_arm_outcome_cell_n_observed, 4)
  expect_equal(diagnostic$min_target_arm_outcome_cell_n_observed, 10)
  expect_equal(diagnostic$sparse_binary_cell_replications, 2L)
  expect_equal(diagnostic$mean_dr_weight_fraction_clipped, 0.01)
  expect_equal(diagnostic$max_dr_weight_site_fraction_clipped, 0.04)
  expect_equal(diagnostic$dr_weight_clipped_replications, 2L)
  expect_equal(diagnostic$min_dr_weight_before_clipping, 0.04)
  expect_equal(diagnostic$max_dr_weight_before_clipping, 14)
  expect_match(diagnostic$diagnostic_status, "sparse_binary_cells_detected")
  expect_match(diagnostic$diagnostic_status, "density_ratio_clipping_detected")
})

test_that("setting diagnostics reject internally inconsistent design diagnostics", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "federated_dr_ate")
  fixture$min_site_arm_outcome_cell_n <- c(4, Inf, 9)
  fixture$min_target_arm_outcome_cell_n <- c(12, 11, 10)
  fixture$n_site_arm_outcome_cells_below_8 <- c(2, 1, 0)
  fixture$dr_weight_n <- 100L
  fixture$dr_weight_n_clipped <- c(2L, 101L, 0L)
  fixture$dr_weight_fraction_clipped <- c(0.02, 1.01, 0)
  fixture$dr_weight_max_site_fraction_clipped <- c(0.03, 1.01, 0)
  fixture$dr_weight_min_before_clipping <- c(0.1, 0.1, 0.1)
  fixture$dr_weight_max_before_clipping <- c(10, 10, 10)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_match(diagnostic$diagnostic_status, "invalid_cell_diagnostics")
  expect_match(diagnostic$diagnostic_status, "invalid_dr_weight_diagnostics")
})

test_that("binary and DR diagnostics fail closed when metadata are incomplete", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "federated_dr_ate")
  fixture$outcome_family <- "binomial"
  fixture$min_site_arm_outcome_cell_n <- c(4L, NA, 9L)
  fixture$min_target_arm_outcome_cell_n <- c(12L, 11L, 10L)
  fixture$n_site_arm_outcome_cells_below_8 <- c(2L, 1L, 0L)
  fixture$dr_weight_n <- 2000L
  fixture$dr_weight_n_clipped <- c(20L, 0L, 40L)
  fixture$dr_weight_fraction_clipped <- c(0.01, 0, NA)
  fixture$dr_weight_max_site_fraction_clipped <- c(0.02, 0, 0.04)
  fixture$dr_weight_min_before_clipping <- c(0.08, 0.11, 0.04)
  fixture$dr_weight_max_before_clipping <- c(11, 8, 14)
  fixture$n_bootstrap <- 5000L
  fixture$comparison_variance_method <- "bootstrap"
  fixture$comparison_se_analytic <- fixture$se
  fixture$comparison_se_bootstrap <- fixture$se
  fixture$comparison_n_bootstrap <- 5000L

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_match(diagnostic$diagnostic_status, "invalid_cell_diagnostics")
  expect_match(diagnostic$diagnostic_status, "invalid_dr_weight_diagnostics")
})

test_that("DR diagnostics fail closed when every DR metadata column is absent", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "federated_dr_ate")

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_match(diagnostic$diagnostic_status, "invalid_dr_weight_diagnostics")
  expect_true(is.na(diagnostic$dr_weight_clipped_replications))
})

test_that("comparison bootstrap diagnostics fail closed when absent", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "sample_size_ate")
  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )
  expect_match(diagnostic$diagnostic_status, "invalid_bootstrap_diagnostics")
})

test_that("simulation diagnostics reject ambiguous numeric metadata", {
  fixture <- make_diagnostic_fixture(n = 3L)

  expect_error(
    RoCE:::diagnose_simulation_results(
      fixture, expected_replications = 2.5
    ),
    "positive integer"
  )
  expect_error(
    RoCE:::diagnose_simulation_results(
      fixture, identity_tolerance = -1
    ),
    "non-negative"
  )
  fixture$estimate <- as.character(fixture$estimate)
  expect_error(
    RoCE:::diagnose_simulation_results(fixture),
    "nonnumeric storage"
  )
})

test_that("TATE benchmark labels are defined once and remain complete", {
  expect_identical(
    RoCE:::.tate_benchmark_methods(),
    c(
      "target_only_ate", "sample_size_ate", "inverse_variance_ate",
      "federated_dr_ate", "pooled_dr_ate"
    )
  )
})

test_that("TATE aggregation diagnostic schema is defined once", {
  columns <- RoCE:::.direct_tate_aggregation_diagnostic_columns()
  expect_true(all(c(
    "target_anchor_weight", "max_wald_statistic",
    "max_weight_optimizer_iterations", "inference_safety_clip_count"
  ) %in% columns))
  expect_identical(anyDuplicated(columns), 0L)
})

test_that("hard-threshold diagnostics enforce exclusion and zero weight", {
  fixture <- make_diagnostic_fixture(
    n = 2L, method = "one_round_crossfit_ate_hard_threshold"
  )
  defaults <- c(
    target_anchor_weight = 0.8, mean_abs_source_weight = 0.1,
    max_abs_source_weight = 0.2, max_wald_statistic = 2,
    mean_wald_statistic = 1, penalized_source_fold_fraction = 0.5,
    max_weight_optimizer_iterations = 3,
    max_weight_psd_ridge = 0, weight_psd_ridge_fold_fraction = 0,
    se_fixed_weights = 0.01, weight_layer_indirect_variance = 1e-6,
    weight_layer_cross_term = 0, weight_layer_kink_cells = 0,
    inference_logit_truncated = 0,
    inference_logit_truncation_fraction = 0,
    inference_max_abs_logit = 2, inference_safety_clip_count = 0
  )
  for (name in names(defaults)) fixture[[name]] <- defaults[[name]]
  fixture$source_s1_inclusion_fraction <- 0.5
  fixture$source_s1_fold_1_included <- 0
  fixture$source_s1_fold_1_weight <- 0
  fixture$source_s1_fold_2_included <- 1
  fixture$source_s1_fold_2_weight <- 0.2

  valid <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 2L
  )
  expect_false(grepl(
    "invalid_direct_tate_diagnostics", valid$diagnostic_status
  ))

  fixture$source_s1_fold_1_weight[[2L]] <- 0.01
  invalid <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 2L
  )
  expect_match(
    invalid$diagnostic_status, "invalid_direct_tate_diagnostics"
  )
})

test_that("production scientific metadata enforces the audited p=100 design", {
  fixture <- make_diagnostic_fixture(n = 3L)
  fixture$dgp_type <- "face"
  fixture$outcome_family <- "binomial"
  fixture$estimand_type <- "superpopulation"
  fixture$heterogeneity_type <- "none"

  expect_invisible(
    RoCE:::.validate_face_production_scientific_metadata(fixture)
  )
  fixture$rho <- 1
  expect_error(
    RoCE:::.validate_face_production_scientific_metadata(fixture),
    "rho-consistent"
  )
  fixture$heterogeneity_type <- "one_deviated_source"
  expect_invisible(
    RoCE:::.validate_face_production_scientific_metadata(fixture)
  )
  fixture$p <- 50L
  expect_error(
    RoCE:::.validate_face_production_scientific_metadata(fixture),
    "audited p=100"
  )
  fixture$p <- 100L
  fixture$estimand_scope <- "treated_mean"
  expect_error(
    RoCE:::.validate_face_production_scientific_metadata(fixture),
    "method-consistent"
  )
  fixture$estimand_scope <- "tate"
  fixture$rho <- 0.25
  fixture$heterogeneity_type <- "one_deviated_source"
  expect_error(
    RoCE:::.validate_face_production_scientific_metadata(fixture),
    "audited p=100"
  )
})

test_that("simulation diagnostics never merge nuisance lambda grids", {
  grid_50 <- make_diagnostic_fixture(n = 3L)
  grid_50$nlambda_init <- 50L
  grid_100 <- grid_50
  grid_100$nlambda_init <- 100L

  diagnostic <- RoCE:::diagnose_simulation_results(rbind(grid_50, grid_100))

  expect_equal(nrow(diagnostic), 2L)
  expect_equal(sort(diagnostic$nlambda_init), c(50L, 100L))
  expect_true(all(diagnostic$n_unique_replications == 3L))
})

test_that("simulation diagnostics never merge nuisance CV selection rules", {
  minimum <- make_diagnostic_fixture(n = 3L)
  minimum$nuisance_lambda_rule <- "min"
  one_se <- minimum
  one_se$nuisance_lambda_rule <- "1se"

  diagnostic <- RoCE:::diagnose_simulation_results(rbind(minimum, one_se))

  expect_equal(nrow(diagnostic), 2L)
  expect_setequal(diagnostic$nuisance_lambda_rule, c("min", "1se"))
  expect_true(all(diagnostic$n_unique_replications == 3L))
})

test_that("simulation diagnostics never merge scientific DGP settings", {
  compatible <- make_diagnostic_fixture(n = 3L)
  compatible$dgp_type <- "face"
  compatible$heterogeneity_type <- "none"
  compatible$estimand_type <- "superpopulation"
  deviated <- compatible
  deviated$heterogeneity_type <- "one_deviated_source"

  diagnostic <- RoCE:::diagnose_simulation_results(
    rbind(compatible, deviated)
  )

  expect_equal(nrow(diagnostic), 2L)
  expect_setequal(
    diagnostic$heterogeneity_type,
    c("none", "one_deviated_source")
  )
})

test_that("simulation diagnostics retain nuisance convergence review signals", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "one_round_crossfit_ate")
  fixture$face_initial_dr_nonconverged <- c(0, 1, 0)
  fixture$face_calibrated_dr_nonconverged <- c(2, 0, 0)
  fixture$face_calibrated_outcome_nonconverged <- c(0, 0, 0)
  fixture$face_initial_dr_line_search_failures <- c(0, 3, 0)
  fixture$face_calibrated_dr_line_search_failures <- c(2, 0, 0)
  fixture$face_calibrated_outcome_line_search_failures <- c(0, 0, 4)
  fixture$face_initial_dr_cv_invalid_fold_fits <- c(0, 4, 0)
  fixture$face_initial_dr_cv_invalid_lambdas <- c(0, 1, 0)
  fixture$face_initial_dr_cv_path_tail_skipped_fold_fits <- c(0, 2, 0)
  fixture$face_calibrated_dr_cv_invalid_fold_fits <- c(2, 0, 0)
  fixture$face_calibrated_dr_cv_invalid_lambdas <- c(1, 0, 0)
  fixture$face_calibrated_dr_cv_path_tail_skipped_fold_fits <- c(1, 0, 0)
  fixture$face_calibrated_outcome_cv_invalid_fold_fits <- c(0, 0, 3)
  fixture$face_calibrated_outcome_cv_invalid_lambdas <- c(0, 0, 1)
  fixture$face_calibrated_outcome_cv_path_tail_skipped_fold_fits <- c(0, 0, 2)
  fixture$face_mu1_initial_dr_nonconverged <- c(0, 1, 0)
  fixture$face_mu0_initial_dr_nonconverged <- c(0, 0, 0)
  fixture$face_mu1_calibrated_dr_nonconverged <- c(1, 0, 0)
  fixture$face_mu0_calibrated_dr_nonconverged <- c(1, 0, 0)
  fixture$face_mu1_calibrated_outcome_nonconverged <- c(0, 0, 0)
  fixture$face_mu0_calibrated_outcome_nonconverged <- c(0, 0, 0)
  fixture$face_mu1_initial_dr_line_search_failures <- c(0, 3, 0)
  fixture$face_mu0_initial_dr_line_search_failures <- c(0, 0, 0)
  fixture$face_mu1_calibrated_dr_line_search_failures <- c(2, 0, 0)
  fixture$face_mu0_calibrated_dr_line_search_failures <- c(0, 0, 0)
  fixture$face_mu1_calibrated_outcome_line_search_failures <- c(0, 0, 4)
  fixture$face_mu0_calibrated_outcome_line_search_failures <- c(0, 0, 0)
  fixture$face_initial_dr_support_floor_applied <- c(1, 0, 2)
  fixture$face_calibrated_dr_support_floor_applied <- c(0, 1, 0)
  fixture$face_mu1_initial_dr_support_floor_applied <- c(1, 0, 0)
  fixture$face_mu0_initial_dr_support_floor_applied <- c(0, 0, 2)
  fixture$face_mu1_calibrated_dr_support_floor_applied <- c(0, 1, 0)
  fixture$face_mu0_calibrated_dr_support_floor_applied <- c(0, 0, 0)
  fixture$face_max_initial_dr_iterations <- c(20, 10000, 30)
  fixture$face_max_calibrated_dr_iterations <- c(10000, 50, 40)
  fixture$face_max_calibrated_outcome_iterations <- c(10, 12, 11)
  fixture$face_max_initial_dr_update_ratio <- c(0.9, 1.8, 0.7)
  fixture$face_max_initial_dr_abs_coefficient <- c(1.0, 1.7, 1.2)
  fixture$face_max_calibrated_dr_update_ratio <- c(2.4, 0.8, 0.6)
  fixture$face_max_calibrated_dr_abs_coefficient <- c(1.2, 1.5, 1.1)
  fixture$face_max_calibrated_outcome_update_ratio <- c(0.7, 1.3, 0.9)
  fixture$face_max_calibrated_outcome_abs_coefficient <- c(1.9, 1.4, 2.1)
  fixture$face_max_initial_dr_support_floor <- c(0.1, 0, 0.3)
  fixture$face_max_calibrated_dr_support_floor <- c(0, 0.2, 0)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$nuisance_diagnostic_replications, 3L)
  expect_equal(diagnostic$nuisance_nonconverged_replications, 2L)
  expect_equal(diagnostic$nuisance_nonconverged_fits_total, 3)
  expect_equal(diagnostic$face_initial_dr_nonconverged_total, 1)
  expect_equal(diagnostic$face_calibrated_dr_nonconverged_total, 2)
  expect_equal(diagnostic$face_calibrated_outcome_nonconverged_total, 0)
  expect_equal(diagnostic$face_initial_dr_cv_invalid_fold_fits_total, 4)
  expect_equal(diagnostic$face_initial_dr_cv_invalid_lambdas_total, 1)
  expect_equal(
    diagnostic$face_initial_dr_cv_path_tail_skipped_fold_fits_total, 2
  )
  expect_equal(
    diagnostic$face_calibrated_dr_cv_invalid_fold_fits_total, 2
  )
  expect_equal(diagnostic$face_calibrated_dr_cv_invalid_lambdas_total, 1)
  expect_equal(
    diagnostic$face_calibrated_dr_cv_path_tail_skipped_fold_fits_total, 1
  )
  expect_equal(
    diagnostic$face_calibrated_outcome_cv_invalid_fold_fits_total, 3
  )
  expect_equal(
    diagnostic$face_calibrated_outcome_cv_invalid_lambdas_total, 1
  )
  expect_equal(
    diagnostic$face_calibrated_outcome_cv_path_tail_skipped_fold_fits_total,
    2
  )
  expect_equal(diagnostic$face_initial_dr_line_search_failures_total, 3)
  expect_equal(diagnostic$face_calibrated_dr_line_search_failures_total, 2)
  expect_equal(
    diagnostic$face_calibrated_outcome_line_search_failures_total, 4
  )
  expect_equal(diagnostic$face_initial_dr_support_floor_applied_total, 3)
  expect_equal(diagnostic$face_calibrated_dr_support_floor_applied_total, 1)
  expect_equal(diagnostic$face_mu1_initial_dr_nonconverged_total, 1)
  expect_equal(diagnostic$face_mu0_initial_dr_nonconverged_total, 0)
  expect_equal(diagnostic$face_mu1_calibrated_dr_nonconverged_total, 1)
  expect_equal(diagnostic$face_mu0_calibrated_dr_nonconverged_total, 1)
  expect_equal(diagnostic$face_mu1_initial_dr_line_search_failures_total, 3)
  expect_equal(diagnostic$face_mu0_initial_dr_line_search_failures_total, 0)
  expect_equal(
    diagnostic$face_mu1_calibrated_dr_line_search_failures_total, 2
  )
  expect_equal(
    diagnostic$face_mu0_calibrated_dr_line_search_failures_total, 0
  )
  expect_equal(
    diagnostic$face_mu1_calibrated_outcome_line_search_failures_total, 4
  )
  expect_equal(
    diagnostic$face_mu0_calibrated_outcome_line_search_failures_total, 0
  )
  expect_equal(
    diagnostic$face_mu1_initial_dr_support_floor_applied_total, 1
  )
  expect_equal(
    diagnostic$face_mu0_initial_dr_support_floor_applied_total, 2
  )
  expect_equal(
    diagnostic$face_mu1_calibrated_dr_support_floor_applied_total, 1
  )
  expect_equal(
    diagnostic$face_mu0_calibrated_dr_support_floor_applied_total, 0
  )
  expect_equal(diagnostic$face_max_initial_dr_iterations_observed, 10000)
  expect_equal(
    diagnostic$face_max_calibrated_dr_iterations_observed, 10000
  )
  expect_equal(
    diagnostic$face_max_calibrated_outcome_iterations_observed, 12
  )
  expect_equal(diagnostic$face_max_initial_dr_update_ratio_observed, 1.8)
  expect_equal(diagnostic$face_max_initial_dr_abs_coefficient_observed, 1.7)
  expect_equal(diagnostic$face_max_calibrated_dr_update_ratio_observed, 2.4)
  expect_equal(
    diagnostic$face_max_calibrated_dr_abs_coefficient_observed, 1.5
  )
  expect_equal(
    diagnostic$face_max_calibrated_outcome_update_ratio_observed, 1.3
  )
  expect_equal(
    diagnostic$face_max_calibrated_outcome_abs_coefficient_observed, 2.1
  )
  expect_equal(diagnostic$face_max_initial_dr_support_floor_observed, 0.3)
  expect_equal(diagnostic$face_max_calibrated_dr_support_floor_observed, 0.2)
  expect_match(
    diagnostic$diagnostic_status, "nuisance_nonconvergence_detected"
  )
  expect_match(
    diagnostic$diagnostic_status, "nuisance_line_search_failures_detected"
  )
  expect_match(
    diagnostic$diagnostic_status, "nuisance_cv_candidates_excluded"
  )
})

test_that("support-floor activation is reported but is not nonconvergence", {
  fixture <- make_diagnostic_fixture(n = 2L, method = "one_round_crossfit_ate")
  fixture$face_initial_dr_nonconverged <- 0
  fixture$face_calibrated_dr_nonconverged <- 0
  fixture$face_calibrated_outcome_nonconverged <- 0
  fixture$face_initial_dr_support_floor_applied <- c(1, 0)
  fixture$face_calibrated_dr_support_floor_applied <- c(0, 1)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 2L
  )

  expect_equal(diagnostic$nuisance_nonconverged_replications, 0)
  expect_equal(diagnostic$nuisance_nonconverged_fits_total, 0)
  expect_equal(diagnostic$face_initial_dr_support_floor_applied_total, 1)
  expect_equal(diagnostic$face_calibrated_dr_support_floor_applied_total, 1)
  expect_false(grepl(
    "nuisance_nonconvergence_detected", diagnostic$diagnostic_status
  ))
})

test_that("nuisance count diagnostics reject fractional values", {
  fixture <- make_diagnostic_fixture(
    n = 2L, method = "one_round_crossfit_ate"
  )
  fixture$face_initial_dr_nonconverged <- c(0, 0.5)
  fixture$face_max_initial_dr_iterations <- c(20, 20.5)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 2L
  )

  expect_match(diagnostic$diagnostic_status, "invalid_nuisance_diagnostics")
})

test_that("TATE nuisance diagnostics fail closed when absent", {
  fixture <- make_diagnostic_fixture(
    n = 2L, method = "one_round_crossfit_ate"
  )

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 2L
  )

  expect_match(diagnostic$diagnostic_status, "invalid_nuisance_diagnostics")
  expect_match(diagnostic$diagnostic_status, "invalid_direct_tate_diagnostics")
})

test_that("simulation diagnostics flag extreme unconstrained aggregation weights", {
  fixture <- make_diagnostic_fixture(
    n = 3L, method = "one_round_crossfit_ate"
  )
  fixture$max_abs_source_weight <- c(0.5, 12, 2)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$max_abs_source_weight_observed, 12)
  expect_match(diagnostic$diagnostic_status, "extreme_aggregation_weight")
})

test_that("simulation diagnostics retain aggregation-optimizer safeguards", {
  fixture <- make_diagnostic_fixture(
    n = 3L, method = "one_round_crossfit_ate"
  )
  fixture$max_weight_optimizer_iterations <- c(4, 7, 5)
  fixture$max_weight_psd_ridge <- c(0, 2e-8, 0)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$max_weight_optimizer_iterations_observed, 7)
  expect_equal(diagnostic$max_weight_psd_ridge_observed, 2e-8)
  expect_equal(diagnostic$weight_psd_ridge_replications, 1L)
  expect_false(grepl(
    "invalid_weight_optimizer_diagnostics",
    diagnostic$diagnostic_status
  ))

  fixture$max_weight_psd_ridge[[1L]] <- -1
  invalid <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )
  expect_match(
    invalid$diagnostic_status,
    "invalid_weight_optimizer_diagnostics"
  )
})

test_that("simulation diagnostics retain inference truncation signals", {
  fixture <- make_diagnostic_fixture(
    n = 3L, method = "one_round_crossfit_ate"
  )
  fixture$inference_logit_truncation_fraction <- c(0, 0.03, 0.06)
  fixture$inference_max_abs_logit <- c(4.8, 5.5, 7.2)
  fixture$inference_safety_clip_count <- c(0, 0, 2)

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(
    diagnostic$mean_inference_logit_truncation_fraction, 0.03
  )
  expect_equal(diagnostic$max_inference_abs_logit_observed, 7.2)
  expect_equal(diagnostic$inference_safety_clip_total, 2)
  expect_match(
    diagnostic$diagnostic_status, "inference_safety_clipping_detected"
  )
})

test_that("simulation diagnostics retain and flag incomplete non-finite rows", {
  fixture <- make_diagnostic_fixture()
  fixture$estimate[[1L]] <- NA_real_
  fixture$se[[2L]] <- -1
  fixture$coverage[[3L]] <- NA
  fixture <- fixture[-4L, ]
  diagnostic <- RoCE:::diagnose_simulation_results(fixture)

  expect_equal(diagnostic$n_rows, 199L)
  expect_equal(diagnostic$n_unique_replications, 199L)
  expect_equal(diagnostic$n_complete_replications, 197L)
  expect_match(diagnostic$diagnostic_status, "incomplete_replications")
  expect_match(diagnostic$diagnostic_status, "nonfinite_values")
  expect_match(diagnostic$diagnostic_status, "nonpositive_standard_error")
})

test_that("simulation diagnostics reject nonbinary numeric coverage codes", {
  fixture <- make_diagnostic_fixture(n = 3L)
  fixture$coverage <- as.numeric(fixture$coverage)
  fixture$coverage[[1L]] <- 2

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$n_complete_replications, 2L)
  expect_match(diagnostic$diagnostic_status, "nonfinite_values")
})

test_that("simulation diagnostics audit comparison bootstrap metadata", {
  fixture <- make_diagnostic_fixture(n = 3L, method = "sample_size_ate")
  fixture$n_bootstrap <- 5000L
  fixture$comparison_variance_method <- "bootstrap"
  fixture$comparison_se_analytic <- fixture$se * 0.98
  fixture$comparison_se_bootstrap <- fixture$se
  fixture$comparison_n_bootstrap <- 5000L

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$mean_comparison_se_bootstrap, mean(fixture$se))
  expect_equal(
    diagnostic$bootstrap_to_analytic_se,
    mean(fixture$se) / mean(fixture$se * 0.98)
  )
  expect_equal(diagnostic$bootstrap_report_identity_error, 0)
  expect_false(grepl(
    "invalid_bootstrap_diagnostics", diagnostic$diagnostic_status
  ))

  fixture$comparison_n_bootstrap[[1L]] <- 1L
  invalid <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )
  expect_match(invalid$diagnostic_status, "invalid_bootstrap_diagnostics")
})

test_that("simulation diagnostics detect duplicated replications and CI mismatch", {
  fixture <- make_diagnostic_fixture()
  fixture <- rbind(fixture, fixture[1L, ])
  fixture$coverage[[1L]] <- !fixture$coverage[[1L]]
  diagnostic <- RoCE:::diagnose_simulation_results(fixture)

  expect_equal(diagnostic$n_duplicate_replications, 1L)
  expect_match(diagnostic$diagnostic_status, "duplicated_replications")
  expect_match(diagnostic$diagnostic_status, "coverage_ci_inconsistent")
})

test_that("simulation diagnostics reject inconsistent statistical identities", {
  fixture <- make_diagnostic_fixture(n = 4L)
  fixture$bias[[1L]] <- fixture$bias[[1L]] + 0.01
  fixture$ci_upper[[2L]] <- fixture$ci_upper[[2L]] + 0.01
  fixture$ci_lower[[3L]] <- fixture$ci_upper[[3L]] + 0.01
  fixture$sim_id[[4L]] <- NA_real_

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 4L
  )

  expect_equal(diagnostic$n_invalid_replication_ids, 1L)
  expect_equal(diagnostic$n_complete_replications, 2L)
  expect_gt(diagnostic$bias_identity_error, 0)
  expect_gt(diagnostic$ci_arithmetic_error, 0)
  expect_match(diagnostic$diagnostic_status, "invalid_replication_id")
  expect_match(diagnostic$diagnostic_status, "invalid_confidence_interval")
  expect_match(diagnostic$diagnostic_status, "bias_truth_inconsistent")
  expect_match(diagnostic$diagnostic_status, "ci_arithmetic_inconsistent")
})

test_that("simulation diagnostics require finite truth and CI endpoints", {
  fixture <- make_diagnostic_fixture(n = 3L)
  fixture$truth[[1L]] <- NA_real_
  fixture$ci_lower[[2L]] <- Inf

  diagnostic <- RoCE:::diagnose_simulation_results(
    fixture, expected_replications = 3L
  )

  expect_equal(diagnostic$n_complete_replications, 1L)
  expect_match(diagnostic$diagnostic_status, "nonfinite_values")
  expect_match(diagnostic$diagnostic_status, "invalid_confidence_interval")
})

test_that("RMSE comparisons are computed within each setting", {
  target <- make_diagnostic_fixture(method = "target_only_ate")
  face <- make_diagnostic_fixture(method = "one_round_crossfit_ate")
  face$bias <- face$bias / 2
  face$estimate <- face$truth + face$bias
  z <- stats::qnorm(0.975)
  face$ci_lower <- face$estimate - z * face$se
  face$ci_upper <- face$estimate + z * face$se
  face$coverage <- face$truth >= face$ci_lower & face$truth <= face$ci_upper

  diagnostic <- RoCE:::diagnose_simulation_results(rbind(target, face))
  diagnostic <- RoCE:::add_simulation_rmse_comparisons(diagnostic)

  face_row <- diagnostic[diagnostic$method == "one_round_crossfit_ate", ]
  target_row <- diagnostic[diagnostic$method == "target_only_ate", ]
  expect_equal(face_row$rmse_rank, 1)
  expect_equal(target_row$rmse_rank, 2)
  expect_equal(face_row$rmse_relative_to_target, 0.5)
  expect_equal(target_row$rmse_relative_to_target, 1)
})

test_that("RMSE rankings never mix potential-outcome means with TATE", {
  target_tate <- make_diagnostic_fixture(method = "target_only_ate")
  face_tate <- make_diagnostic_fixture(method = "one_round_crossfit_ate")
  target_mean <- make_diagnostic_fixture(method = "target_only")
  target_mean$estimand_scope <- "treated_mean"
  treated_mean <- make_diagnostic_fixture(method = "sample_size")
  treated_mean$estimand_scope <- "treated_mean"
  treated_mean$bias <- treated_mean$bias / 100
  treated_mean$estimate <- treated_mean$truth + treated_mean$bias

  diagnostic <- RoCE:::diagnose_simulation_results(
    rbind(target_tate, face_tate, target_mean, treated_mean)
  )
  diagnostic <- RoCE:::add_simulation_rmse_comparisons(diagnostic)

  expect_equal(
    diagnostic$rmse_rank[diagnostic$method == "sample_size"],
    1
  )
  expect_equal(
    diagnostic$rmse_relative_to_target[diagnostic$method == "sample_size"],
    0.01
  )
  expect_equal(
    diagnostic$rmse_relative_to_target[
      diagnostic$method == "target_only_ate"
    ],
    1
  )
})

test_that("legacy rows without estimand_scope are inferred and retained", {
  fixture <- make_diagnostic_fixture()
  fixture$estimand_scope <- NULL

  diagnostic <- RoCE:::diagnose_simulation_results(fixture)

  expect_equal(nrow(diagnostic), 1L)
  expect_equal(diagnostic$estimand_scope, "tate")
  expect_equal(diagnostic$n_complete_replications, nrow(fixture))
})

test_that("mixed legacy truncation metadata is grouped without dropping rows", {
  specified <- make_diagnostic_fixture(n = 3L)
  specified$M_tau <- 5
  specified$M_tau_inference <- 5
  legacy <- make_diagnostic_fixture(n = 2L)
  legacy$sim_id <- 4:5
  legacy$M_tau <- NA_real_
  legacy$M_tau_inference <- NA_real_

  diagnostic <- RoCE:::diagnose_simulation_results(
    rbind(specified, legacy), expected_replications = 3L
  )

  expect_equal(nrow(diagnostic), 2L)
  expect_equal(sum(diagnostic$n_rows), 5L)

  ranked <- RoCE:::add_simulation_rmse_comparisons(diagnostic)
  expect_equal(nrow(ranked), 2L)
  expect_equal(sum(ranked$n_rows), 5L)
})

test_that("simulation CSV binding fills columns added by newer tasks", {
  legacy <- make_diagnostic_fixture(n = 2L)
  modern <- make_diagnostic_fixture(n = 2L)
  modern$sim_id <- 3:4
  modern$M_tau <- 5
  modern$M_tau_inference <- 5

  bound <- RoCE:::.bind_simulation_result_frames(list(legacy, modern))

  expect_equal(nrow(bound), 4L)
  expect_true(all(c("M_tau", "M_tau_inference") %in% names(bound)))
  expect_true(all(is.na(bound$M_tau[1:2])))
  expect_equal(bound$M_tau[3:4], c(5, 5))
})

test_that("coverage diagnosis identifies bias-driven undercoverage", {
  fixture <- make_diagnostic_fixture()
  fixture$bias <- seq(0.01, 0.03, length.out = nrow(fixture))
  fixture$estimate <- fixture$truth + fixture$bias
  fixture$se <- rep(stats::sd(fixture$bias), nrow(fixture))
  z <- stats::qnorm(0.975)
  fixture$ci_lower <- fixture$estimate - z * fixture$se
  fixture$ci_upper <- fixture$estimate + z * fixture$se
  fixture$coverage <- fixture$truth >= fixture$ci_lower &
    fixture$truth <= fixture$ci_upper

  diagnostic <- RoCE:::diagnose_simulation_results(fixture)

  expect_lt(diagnostic$coverage, diagnostic$coverage_mc_lower)
  expect_equal(diagnostic$coverage_diagnosis, "primarily_centering_bias")
})

test_that("paired MSE summaries retain common-random-number pairing", {
  focal <- make_diagnostic_fixture(n = 4L, method = "one_round_crossfit_ate")
  benchmark <- make_diagnostic_fixture(n = 4L, method = "target_only_ate")
  focal$bias <- c(1, 2, 3, 4) / 100
  benchmark$bias <- c(2, 1, 4, 3) / 100
  expected_difference <- focal$bias^2 - benchmark$bias^2

  result <- RoCE:::summarize_paired_mse_differences(
    rbind(focal, benchmark[c(4, 2, 1, 3), ]),
    method = "one_round_crossfit_ate",
    benchmark_method = "target_only_ate"
  )

  expect_equal(result$n_paired_replications, 4L)
  expect_equal(
    result$mean_squared_error_difference,
    mean(expected_difference)
  )
  expect_equal(
    result$squared_error_difference_mcse,
    stats::sd(expected_difference) / 2
  )
  expect_equal(
    result$squared_error_difference_z,
    mean(expected_difference) / (stats::sd(expected_difference) / 2)
  )
  expect_equal(
    result$method_lower_squared_error_fraction,
    mean(expected_difference < 0)
  )
})

test_that("TATE paired MSE wrapper covers every standard benchmark", {
  direct <- make_diagnostic_fixture(
    n = 3L, method = "one_round_crossfit_ate"
  )
  benchmarks <- lapply(
    RoCE:::.tate_benchmark_methods(),
    function(method) make_diagnostic_fixture(n = 3L, method = method)
  )

  result <- RoCE:::summarize_tate_paired_mse_benchmarks(
    do.call(rbind, c(list(direct), benchmarks))
  )

  expect_setequal(
    result$benchmark_method,
    RoCE:::.tate_benchmark_methods()
  )
  expect_true(all(result$n_paired_replications == 3L))
})

test_that("paired MSE summaries reject incomplete or ambiguous pairs", {
  focal <- make_diagnostic_fixture(n = 3L, method = "one_round_crossfit_ate")
  benchmark <- make_diagnostic_fixture(n = 3L, method = "target_only_ate")

  expect_error(
    RoCE:::summarize_paired_mse_differences(
      rbind(focal, benchmark[-1L, ]),
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "sim_id sets must match"
  )
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      rbind(focal, focal[1L, ], benchmark),
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "one unique row"
  )
  focal$bias[1L] <- NA_real_
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      rbind(focal, benchmark),
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "must be finite"
  )
  focal <- make_diagnostic_fixture(n = 3L, method = "one_round_crossfit_ate")
  benchmark$sim_id[1L] <- NA_integer_
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      rbind(focal, benchmark),
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "sim_id values must be nonmissing"
  )
  character_bias <- rbind(focal, make_diagnostic_fixture(
    n = 3L, method = "target_only_ate"
  ))
  character_bias$bias <- as.character(character_bias$bias)
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      character_bias,
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "bias must be numeric"
  )
  character_ids <- rbind(focal, make_diagnostic_fixture(
    n = 3L, method = "target_only_ate"
  ))
  character_ids$sim_id <- as.character(character_ids$sim_id)
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      character_ids,
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "sim_id values must be positive finite integers"
  )
  invalid_ids <- rbind(focal, make_diagnostic_fixture(
    n = 3L, method = "target_only_ate"
  ))
  invalid_ids$sim_id[invalid_ids$method == "target_only_ate"] <- 0:2
  expect_error(
    RoCE:::summarize_paired_mse_differences(
      invalid_ids,
      "one_round_crossfit_ate",
      "target_only_ate"
    ),
    "positive finite integers"
  )
})
test_that("source-specific TATE diagnostics preserve source and fold identity", {
  result <- list(
    weights = c(s1 = 0.2, s2 = -0.1),
    target_only = list(estimate = 1.5),
    source_estimates = c(s1 = 1.7, s2 = 1.2),
    fold_weights = matrix(
      c(0.1, -0.2, 0.3, 0), nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2"))
    ),
    fold_wald_statistics = matrix(
      c(0.5, 2, 1.5, 3), nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2"))
    ),
    fold_penalty_coefficients = matrix(
      c(0, 1, 0.5, 2), nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2"))
    ),
    fold_source_included = matrix(
      c(TRUE, FALSE, TRUE, TRUE), nrow = 2L, byrow = TRUE,
      dimnames = list(NULL, c("s1", "s2"))
    ),
    intermediates = list(
      target_estimates = c(1, 2),
      source_estimates_matrix = matrix(
        c(1.1, 0.8, 2.3, 1.7), nrow = 2L, byrow = TRUE,
        dimnames = list(NULL, c("s1", "s2"))
      ),
      phase1b = list(
        avg_target_est = c(10, 20),
        avg_source_est = matrix(
          c(11, 8, 23, 17), nrow = 2L, byrow = TRUE,
          dimnames = list(NULL, c("s1", "s2"))
        ),
        V_ot = c(0, 0),
        V_t = matrix(c(4, 1, 4, 1), nrow = 2L, byrow = TRUE),
        V_s = matrix(0, nrow = 2L, ncol = 2L),
        C_ot = matrix(0, nrow = 2L, ncol = 2L),
        n_t = c(1, 1),
        n_s = matrix(1, nrow = 2L, ncol = 2L)
      )
    )
  )

  diagnostics <- .named_source_diagnostics(result)

  expect_equal(diagnostics[["source_s1_weight"]], 0.2)
  expect_equal(diagnostics[["source_s1_target_estimate"]], 1.5)
  expect_equal(diagnostics[["source_s1_source_estimate"]], 1.7)
  expect_equal(diagnostics[["source_s2_weight"]], -0.1)
  expect_equal(diagnostics[["source_s1_fold_2_weight"]], 0.3)
  expect_equal(diagnostics[["source_s2_fold_1_wald"]], 2)
  expect_equal(
    diagnostics[["source_s2_fold_2_screening_discrepancy"]], -3
  )
  expect_equal(
    diagnostics[["source_s2_fold_2_evaluation_discrepancy"]], -0.3
  )
  expect_equal(
    diagnostics[["source_s1_penalty_activation_fraction"]], 0.5
  )
  expect_equal(
    diagnostics[["source_s2_screening_discrepancy_mean_abs"]], 2.5
  )
  expect_equal(
    diagnostics[["source_s2_evaluation_discrepancy_mean_abs"]], 0.25
  )
  expect_equal(diagnostics[["source_s2_inclusion_fraction"]], 0.5)
  expect_equal(diagnostics[["source_s2_fold_1_included"]], 0)
  expect_equal(
    diagnostics[["source_s2_fold_2_screening_discrepancy_se"]], 1
  )
  expect_equal(diagnostics[["source_s2_wald_identity_error_max"]], 0)

  prefixed <- .named_source_diagnostics(result, prefix = "mu1_")
  expect_true("mu1_source_s1_weight" %in% names(prefixed))

  unnamed <- result
  names(unnamed$weights) <- NULL
  expect_error(
    .named_source_diagnostics(unnamed), "unique source names"
  )
  collided <- result
  names(collided$weights) <- names(collided$source_estimates) <- c("s-1", "s 1")
  colnames(collided$fold_weights) <- colnames(collided$fold_wald_statistics) <-
    colnames(collided$fold_penalty_coefficients) <-
    colnames(collided$fold_source_included) <- c("s-1", "s 1")
  colnames(collided$intermediates$source_estimates_matrix) <-
    colnames(collided$intermediates$phase1b$avg_source_est) <- c("s-1", "s 1")
  expect_error(
    .named_source_diagnostics(collided), "ambiguous diagnostic keys"
  )
  malformed <- result
  malformed$fold_wald_statistics <- malformed$fold_wald_statistics[1, , drop = FALSE]
  expect_error(
    .named_source_diagnostics(malformed), "dimensions disagree"
  )
})
