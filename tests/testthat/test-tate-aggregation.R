library(testthat)

test_that("new aggregation-grid arguments preserve the v32 positional API", {
  v32_formals <- list(
    run_crossfit = c(
      "data_split", "n_folds", "communication_mode", "lambda_selection",
      "lambda_rule", "verbose", "M_tau", "M_tau_inference", "n_cores",
      "nlambda_init", "family", "A_val", "use_lambda_cache",
      "precomputed_folds", "target_only_ps_cache", "target_only_fit_cache",
      "nuisance_lambda_rule"
    ),
    run_tate_crossfit = c(
      "data_split", "n_folds", "communication_mode", "lambda_selection",
      "lambda_rule", "verbose", "M_tau", "M_tau_inference", "n_cores",
      "nlambda_init", "family", "use_lambda_cache", "precomputed_folds",
      "target_only_ps_cache", "target_only_fit_cache",
      "nuisance_lambda_rule", "parallel_arms"
    ),
    calculate_crossfit_aggregation = c(
      "data_split", "target_data", "source_sites", "K", "target_folds",
      "source_folds", "fold_results", "n_folds", "M_tau",
      "M_tau_inference", "lambda_selection", "verbose", "lambda_rule",
      "final_target_estimate", "target_estimates", "source_estimates",
      "source_estimates_matrix", "crossfit_type", "algorithm_label",
      "family_int", "link_int", "A_val"
    ),
    calculate_tate_crossfit_aggregation = c(
      "data_split", "mu1_result", "mu0_result", "lambda_selection",
      "lambda_rule", "verbose"
    ),
    reaggregate_tate_crossfit = c(
      "data_split", "fitted_tate", "M_tau_inference",
      "lambda_selection", "lambda_rule", "verbose"
    )
  )

  for (function_name in names(v32_formals)) {
    current <- names(formals(getExportedValue("RoCE", function_name)))
    expected_prefix <- v32_formals[[function_name]]
    expect_identical(
      current[seq_along(expected_prefix)],
      expected_prefix,
      info = function_name
    )
  }
})

test_that("TATE aggregation grids require inner validation", {
  expect_error(
    run_tate_crossfit(
      data_split = NULL,
      lambda_selection = 0.5,
      aggregation_lambda_grid = c(1, 0.5),
      verbose = FALSE
    ),
    "used only when lambda_selection = 'cv'"
  )
})

test_that("multisite pseudo-value variance removes between-site means", {
  phi <- c(10, 12, -10, -8)

  variance <- RoCE:::.multisite_pseudovalue_variance(phi, c(2L, 2L))

  expect_equal(variance, 4 / 16)
  expect_gt(mean((phi - mean(phi))^2) / length(phi), variance)
})

test_that("pseudo-value assembly requires an exact fold partition", {
  make_info <- function(indices, values) {
    list(
      target_idx = indices,
      varphi_ot = values,
      fold_target_estimate = 0
    )
  }
  valid <- list(
    make_info(c(1, 3), c(10, 30)),
    make_info(c(2, 4), c(20, 40))
  )

  expect_equal(
    RoCE:::.assemble_target_pseudovalues(valid, 4L, "test"),
    c(10, 20, 30, 40)
  )
  expect_error(
    RoCE:::.assemble_target_pseudovalues(
      list(make_info(c(1, 2), c(10, 20)), make_info(c(2, 3), c(20, 30))),
      4L,
      "test"
    ),
    "exactly once"
  )
  expect_error(
    RoCE:::.assemble_target_pseudovalues(
      list(make_info(c(1, 5), c(10, 50))), 4L, "test"
    ),
    "matching valid lengths"
  )
})

test_that("inner-fold summaries pool raw pseudo-values with exact fold sizes", {
  make_component <- function(
      varphi, zeta, xi, target_est, mu_pred, delta) {
    list(
      V_ot = RoCE:::.centered_second_moment(varphi),
      V_t = RoCE:::.centered_second_moment(zeta),
      V_s = RoCE:::.centered_second_moment(xi),
      C_ot = mean((varphi - mean(varphi)) * (zeta - mean(zeta))),
      C_cross = matrix(0, 1, 1),
      avg_target_est = target_est,
      avg_source_est = mu_pred + delta,
      mu_pred_ts = mu_pred,
      delta_ts = delta,
      n_t = length(varphi),
      n_s = length(xi),
      varphi_ot = varphi,
      zeta_components = list(zeta),
      xi_components = list(xi)
    )
  }
  components <- list(
    make_component(c(0, 2), c(1, 3), 1, 10, 12, 8),
    make_component(c(4, 6, 8), c(5, 7, 9), c(3, 5), 20, 25, 15)
  )

  pooled <- RoCE:::.average_aggregation_components(components)
  pooled_varphi <- c(0, 2, 4, 6, 8)
  pooled_zeta <- c(1, 3, 5, 7, 9)

  expect_equal(pooled$n_t, 5)
  expect_equal(pooled$n_s, 3)
  expect_equal(pooled$avg_target_est, 16)
  expected_mu_pred <- (2 * 12 + 3 * 25) / 5
  expected_delta <- (8 + 2 * 15) / 3
  expect_equal(pooled$mu_pred_ts, expected_mu_pred)
  expect_equal(pooled$delta_ts, expected_delta)
  expect_equal(pooled$avg_source_est, expected_mu_pred + expected_delta)
  expect_equal(
    RoCE:::.pool_fold_pairwise_estimates(
      rbind(12, 25), rbind(8, 15), c(2, 3), rbind(1, 2)
    ),
    expected_mu_pred + expected_delta
  )
  expect_equal(
    pooled$V_ot, RoCE:::.centered_second_moment(pooled_varphi)
  )
  expect_equal(
    pooled$V_t, RoCE:::.centered_second_moment(pooled_zeta)
  )
  expect_equal(
    pooled$C_ot,
    mean(
      (pooled_varphi - mean(pooled_varphi)) *
        (pooled_zeta - mean(pooled_zeta))
    )
  )
})

test_that("single-source fold extraction preserves fold-by-source shape", {
  fold_results <- lapply(1:3, function(k) {
    list(source_results = list(s1 = list(mu_ts = k)))
  })
  extracted <- RoCE:::extract_source_result_matrix(
    fold_results, "s1", n_folds = 3L, field = "mu_ts"
  )

  expect_equal(dim(extracted), c(3L, 1L))
  expect_equal(as.numeric(extracted), 1:3)
})

test_that("TATE fold information is formed before variance calculation", {
  make_info <- function(sign) {
    list(
      fold_target_estimate = sign * 2,
      fold_source_estimates = sign * c(1.5, 1.8),
      mu_pred_ts = sign * c(1.1, 1.3),
      delta_ts = sign * c(0.4, 0.5),
      varphi_ot = sign * c(-1, 0, 1),
      zeta_components = list(
        sign * c(-0.5, 0, 0.5),
        sign * c(-1, 0.5, 0.5)
      ),
      xi_components = list(
        sign * c(-0.25, 0, 0.25),
        sign * c(-0.5, 0.25, 0.25)
      ),
      source_variance_clip_diagnostics = vector("list", 2L),
      target_idx = 1:3,
      source_idx = list(1:3, 1:3)
    )
  }

  treated <- make_info(1)
  control <- make_info(-1)
  tate <- RoCE:::.combine_tate_fold_info(treated, control, K = 2L)

  expect_equal(tate$fold_target_estimate, 4)
  expect_equal(tate$fold_source_estimates, c(3, 3.6))
  expect_equal(tate$varphi_ot, c(-2, 0, 2))
  expect_equal(tate$V_ot, mean(c(-2, 0, 2)^2))
  expect_equal(tate$zeta_components[[1]], c(-1, 0, 1))
  expect_equal(tate$xi_components[[1]], c(-0.5, 0, 0.5))
  expect_true(all(is.finite(tate$C_cross)))
})

test_that("TATE variance components retain cross-arm covariance", {
  make_info <- function(varphi, zeta, xi) {
    list(
      fold_target_estimate = mean(varphi),
      fold_source_estimates = mean(zeta) + mean(xi),
      mu_pred_ts = mean(zeta),
      delta_ts = mean(xi),
      varphi_ot = varphi,
      zeta_components = list(zeta),
      xi_components = list(xi),
      source_variance_clip_diagnostics = list(NULL),
      target_idx = seq_along(varphi),
      source_idx = list(seq_along(xi))
    )
  }
  treated <- make_info(
    varphi = c(-2, 0, 2),
    zeta = c(-1, 0, 1),
    xi = c(-3, 0, 3)
  )
  control <- make_info(
    varphi = c(-1, 0, 1),
    zeta = c(1, 0, -1),
    xi = c(-1, 0, 1)
  )

  tate <- RoCE:::.combine_tate_fold_info(treated, control, K = 1L)
  contrast_variance <- function(treated_values, control_values) {
    treated_centered <- treated_values - mean(treated_values)
    control_centered <- control_values - mean(control_values)
    RoCE:::.centered_second_moment(treated_values) +
      RoCE:::.centered_second_moment(control_values) -
      2 * mean(treated_centered * control_centered)
  }

  expect_equal(
    tate$V_ot,
    contrast_variance(treated$varphi_ot, control$varphi_ot)
  )
  expect_equal(
    tate$V_t[[1L]],
    contrast_variance(
      treated$zeta_components[[1L]], control$zeta_components[[1L]]
    )
  )
  expect_equal(
    tate$V_s[[1L]],
    contrast_variance(
      treated$xi_components[[1L]], control$xi_components[[1L]]
    )
  )
  expect_equal(tate$varphi_ot, treated$varphi_ot - control$varphi_ot)
})

test_that("expanded TATE aggregation variance matches sitewise pseudo-values", {
  varphi <- c(-1.5, -0.5, 0.5, 1.5)
  zeta <- list(
    c(-0.6, 0.2, -0.2, 0.6),
    c(0.4, -0.8, 0.8, -0.4)
  )
  xi <- list(c(-1, 0, 1), c(-0.5, 0.5))
  eta <- c(0.25, -0.1)
  n_t <- length(varphi)
  n_s <- vapply(xi, length, integer(1L))

  V_ot <- RoCE:::.centered_second_moment(varphi)
  V_t <- vapply(zeta, RoCE:::.centered_second_moment, numeric(1L))
  V_s <- vapply(xi, RoCE:::.centered_second_moment, numeric(1L))
  C_ot <- vapply(zeta, function(value) mean(
    (varphi - mean(varphi)) * (value - mean(value))
  ), numeric(1L))
  C_cross <- matrix(0, 2L, 2L)
  C_cross[1L, 2L] <- C_cross[2L, 1L] <- mean(
    (zeta[[1L]] - mean(zeta[[1L]])) *
      (zeta[[2L]] - mean(zeta[[2L]]))
  )

  expanded <- calculate_aggregated_variance_cpp(
    eta, V_t, V_s, n_s, V_ot, n_t, C_ot, C_cross, 0, 0
  )
  target_pseudo <- (1 - sum(eta)) * varphi +
    eta[[1L]] * zeta[[1L]] + eta[[2L]] * zeta[[2L]]
  direct <- sum((target_pseudo - mean(target_pseudo))^2) / n_t^2
  for (j in seq_along(xi)) {
    source_pseudo <- eta[[j]] * xi[[j]]
    direct <- direct +
      sum((source_pseudo - mean(source_pseudo))^2) / n_s[[j]]^2
  }

  expect_equal(expanded, direct, tolerance = 1e-14)
  for (j in seq_along(xi)) {
    discrepancy_formula <-
      (V_ot + V_t[[j]] - 2 * C_ot[[j]]) / n_t + V_s[[j]] / n_s[[j]]
    target_discrepancy <- zeta[[j]] - varphi
    discrepancy_direct <-
      sum((target_discrepancy - mean(target_discrepancy))^2) / n_t^2 +
      sum((xi[[j]] - mean(xi[[j]]))^2) / n_s[[j]]^2
    expect_equal(
      discrepancy_formula, discrepancy_direct, tolerance = 1e-14
    )
  }
  expect_equal(
    calculate_aggregated_variance_cpp(
      c(0, 0), V_t, V_s, n_s, V_ot, n_t, C_ot, C_cross, 0, 0
    ),
    V_ot / n_t,
    tolerance = 1e-14
  )
})

test_that("TATE cross-fitting returns one common source-weight vector", {
  skip_on_cran()

  set.seed(812)
  data <- generate_simulation_data(
    n_total = 360, K = 1, p = 3, config = "C1",
    dgp_type = "face",
    n_target = 120, n_source_sizes = 240,
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  target_only_fit_cache <- new.env(hash = TRUE, parent = emptyenv())

  result <- run_tate_crossfit(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    family = "gaussian",
    nlambda_init = 5,
    n_cores = 1,
    target_only_fit_cache = target_only_fit_cache,
    verbose = FALSE
  )

  expect_identical(result$estimand, "TATE")
  expect_length(result$weights, 1L)
  expect_equal(result$fold_weights[, 1L], as.numeric(result$fold_weights[, 1L]))
  expect_true(all(is.finite(c(result$estimate, result$se, result$weights))))
  expect_gt(result$se, 0)
  expect_true(is.finite(result$timing$total_seconds))
  expect_gt(result$timing$total_seconds, 0)
  expect_equal(nrow(result$timing$mu1$folds), result$n_folds)
  expect_equal(nrow(result$timing$mu0$folds), result$n_folds)
  expect_true(all(vapply(
    result$timing[c("mu1", "mu0")],
    function(arm_timing) all(is.finite(as.matrix(arm_timing$folds))),
    logical(1L)
  )))
  expect_equal(nrow(result$nuisance_fit_diagnostics), 2L * result$n_folds)
  expect_setequal(unique(result$nuisance_fit_diagnostics$A_val), c(0L, 1L))
  expect_true(all(
    result$nuisance_fit_diagnostics$initial_dr_nonconverged >= 0
  ))
  expect_true(all(c(
    "calibrated_outcome_line_search_failures",
    "max_calibrated_outcome_update_ratio",
    "max_calibrated_outcome_abs_coefficient"
  ) %in% names(result$nuisance_fit_diagnostics)))
  expect_true(!is.null(result$clip_diagnostics$by_arm$treated))
  expect_true(!is.null(result$clip_diagnostics$by_arm$control))
  expect_equal(
    result$clip_diagnostics$total$logit_truncated,
    sum(vapply(
      result$clip_diagnostics$by_arm,
      function(arm) arm$total$logit_truncated,
      numeric(1L)
    ))
  )
  expect_equal(
    result$clip_diagnostics$total$logit_truncation_fraction,
    result$clip_diagnostics$total$logit_truncated /
      result$clip_diagnostics$total$n_obs
  )
  expect_length(ls(target_only_fit_cache, all.names = TRUE), 18L)
  expect_equal(
    result$fold_lambdas,
    rep(RoCE:::AGG_WALD_LAMBDA, result$n_folds)
  )
  expect_equal(
    dim(result$fold_wald_statistics),
    c(result$n_folds, length(result$weights))
  )
  expect_equal(
    result$fold_penalty_coefficients,
    pmax(
      RoCE:::AGG_WALD_LAMBDA * result$fold_wald_statistics - 1,
      0
    )
  )

  hard_screened <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = result$arm_results$mu1,
    mu0_result = result$arm_results$mu0,
    lambda_selection = RoCE:::AGG_WALD_LAMBDA,
    screening_rule = "hard_threshold",
    verbose = FALSE
  )
  expect_identical(
    hard_screened$aggregation_screening_rule, "hard_threshold"
  )
  expect_identical(
    hard_screened$fold_source_included,
    hard_screened$fold_wald_statistics <=
      1 / hard_screened$fold_lambdas
  )
  expect_true(all(
    hard_screened$fold_weights[!hard_screened$fold_source_included] == 0
  ))
  expect_error(
    calculate_tate_crossfit_aggregation(
      data_split = data_split,
      mu1_result = result$arm_results$mu1,
      mu0_result = result$arm_results$mu0,
      lambda_selection = 0,
      screening_rule = "hard_threshold",
      verbose = FALSE
    ),
    "strictly positive"
  )

  quadratic <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = result$arm_results$mu1,
    mu0_result = result$arm_results$mu0,
    lambda_selection = RoCE:::AGG_WALD_LAMBDA,
    screening_rule = "quadratic_bias",
    verbose = FALSE
  )
  expect_identical(quadratic$aggregation_screening_rule, "quadratic_bias")
  expect_true(all(quadratic$fold_source_included))
  expect_true(all(quadratic$fold_weight_psd_ridge == 0))
  expect_identical(quadratic$weight_layer$kink_cells, 0L)
  # The Wald diagnostics are descriptive and identical across rules.
  expect_equal(quadratic$fold_wald_statistics, result$fold_wald_statistics)
  expect_equal(quadratic$fold_penalty_coefficients, result$fold_penalty_coefficients)
  expect_true(all(is.finite(c(quadratic$estimate, quadratic$se, quadratic$se_fixed_weights))))
  # The pseudo-value variance is the fixed-weight part; `variance` adds the
  # weight layer (HISTORY #0005).
  expect_equal(
    hard_screened$variance_fixed_weights,
    RoCE:::.multisite_pseudovalue_variance(
      hard_screened$all_phi_tau,
      vapply(data_split, function(site) site$n, integer(1L))
    )
  )

  # A caller-supplied aggregation grid must be evaluated only by the
  # inner-validation rule. Reuse the already fitted arm nuisances so this test
  # isolates the low-dimensional TATE aggregation layer.
  candidate_grid <- c(1, 0.5, 1 / 3)
  grid_selected <- calculate_tate_crossfit_aggregation(
    data_split = data_split,
    mu1_result = result$arm_results$mu1,
    mu0_result = result$arm_results$mu0,
    lambda_selection = "cv",
    lambda_rule = "min",
    aggregation_lambda_grid = candidate_grid,
    verbose = FALSE
  )
  expect_true(all(grid_selected$fold_lambdas %in% candidate_grid))
  expect_identical(grid_selected$aggregation_lambda_selection, "cv")
  expect_equal(grid_selected$aggregation_lambda_grid, candidate_grid)
  expect_true(all(vapply(
    grid_selected$fold_lambda_info,
    function(info) length(info$cv_scores) == length(candidate_grid),
    logical(1L)
  )))
  expect_error(
    calculate_tate_crossfit_aggregation(
      data_split = data_split,
      mu1_result = result$arm_results$mu1,
      mu0_result = result$arm_results$mu0,
      lambda_selection = 0.5,
      aggregation_lambda_grid = candidate_grid,
      verbose = FALSE
    ),
    "used only when lambda_selection = 'cv'"
  )
  expect_error(
    calculate_tate_crossfit_aggregation(
      data_split = data_split,
      mu1_result = result$arm_results$mu1,
      mu0_result = result$arm_results$mu0,
      lambda_selection = "cv",
      aggregation_lambda_grid = c(0.5, 0.5),
      verbose = FALSE
    ),
    "must remain unique"
  )
  expected_target_training_n <- vapply(
    result$arm_results$mu1$intermediates$fold_info,
    function(info) data_split$t$n - length(info$target_idx),
    numeric(1L)
  )
  expect_equal(result$intermediates$phase1b$n_t, expected_target_training_n)
  expect_equal(
    vapply(
      result$fold_training_sample_sizes,
      function(sizes) sizes$n_t,
      numeric(1L)
    ),
    expected_target_training_n
  )
  expect_equal(
    do.call(rbind, lapply(
      result$fold_training_sample_sizes,
      function(sizes) sizes$n_s
    )),
    result$intermediates$phase1b$n_s
  )
  for (k1 in seq_len(result$n_folds)) {
    phase1b <- result$intermediates$phase1b
    discrepancy_variance <-
      phase1b$V_ot[k1] / phase1b$n_t[k1] +
      phase1b$V_t[k1, ] / phase1b$n_t[k1] +
      phase1b$V_s[k1, ] / phase1b$n_s[k1, ] -
      2 * phase1b$C_ot[k1, ] / phase1b$n_t[k1]
    expected_wald <- abs(
      phase1b$avg_target_est[k1] - phase1b$avg_source_est[k1, ]
    ) / sqrt(pmax(discrepancy_variance, RoCE:::VARIANCE_MIN))
    expect_equal(result$fold_wald_statistics[k1, ], expected_wald)
  }
  expect_equal(
    result$variance_fixed_weights,
    RoCE:::.multisite_pseudovalue_variance(
      result$all_phi_tau,
      vapply(data_split, function(site) site$n, integer(1L))
    )
  )
  expect_equal(
    result$variance,
    result$variance_fixed_weights + result$weight_layer$indirect_variance +
      result$weight_layer$cross_term,
    tolerance = 1e-12
  )
  expect_equal(result$estimate, mean(result$all_phi_tau))
  expect_equal(
    result$intermediates$target_only_phi,
    result$arm_results$mu1$intermediates$target_only_phi -
      result$arm_results$mu0$intermediates$target_only_phi
  )
  expect_equal(
    result$target_only$estimate,
    mean(result$intermediates$target_only_phi)
  )
  for (arm_name in c("mu1", "mu0")) {
    arm <- result$arm_results[[arm_name]]
    expect_equal(arm$target_only$estimate, mean(arm$intermediates$target_only_phi))
    expect_equal(
      arm$target_only$variance,
      RoCE:::.multisite_pseudovalue_variance(
        arm$intermediates$target_only_phi,
        data_split$t$n
      )
    )
    expect_equal(arm$target_only$se^2, arm$target_only$variance)
  }

  # Reaggregation at the fitted inference radius must be an identity. This is
  # the safety gate that lets cutoff and inference-radius sensitivities reuse
  # the expensive high-dimensional nuisance fits.
  reaggregated <- reaggregate_tate_crossfit(
    data_split = data_split,
    fitted_tate = result,
    M_tau_inference = result$M_tau_inference,
    verbose = FALSE
  )
  expect_true(reaggregated$reused_nuisance_fits)
  expect_equal(reaggregated$estimate, result$estimate, tolerance = 1e-12)
  expect_equal(reaggregated$se, result$se, tolerance = 1e-12)
  expect_equal(reaggregated$weights, result$weights, tolerance = 1e-12)
  expect_equal(
    reaggregated$fold_wald_statistics,
    result$fold_wald_statistics,
    tolerance = 1e-12
  )
  expect_equal(reaggregated$all_phi_tau, result$all_phi_tau,
               tolerance = 1e-12)
  expect_equal(reaggregated$timing$nuisance_refit_seconds, 0)

  # Changing neither radius nor cutoff must preserve every supported rule,
  # including fitted objects from before the screening metadata was added.
  for (rule in c("soft_penalty", "hard_threshold", "quadratic_bias")) {
    original <- result
    aggregated <- calculate_tate_crossfit_aggregation(
      data_split, result$arm_results$mu1, result$arm_results$mu0,
      lambda_selection = result$aggregation_lambda_selection,
      screening_rule = rule, verbose = FALSE
    )
    original[names(aggregated)] <- aggregated
    refreshed <- reaggregate_tate_crossfit(
      data_split, original, original$M_tau_inference, verbose = FALSE
    )
    expect_identical(refreshed$aggregation_screening_rule, rule)
    for (field in c("estimate", "se", "weights", "fold_weights")) {
      expect_equal(refreshed[[field]], original[[field]], tolerance = 1e-12)
    }
    grid <- reaggregate_tate_sensitivity_grid(
      data_split, original,
      data.frame(cutoff = c(1, 2, 2), M_tau_inference = c(5, 5, Inf))
    )
    for (i in seq_along(grid$results)) {
      candidate <- grid$results[[i]]
      expected <- calculate_tate_crossfit_aggregation(
        data_split, candidate$arm_results$mu1, candidate$arm_results$mu0,
        lambda_selection = 1 / grid$grid$cutoff[i], screening_rule = rule,
        verbose = FALSE
      )
      expect_identical(candidate$aggregation_screening_rule, rule)
      expect_equal(candidate$estimate, expected$estimate, tolerance = 1e-12)
      expect_equal(candidate$se, expected$se, tolerance = 1e-12)
      expect_equal(candidate$weights, expected$weights, tolerance = 1e-12)
    }
  }
  legacy <- result
  legacy$aggregation_screening_rule <- NULL
  legacy_refreshed <- reaggregate_tate_crossfit(
    data_split, legacy, legacy$M_tau_inference, verbose = FALSE
  )
  expect_identical(legacy_refreshed$aggregation_screening_rule, "soft_penalty")
  expect_equal(legacy_refreshed$estimate, result$estimate, tolerance = 1e-12)
  for (arm_name in c("mu1", "mu0")) {
    for (fold_index in seq_len(result$n_folds)) {
      original_sources <-
        result$arm_results[[arm_name]]$fold_results[[fold_index]]$source_results
      reused_arm <- reaggregated$arm_results[[arm_name]]
      reused_sources <- reused_arm$fold_results[[fold_index]]$source_results
      for (source_name in names(original_sources)) {
        expect_identical(
          reused_sources[[source_name]]$gamma_s,
          original_sources[[source_name]]$gamma_s
        )
        expect_identical(
          reused_sources[[source_name]]$alpha_ts,
          original_sources[[source_name]]$alpha_ts
        )
      }
    }
  }

  untruncated_inference <- reaggregate_tate_crossfit(
    data_split = data_split,
    fitted_tate = result,
    M_tau_inference = Inf,
    lambda_selection = 1 / 2.5,
    verbose = FALSE
  )
  expect_true(all(is.finite(c(
    untruncated_inference$estimate,
    untruncated_inference$se,
    untruncated_inference$weights
  ))))
  expect_equal(untruncated_inference$M_tau, result$M_tau)
  expect_equal(untruncated_inference$M_tau_inference, Inf)
  expect_equal(untruncated_inference$aggregation_lambda_selection, 1 / 2.5)
  expect_null(untruncated_inference$aggregation_lambda_grid)

  cv_grid <- c(1, 0.5)
  cv_template <- result
  cv_template$aggregation_lambda_selection <- "cv"
  cv_template$aggregation_lambda_grid <- cv_grid
  cv_reaggregated <- reaggregate_tate_crossfit(
    data_split = data_split,
    fitted_tate = cv_template,
    M_tau_inference = result$M_tau_inference,
    verbose = FALSE
  )
  expect_equal(cv_reaggregated$aggregation_lambda_grid, cv_grid)
  expect_true(all(cv_reaggregated$fold_lambdas %in% cv_grid))
  expect_error(
    reaggregate_tate_crossfit(
      data_split = data_split,
      fitted_tate = result,
      M_tau_inference = result$M_tau_inference,
      lambda_selection = 0.5,
      aggregation_lambda_grid = cv_grid,
      verbose = FALSE
    ),
    "used only when lambda_selection = 'cv'"
  )

  sensitivity <- reaggregate_tate_sensitivity_grid(
    data_split = data_split,
    fitted_tate = result,
    sensitivity_grid = data.frame(
      cutoff = c(AGG_WALD_CUTOFF, 1.5, 2, 2.5),
      M_tau_inference = c(5, 5, Inf, Inf)
    ),
    verbose = FALSE
  )
  expect_equal(sensitivity$n_nuisance_refits, 0L)
  expect_equal(sensitivity$n_inference_refreshes, 1L)
  expect_identical(sensitivity$results[[1L]], result)
  expect_equal(
    sensitivity$results[[3L]]$estimate,
    reaggregate_tate_crossfit(
      data_split, result, M_tau_inference = Inf,
      lambda_selection = 0.5, verbose = FALSE
    )$estimate,
    tolerance = 1e-12
  )
  expect_equal(
    vapply(
      sensitivity$results,
      function(fit) fit$timing$nuisance_refit_seconds %||% 0,
      numeric(1L)
    ),
    rep(0, 4L)
  )
  expect_error(
    reaggregate_tate_sensitivity_grid(
      data_split, result,
      data.frame(cutoff = c(2, 2), M_tau_inference = c(5, 5))
    ),
    "rows must be unique"
  )
})

test_that("arm-specific aggregation preserves source identity on diagnostics", {
  skip_on_cran()

  set.seed(913)
  data <- generate_simulation_data(
    n_total = 360, K = 1, p = 2, config = "C1",
    dgp_type = "face", n_target = 120, n_source_sizes = 240,
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  fit <- run_crossfit(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    lambda_selection = AGG_WALD_LAMBDA,
    verbose = FALSE,
    n_cores = 1,
    nlambda_init = 5,
    A_val = 1L
  )

  expect_identical(names(fit$weights), "s1")
  expect_identical(colnames(fit$fold_weights), "s1")
  expect_identical(colnames(fit$fold_wald_statistics), "s1")
  expect_identical(colnames(fit$fold_penalty_coefficients), "s1")
  diagnostics <- .named_source_diagnostics(fit, prefix = "mu1_")
  expect_true(length(diagnostics) > 0L)
  expect_true(all(is.finite(diagnostics)))
  expect_true("mu1_source_s1_weight" %in% names(diagnostics))
})

test_that("concurrent treatment arms are reproducible and resource-safe", {
  skip_on_cran()
  skip_on_os("windows")

  set.seed(813)
  data <- generate_simulation_data(
    n_total = 360, K = 1, p = 3, config = "C1",
    dgp_type = "face",
    n_target = 120, n_source_sizes = 240,
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)

  sequential <- run_tate_crossfit(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    family = "gaussian",
    nlambda_init = 3,
    n_cores = 1,
    parallel_arms = FALSE,
    verbose = FALSE
  )
  concurrent <- run_tate_crossfit(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    family = "gaussian",
    nlambda_init = 3,
    n_cores = 1,
    parallel_arms = TRUE,
    verbose = FALSE
  )

  expect_true(concurrent$parallel_arms)
  expect_true(all(is.finite(c(
    concurrent$estimate, concurrent$se, concurrent$weights
  ))))
  expect_equal(concurrent$estimate, sequential$estimate, tolerance = 1e-12)
  expect_equal(concurrent$se, sequential$se, tolerance = 1e-12)
  expect_equal(concurrent$weights, sequential$weights, tolerance = 1e-12)
  expect_equal(
    concurrent$fold_wald_statistics,
    sequential$fold_wald_statistics,
    tolerance = 1e-12
  )
  expect_error(
    run_tate_crossfit(
      data_split = data_split,
      n_folds = 3,
      n_cores = -1,
      parallel_arms = TRUE,
      verbose = FALSE
    ),
    "explicit positive n_cores"
  )
})

test_that("source-site parallelism preserves the cross-fit estimator", {
  skip_on_cran()
  skip_on_os("windows")

  set.seed(814)
  data <- generate_simulation_data(
    n_total = 450, K = 2, p = 2, config = "C1",
    dgp_type = "face", outcome_type = "continuous",
    n_target = 150, n_source_sizes = c(150, 150),
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  common_args <- list(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    family = "gaussian",
    A_val = 1L,
    nlambda_init = 3,
    verbose = FALSE
  )

  set.seed(815)
  sequential <- do.call(run_crossfit, c(common_args, list(n_cores = 1)))
  set.seed(815)
  concurrent <- do.call(run_crossfit, c(common_args, list(n_cores = 2)))

  expect_equal(concurrent$estimate, sequential$estimate, tolerance = 1e-12)
  expect_equal(concurrent$se, sequential$se, tolerance = 1e-12)
  expect_equal(concurrent$weights, sequential$weights, tolerance = 1e-12)
  expect_equal(
    concurrent$fold_wald_statistics,
    sequential$fold_wald_statistics,
    tolerance = 1e-12
  )
})

test_that("comparison TATE variance retains within-site arm covariance", {
  mu1 <- list(
    estimate = 0.7,
    n = 2L,
    influence_blocks = list(
      t = list(influence = c(-1, 1), weight = 1)
    )
  )
  mu0 <- list(
    estimate = 0.4,
    n = 2L,
    influence_blocks = list(
      t = list(influence = c(-0.5, 0.5), weight = 1)
    )
  )

  result <- RoCE:::.combine_comparison_arms(
    mu1, mu0, method = "target_only", variance_method = "analytic"
  )

  expect_equal(result$estimate, 0.3)
  expect_equal(result$variance, 0.125)
  expect_equal(result$se, sqrt(0.125))
  expect_equal(
    result$influence_blocks$t$influence,
    c(-0.5, 0.5)
  )
})

test_that("comparison TATE results retain density-ratio QC diagnostics", {
  dr <- list(
    dr_weight_n = 20L,
    dr_weight_n_clipped = 2L,
    dr_weight_fraction_clipped = 0.1,
    dr_weight_max_site_fraction_clipped = 0.2,
    dr_weight_min_before_clipping = 0.03,
    dr_weight_max_before_clipping = 17
  )
  make_arm <- function(estimate) {
    list(
      estimate = estimate,
      n = 2L,
      influence_blocks = list(
        t = list(influence = c(-1, 1), weight = 1)
      ),
      components = list(dr_weight_diagnostics = dr)
    )
  }

  result <- RoCE:::.combine_comparison_arms(
    make_arm(0.7), make_arm(0.4), method = "federated_dr",
    variance_method = "analytic"
  )
  diagnostic <- RoCE:::.comparison_variance_diagnostics(result)

  expect_identical(result$components$dr_weight_diagnostics, dr)
  expect_equal(diagnostic$dr_weight_n, 20L)
  expect_equal(diagnostic$dr_weight_fraction_clipped, 0.1)
  expect_equal(diagnostic$dr_weight_max_before_clipping, 17)

  mu0 <- make_arm(0.4)
  mu0$components$dr_weight_diagnostics$dr_weight_n_clipped <- 3L
  expect_error(
    RoCE:::.combine_comparison_arms(
      make_arm(0.7), mu0, method = "federated_dr",
      variance_method = "analytic"
    ),
    "inconsistent density-ratio diagnostics"
  )
})

test_that("arm timing summaries are additive and retain stage names", {
  make_timing_result <- function(total, multiplier) {
    list(timing = list(
      total_seconds = total,
      folds = data.frame(
        fold = 1:2,
        initial_nuisance_seconds = multiplier * c(1, 2),
        source_processing_wall_seconds = multiplier * c(3, 4),
        target_only_seconds = multiplier * c(5, 6)
      )
    ))
  }

  timing <- RoCE:::.summarize_tate_crossfit_timing(
    make_timing_result(10, 1), make_timing_result(20, 2)
  )

  expect_equal(timing[["face_sum_arm_seconds"]], 30)
  expect_equal(timing[["face_initial_nuisance_seconds"]], 9)
  expect_equal(timing[["face_source_processing_wall_seconds"]], 21)
  expect_equal(timing[["face_target_only_seconds"]], 33)
})

test_that("TATE nuisance diagnostics combine counts and iteration maxima", {
  make_diagnostic_result <- function(multiplier) {
    list(nuisance_fit_diagnostics = data.frame(
      fold = 1:2,
      initial_outcome_degenerate = multiplier * c(1, 0),
      target_only_outcome_degenerate = multiplier * c(0, 1),
      initial_dr_nonconverged = multiplier * c(0, 1),
      calibrated_dr_nonconverged = multiplier * c(1, 0),
      calibrated_outcome_nonconverged = multiplier * c(0, 2),
      initial_dr_cv_invalid_fold_fits = multiplier * c(1, 2),
      initial_dr_cv_invalid_lambdas = multiplier * c(0, 1),
      initial_dr_cv_path_tail_skipped_fold_fits = multiplier * c(0, 2),
      calibrated_dr_cv_invalid_fold_fits = multiplier * c(3, 0),
      calibrated_dr_cv_invalid_lambdas = multiplier * c(1, 0),
      calibrated_dr_cv_path_tail_skipped_fold_fits = multiplier * c(2, 0),
      calibrated_outcome_cv_invalid_fold_fits = multiplier * c(0, 4),
      calibrated_outcome_cv_invalid_lambdas = multiplier * c(0, 2),
      calibrated_outcome_cv_path_tail_skipped_fold_fits =
        multiplier * c(0, 3),
      initial_dr_line_search_failures = multiplier * c(0, 3),
      calibrated_dr_line_search_failures = multiplier * c(2, 0),
      calibrated_outcome_line_search_failures = multiplier * c(1, 4),
      initial_dr_support_floor_applied = multiplier * c(1, 0),
      calibrated_dr_support_floor_applied = multiplier * c(0, 1),
      max_initial_dr_iterations = multiplier * c(5, 8),
      max_calibrated_dr_iterations = multiplier * c(7, 9),
      max_calibrated_outcome_iterations = multiplier * c(4, 6),
      max_initial_dr_update_ratio = multiplier * c(0.4, 0.8),
      max_initial_dr_abs_coefficient = multiplier * c(1.0, 1.3),
      max_calibrated_dr_update_ratio = multiplier * c(0.5, 0.9),
      max_calibrated_dr_abs_coefficient = multiplier * c(1.1, 1.4),
      max_calibrated_outcome_update_ratio = multiplier * c(0.6, 0.7),
      max_calibrated_outcome_abs_coefficient = multiplier * c(1.2, 1.5),
      max_initial_dr_support_floor = multiplier * c(0.1, 0.2),
      max_calibrated_dr_support_floor = multiplier * c(0.3, 0.4)
    ))
  }

  diagnostics <- RoCE:::.summarize_tate_nuisance_diagnostics(
    make_diagnostic_result(1), make_diagnostic_result(2)
  )

  expect_equal(diagnostics[["face_initial_outcome_degenerate"]], 3)
  expect_equal(diagnostics[["face_target_only_outcome_degenerate"]], 3)
  expect_equal(diagnostics[["face_mu1_initial_outcome_degenerate"]], 1)
  expect_equal(diagnostics[["face_mu0_initial_outcome_degenerate"]], 2)
  expect_equal(
    diagnostics[["face_mu1_target_only_outcome_degenerate"]], 1
  )
  expect_equal(
    diagnostics[["face_mu0_target_only_outcome_degenerate"]], 2
  )
  expect_equal(diagnostics[["face_initial_dr_nonconverged"]], 3)
  expect_equal(diagnostics[["face_calibrated_dr_nonconverged"]], 3)
  expect_equal(diagnostics[["face_calibrated_outcome_nonconverged"]], 6)
  expect_equal(diagnostics[["face_initial_dr_cv_invalid_fold_fits"]], 9)
  expect_equal(diagnostics[["face_initial_dr_cv_invalid_lambdas"]], 3)
  expect_equal(
    diagnostics[["face_initial_dr_cv_path_tail_skipped_fold_fits"]], 6
  )
  expect_equal(diagnostics[["face_calibrated_dr_cv_invalid_fold_fits"]], 9)
  expect_equal(diagnostics[["face_calibrated_dr_cv_invalid_lambdas"]], 3)
  expect_equal(
    diagnostics[["face_calibrated_dr_cv_path_tail_skipped_fold_fits"]], 6
  )
  expect_equal(
    diagnostics[["face_calibrated_outcome_cv_invalid_fold_fits"]], 12
  )
  expect_equal(
    diagnostics[["face_calibrated_outcome_cv_invalid_lambdas"]], 6
  )
  expect_equal(
    diagnostics[["face_calibrated_outcome_cv_path_tail_skipped_fold_fits"]],
    9
  )
  expect_equal(diagnostics[["face_initial_dr_line_search_failures"]], 9)
  expect_equal(diagnostics[["face_calibrated_dr_line_search_failures"]], 6)
  expect_equal(diagnostics[["face_calibrated_outcome_line_search_failures"]], 15)
  expect_equal(diagnostics[["face_initial_dr_support_floor_applied"]], 3)
  expect_equal(diagnostics[["face_calibrated_dr_support_floor_applied"]], 3)
  expect_equal(diagnostics[["face_mu1_initial_dr_nonconverged"]], 1)
  expect_equal(diagnostics[["face_mu0_initial_dr_nonconverged"]], 2)
  expect_equal(diagnostics[["face_mu1_calibrated_dr_nonconverged"]], 1)
  expect_equal(diagnostics[["face_mu0_calibrated_dr_nonconverged"]], 2)
  expect_equal(diagnostics[["face_mu1_calibrated_outcome_nonconverged"]], 2)
  expect_equal(diagnostics[["face_mu0_calibrated_outcome_nonconverged"]], 4)
  expect_equal(diagnostics[["face_mu1_initial_dr_line_search_failures"]], 3)
  expect_equal(diagnostics[["face_mu0_initial_dr_line_search_failures"]], 6)
  expect_equal(diagnostics[["face_mu1_calibrated_dr_line_search_failures"]], 2)
  expect_equal(diagnostics[["face_mu0_calibrated_dr_line_search_failures"]], 4)
  expect_equal(
    diagnostics[["face_mu1_calibrated_outcome_line_search_failures"]], 5
  )
  expect_equal(
    diagnostics[["face_mu0_calibrated_outcome_line_search_failures"]], 10
  )
  expect_equal(diagnostics[["face_max_initial_dr_iterations"]], 16)
  expect_equal(diagnostics[["face_max_calibrated_dr_iterations"]], 18)
  expect_equal(diagnostics[["face_max_calibrated_outcome_iterations"]], 12)
  expect_equal(diagnostics[["face_max_initial_dr_update_ratio"]], 1.6)
  expect_equal(diagnostics[["face_max_initial_dr_abs_coefficient"]], 2.6)
  expect_equal(diagnostics[["face_max_calibrated_dr_update_ratio"]], 1.8)
  expect_equal(diagnostics[["face_max_calibrated_dr_abs_coefficient"]], 2.8)
  expect_equal(diagnostics[["face_max_calibrated_outcome_update_ratio"]], 1.4)
  expect_equal(diagnostics[["face_max_calibrated_outcome_abs_coefficient"]], 3)
  expect_equal(diagnostics[["face_max_initial_dr_support_floor"]], 0.4)
  expect_equal(diagnostics[["face_max_calibrated_dr_support_floor"]], 0.8)
  expect_true(is.na(RoCE:::.max_iterations_or_na(c(NA, -1))))
  expect_true(is.na(RoCE:::.max_numeric_or_na(c(NA, Inf))))
})

test_that("TATE degenerate counts are explicit and additive", {
  make_result <- function(initial, target_only) {
    list(nuisance_fit_diagnostics = data.frame(
      initial_outcome_degenerate = initial,
      target_only_outcome_degenerate = target_only
    ))
  }
  expect_identical(
    RoCE:::.count_direct_tate_degenerate_fits(list(
      one_round = make_result(c(1, 0), c(0, 2)),
      two_round = make_result(c(0, 3), c(1, 0))
    )),
    7L
  )
  expect_true(is.na(
    RoCE:::.count_direct_tate_degenerate_fits(list())
  ))
})

test_that("comparison TATE rejects misaligned arm influence blocks", {
  mu1 <- list(
    estimate = 0.7,
    n = 2L,
    influence_blocks = list(
      t = list(influence = c(-1, 1), weight = 1)
    )
  )
  mu0 <- list(
    estimate = 0.4,
    n = 3L,
    influence_blocks = list(
      t = list(influence = c(-1, 0, 1), weight = 1)
    )
  )

  expect_error(
    RoCE:::.combine_comparison_arms(
      mu1, mu0, method = "target_only", variance_method = "analytic"
    ),
    "inconsistent lengths"
  )
})
