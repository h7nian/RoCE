# Delta-method weight layer of the aggregation variance
# (R/aggregation_weight_influence.R; main.tex sec:adaptive_aggregation,
# supplemental.tex supp:weight_layer).

.weight_layer_fixture <- function(K, seed, screening_rule = "soft_penalty") {
  set.seed(seed)
  data <- generate_simulation_data(
    n_total = 120 * (K + 1), K = K, p = 3, config = "C1",
    dgp_type = "face",
    n_target = 120, n_source_sizes = rep(120, K),
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  fit <- run_tate_crossfit(
    data_split = data_split,
    n_folds = 3,
    communication_mode = "one_round",
    family = "gaussian",
    nlambda_init = 5,
    n_cores = 1,
    verbose = FALSE,
    screening_rule = screening_rule
  )
  list(fit = fit, sizes = c(data_split$t$n, vapply(
    setdiff(names(data_split), "t"), function(s) data_split[[s]]$n, numeric(1L)
  )))
}

# The estimator as a functional of observation masses: weights re-solved by the
# fitted rule on the mass-weighted inner-fold moments, then the mass-weighted
# mean of the fold pseudo-values.
.weight_layer_functional <- function(fit, sizes, mass) {
  inner <- fit$intermediates$inner_fold_info
  K <- length(sizes) - 1L
  weights <- do.call(rbind, lapply(seq_along(inner), function(k) {
    moments <- RoCE:::.inner_fold_moments(RoCE:::.inner_fold_records(inner[[k]]), mass)$moments
    if (identical(fit$aggregation_screening_rule, "quadratic_bias")) {
      return(RoCE:::.quadratic_bias_weights(moments))
    }
    as.numeric(optimize_weights(
      moments$avg_source_est,
      list(V_ot = moments$V_ot, V_t = moments$V_t, V_s = moments$V_s),
      moments$C_ot, list(n_t = moments$n_t, n_s = moments$n_s),
      fit$fold_lambdas[k], moments$avg_target_est, moments$C_cross,
      clip_weights = FALSE
    ))
  }))
  phi <- RoCE:::.compute_phase3_all_phi(
    n_folds = fit$n_folds, fold_weights = weights,
    fold_info = fit$intermediates$fold_info,
    n_t = sizes[1L], n_source_full = sizes[-1L], N_all = sum(sizes),
    source_sites = paste0("s", seq_len(K)), K = K, verbose = FALSE
  )
  ends <- cumsum(sizes)
  sum(vapply(seq_along(sizes), function(s) {
    idx <- (ends[s] - sizes[s] + 1L):ends[s]
    values <- phi[idx] * sizes[s] / sum(sizes)
    sum(mass[[s]] * values) / sum(mass[[s]])
  }, numeric(1L)))
}

test_that("weight-layer variance decomposes and its fixed part is the pseudo-value variance", {
  skip_on_cran()
  fixture <- .weight_layer_fixture(K = 2L, seed = 7101)
  fit <- fixture$fit
  sizes <- fixture$sizes

  expect_equal(
    fit$variance_fixed_weights,
    RoCE:::.multisite_pseudovalue_variance(fit$all_phi_agg, as.integer(sizes)),
    tolerance = 1e-14
  )
  expect_equal(fit$se_fixed_weights, sqrt(fit$variance_fixed_weights), tolerance = 1e-14)
  expect_equal(fit$se, sqrt(fit$variance), tolerance = 1e-14)
  expect_equal(
    fit$variance,
    fit$variance_fixed_weights + fit$weight_layer$indirect_variance +
      fit$weight_layer$cross_term,
    tolerance = 1e-12
  )
  expect_identical(fit$weight_layer$fold_source_cells, length(fit$fold_weights))
  expect_gte(fit$weight_layer$indirect_variance, 0)
  expect_identical(fit$estimate, mean(fit$all_phi_agg))
})

.expect_weight_layer_matches_finite_differences <- function(fixture) {
  fit <- fixture$fit
  sizes <- fixture$sizes
  gradient <- RoCE:::.weight_layer_gradient(
    fit$intermediates$fold_info, fit$intermediates$inner_fold_info,
    fit$fold_weights, fit$fold_lambdas, fit$aggregation_screening_rule,
    fit$fold_weight_psd_ridge, sizes[1L], sizes[-1L]
  )
  ends <- cumsum(sizes)
  direct <- lapply(seq_along(sizes), function(s) {
    phi <- fit$all_phi_agg[(ends[s] - sizes[s] + 1L):ends[s]]
    (phi - mean(phi)) / sum(sizes)
  })
  total <- Map(`+`, direct, c(list(gradient$target), gradient$source))
  expect_lt(max(abs(vapply(total, sum, numeric(1L)))), 1e-10)

  unit <- lapply(sizes, function(n) rep(1, n))
  expect_equal(.weight_layer_functional(fit, sizes, unit), fit$estimate, tolerance = 1e-12)
  epsilon <- 1e-5
  for (direction_id in 1:3) {
    direction <- lapply(seq_along(sizes), function(s) sin(seq_len(sizes[s]) * direction_id + s))
    plus <- Map(function(m, d) m + epsilon * d, unit, direction)
    minus <- Map(function(m, d) m - epsilon * d, unit, direction)
    numerical <- (.weight_layer_functional(fit, sizes, plus) -
                    .weight_layer_functional(fit, sizes, minus)) / (2 * epsilon)
    analytical <- sum(unlist(total) * unlist(direction))
    expect_lt(abs(numerical - analytical), 1e-8 * max(1, abs(analytical)))
  }
}

test_that("weight-layer gradient matches finite differences of the mass functional", {
  skip_on_cran()
  .expect_weight_layer_matches_finite_differences(.weight_layer_fixture(K = 2L, seed = 7102))
})

test_that("quadratic-bias weight layer matches finite differences and has no kinks", {
  skip_on_cran()
  fixture <- .weight_layer_fixture(K = 2L, seed = 7105, screening_rule = "quadratic_bias")
  expect_identical(fixture$fit$weight_layer$kink_cells, 0L)
  expect_true(all(fixture$fit$fold_weight_psd_ridge == 0))
  .expect_weight_layer_matches_finite_differences(fixture)
})

test_that("quadratic-bias weights solve their normal equations", {
  skip_on_cran()
  fit <- .weight_layer_fixture(K = 2L, seed = 7106)$fit
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  mass <- lapply(sizes, function(n) rep(1, n))
  power <- RoCE:::AGG_QUADRATIC_BIAS_POWER
  fold_moments <- lapply(fit$intermediates$inner_fold_info, function(inner) {
    RoCE:::.inner_fold_moments(RoCE:::.inner_fold_records(inner), mass)$moments
  })
  for (moments in fold_moments) {
    form <- RoCE:::.variance_quadratic_form(moments, psd_ridge = 0)
    fold_discrepancy <- moments$avg_source_est - moments$avg_target_est
    fold_weights <- RoCE:::.quadratic_bias_weights(moments)
    curvature <- form$Q + diag(moments$n_t^(power - 1) * fold_discrepancy^2, length(fold_weights))
    expect_lt(max(abs(curvature %*% fold_weights + form$l)), 1e-12 * max(1, max(abs(form$l))))
  }
  moments <- fold_moments[[1L]]
  discrepancy <- moments$avg_source_est - moments$avg_target_est
  weights <- RoCE:::.quadratic_bias_weights(moments)
  # Without discrepancies the rule is the unpenalized variance-optimal weight
  # (coordinate descent at WEIGHT_OPT_TOL, hence the 1e-8 tolerance).
  aligned <- moments
  aligned$avg_source_est <- rep(aligned$avg_target_est, length(weights))
  unpenalized <- optimize_weights(
    aligned$avg_source_est,
    list(V_ot = aligned$V_ot, V_t = aligned$V_t, V_s = aligned$V_s),
    aligned$C_ot, list(n_t = aligned$n_t, n_s = aligned$n_s),
    0, aligned$avg_target_est, aligned$C_cross,
    clip_weights = FALSE
  )
  expect_equal(RoCE:::.quadratic_bias_weights(aligned), as.numeric(unpenalized), tolerance = 1e-8)
  # A larger discrepancy shrinks that source's weight toward zero.
  shifted <- moments
  shifted$avg_source_est[1L] <- shifted$avg_target_est + 10 * max(abs(discrepancy), 0.1)
  expect_lt(abs(RoCE:::.quadratic_bias_weights(shifted)[1L]), abs(weights[1L]))
  # Non-finite moments and an indefinite variance quadratic fail loudly.
  broken <- moments
  broken$C_ot[1L] <- NA_real_
  expect_error(RoCE:::.quadratic_bias_weights(broken), "non-finite")
  indefinite <- moments
  indefinite$V_ot <- 100 * max(moments$V_t, moments$V_s)
  expect_error(RoCE:::.quadratic_bias_weights(indefinite), "not positive definite")
})

test_that("weight layer is reported for the arm-specific and hard-threshold paths", {
  skip_on_cran()
  set.seed(7104)
  data <- generate_simulation_data(
    n_total = 360, K = 2, p = 3, config = "C1", dgp_type = "face",
    n_target = 120, n_source_sizes = c(120, 120), warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  arm <- run_crossfit(
    data_split = data_split, n_folds = 3, communication_mode = "one_round",
    family = "gaussian", nlambda_init = 5, n_cores = 1, verbose = FALSE, A_val = 1
  )
  for (fit in list(arm, run_tate_crossfit(
    data_split = data_split, n_folds = 3, communication_mode = "one_round",
    family = "gaussian", nlambda_init = 5, n_cores = 1, verbose = FALSE,
    screening_rule = "hard_threshold"
  ))) {
    expect_true(all(is.finite(c(fit$se, fit$se_fixed_weights))))
    expect_equal(
      fit$variance,
      fit$variance_fixed_weights + fit$weight_layer$indirect_variance +
        fit$weight_layer$cross_term,
      tolerance = 1e-12
    )
  }
})

test_that("weight layer vanishes when no source is retained", {
  skip_on_cran()
  set.seed(7103)
  data <- generate_simulation_data(
    n_total = 240, K = 1, p = 3, config = "C1", dgp_type = "face",
    n_target = 120, n_source_sizes = 120, warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  fit <- run_tate_crossfit(
    data_split = data_split, n_folds = 3, communication_mode = "one_round",
    family = "gaussian", nlambda_init = 5, n_cores = 1, verbose = FALSE,
    lambda_selection = 1e6
  )
  expect_true(all(abs(fit$fold_weights) < 1e-12))
  expect_equal(fit$weight_layer$indirect_variance, 0)
  expect_equal(fit$variance, fit$variance_fixed_weights, tolerance = 1e-14)
  expect_equal(fit$variance, fit$target_only$variance, tolerance = 1e-12)
})
