.make_variance_formula_data_split <- function(n_by_site) {
  make_site <- function(n, offset) {
    x <- seq_len(n)
    list(
      W_outcome = cbind(x1 = (x - mean(x)) / n, x2 = sin(x + offset)),
      Z_site = cbind(z1 = cos(x + offset), z2 = (x %% 5L) / 5),
      A = rep(c(1L, 0L), length.out = n),
      Y = rep(c(1L, 0L, 1L, 1L), length.out = n),
      n = n
    )
  }
  Map(make_site, n_by_site, seq_along(n_by_site))
}

.baseline_site_fits <- function(estimates, variances, n_by_site) {
  Map(function(est, var, n) {
    list(
      estimate = est,
      variance = var,
      se = sqrt(var),
      n = n,
      V_ot = var * n,
      varphi_ot = rep(0, n)
    )
  }, estimates, variances, n_by_site)
}

.baseline_scalar <- function(x) {
  unname(as.numeric(x))
}

test_that("sample-size baseline variance matches its random-effects formula", {
  n_by_site <- c(t = 100L, s1 = 80L, s2 = 120L)
  data_split <- .make_variance_formula_data_split(n_by_site)
  estimates <- c(t = 0.50, s1 = 0.68, s2 = 0.42)
  variances <- c(t = 0.010, s1 = 0.018, s2 = 0.014)
  site_fits <- .baseline_site_fits(estimates, variances, n_by_site)

  res <- estimate_sample_size_weighted(
    data_split,
    family = "binomial",
    use_crossfit = FALSE,
    site_fits = site_fits,
    variance_method = "analytic"
  )

  total_n <- sum(n_by_site)
  weights <- n_by_site / total_n
  expected_estimate <- sum(weights * estimates)
  het <- calculate_dl_heterogeneity(estimates, variances)
  expected_variance <- sum(weights^2 * (variances + het$tau_sq))

  expect_equal(.baseline_scalar(res$estimate), .baseline_scalar(expected_estimate),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$variance), .baseline_scalar(expected_variance),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$se), .baseline_scalar(sqrt(res$variance)),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$components$var_random_effects),
               .baseline_scalar(expected_variance),
               tolerance = 1e-12)
})

test_that("inverse-variance baseline variance matches the documented max(RE, pooled) formula", {
  n_by_site <- c(t = 90L, s1 = 110L, s2 = 130L)
  data_split <- .make_variance_formula_data_split(n_by_site)
  estimates <- c(t = 0.45, s1 = 0.70, s2 = 0.36)
  variances <- c(t = 0.012, s1 = 0.020, s2 = 0.016)
  site_fits <- .baseline_site_fits(estimates, variances, n_by_site)

  res <- estimate_inverse_variance_weighted(
    data_split,
    family = "binomial",
    use_crossfit = FALSE,
    site_fits = site_fits,
    variance_method = "analytic"
  )

  precisions <- 1 / variances
  fe_estimate <- sum(precisions * estimates) / sum(precisions)
  fe_variance <- 1 / sum(precisions)
  het <- calculate_dl_heterogeneity(estimates, variances)

  if (het$tau_sq > 0) {
    re_precisions <- 1 / (variances + het$tau_sq)
    formula_estimate <- sum(re_precisions * estimates) / sum(re_precisions)
    formula_variance <- 1 / sum(re_precisions)
    weights <- re_precisions / sum(re_precisions)
  } else {
    formula_estimate <- fe_estimate
    formula_variance <- fe_variance
    weights <- precisions / sum(precisions)
  }

  within_site_var <- sum(weights^2 * variances)
  between_site_var <- sum(weights^2 * (estimates - formula_estimate)^2)
  pooled_variance <- within_site_var + between_site_var
  expected_variance <- max(formula_variance, pooled_variance)

  expect_equal(.baseline_scalar(res$estimate), .baseline_scalar(formula_estimate),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$variance), .baseline_scalar(expected_variance),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$components$var_pooled),
               .baseline_scalar(pooled_variance), tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$components$within_site_var),
               .baseline_scalar(within_site_var), tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$components$between_site_var),
               .baseline_scalar(between_site_var), tolerance = 1e-12)
  expect_equal(.baseline_scalar(res$se), .baseline_scalar(sqrt(res$variance)),
               tolerance = 1e-12)
})

test_that("baseline variance formulas reduce to fixed-effect formulas when tau is zero", {
  n_by_site <- c(t = 80L, s1 = 100L, s2 = 120L)
  data_split <- .make_variance_formula_data_split(n_by_site)
  estimates <- c(t = 0.55, s1 = 0.55, s2 = 0.55)
  variances <- c(t = 0.010, s1 = 0.020, s2 = 0.015)
  site_fits <- .baseline_site_fits(estimates, variances, n_by_site)

  ss <- estimate_sample_size_weighted(
    data_split,
    family = "binomial",
    use_crossfit = FALSE,
    site_fits = site_fits,
    variance_method = "analytic"
  )
  ss_weights <- n_by_site / sum(n_by_site)
  ss_expected_variance <- sum(ss_weights^2 * variances)

  expect_equal(.baseline_scalar(ss$components$tau_squared), 0, tolerance = 1e-12)
  expect_equal(.baseline_scalar(ss$variance), .baseline_scalar(ss_expected_variance),
               tolerance = 1e-12)

  ivw <- estimate_inverse_variance_weighted(
    data_split,
    family = "binomial",
    use_crossfit = FALSE,
    site_fits = site_fits,
    variance_method = "analytic"
  )
  ivw_expected_variance <- 1 / sum(1 / variances)

  expect_equal(.baseline_scalar(ivw$components$tau_squared), 0, tolerance = 1e-12)
  expect_equal(.baseline_scalar(ivw$estimate), .baseline_scalar(estimates[[1L]]),
               tolerance = 1e-12)
  expect_equal(.baseline_scalar(ivw$variance),
               .baseline_scalar(ivw_expected_variance), tolerance = 1e-12)
  expect_equal(.baseline_scalar(ivw$components$var_pooled),
               .baseline_scalar(ivw_expected_variance), tolerance = 1e-12)
})

test_that("DR baseline variances equal their exposed influence-function components", {
  skip_if_not(exists("make_small_data_split"), message = "test helper data factory not loaded")

  data_split <- make_small_data_split(
    n_per_site = 160L,
    K = 2L,
    p = 4L,
    seed = 9042L,
    n_folds = 3L,
    outcome_type = "binary"
  )

  tilted <- estimate_tilted_aipw(data_split, family = "binomial",
                                 variance_method = "analytic")
  expect_equal(
    tilted$variance,
    tilted$components$target_variance_component +
      tilted$components$source_variance_component,
    tolerance = 1e-12
  )
  expect_equal(tilted$se, sqrt(tilted$variance), tolerance = 1e-12)

  federated <- estimate_federated_dr(data_split, family = "binomial",
                                     variance_method = "analytic")
  expect_equal(
    federated$variance,
    federated$components$var_target_component +
      federated$components$var_source_component,
    tolerance = 1e-12
  )
  expect_equal(federated$se, sqrt(federated$variance), tolerance = 1e-12)

  pooled <- estimate_pooled_dr(data_split, family = "binomial",
                               variance_method = "analytic")
  expected_pooled_variance <- max(
    pooled$components$variance_wss / (pooled$n^2),
    VARIANCE_MIN
  )
  expect_equal(pooled$variance, expected_pooled_variance, tolerance = 1e-12)
  expect_equal(pooled$se, sqrt(pooled$variance), tolerance = 1e-12)
})
