.joint_test_records <- function() {
  set.seed(994)
  target <- lapply(c(70L, 84L), function(n) {
    b <- rnorm(n)
    cbind(b, .8 - .7 * b + rnorm(n, sd = .4), -.4 - .5 * b + rnorm(n, sd = .5))
  })
  source <- lapply(c(60L, 66L), function(n) {
    a <- rep(0:1, length.out = n)
    cbind(0, a * (.6 + rnorm(n, sd = .3)), -(1 - a) * (.4 + rnorm(n, sd = .3)))
  })
  list(blocks = list(target, source),
       ids = list(list(1:70, 71:154), list(1:60, 61:126)))
}

test_that("common weights preserve zero empirical variance components", {
  weights <- optimize_weights(.5, list(V_ot = .6, V_t = 0, V_s = .4),
    C_ot = 0, n_samples = list(n_t = 500, n_s = 1000), lambda = 0, mu_ot = .5)
  expect_equal(as.numeric(weights), .75, tolerance = 1e-12)
  expect_lte(attr(weights, "kkt_residual"), attr(weights, "kkt_threshold"))
})

test_that("joint moments retain same-source covariance and count each site once", {
  records <- .joint_test_records()
  moments <- .joint_inner_moments(records)
  expect_identical(moments$N_all, 280L)
  centered_source <- do.call(rbind, lapply(records$blocks[[2L]], function(x) {
    sweep(x, 2L, colMeans(x), "-")
  }))
  expected_source <- crossprod(centered_source) / nrow(centered_source)
  expect_equal(moments$sites[[2L]]$covariance, expected_source, tolerance = 1e-14)
  expect_gt(expected_source[2L, 3L], .02)
  projection <- matrix(c(1, 1), 2L, 1L)
  common_variance <- sum(vapply(moments$sites, function(site) {
    d <- rowSums(site$centered[, -1L, drop = FALSE])
    mean(d^2) / site$n
  }, numeric(1L)))
  expect_equal(drop(t(projection) %*% moments$Q %*% projection), common_variance, tolerance = 1e-14)
})

test_that("joint QP agrees with an independent smooth optimum and satisfies L1 KKT", {
  moments <- .joint_inner_moments(.joint_test_records())
  unpenalized <- .solve_joint_weights(moments, 0, "soft_penalty")
  expect_equal(unpenalized$weights, as.numeric(solve(moments$Q, -moments$l)), tolerance = 1e-12)
  fit <- .solve_joint_weights(moments, .2, "soft_penalty")
  gradient <- 2 * moments$N_all * drop(moments$Q %*% fit$weights + moments$l)
  residual <- ifelse(fit$weights != 0,
    gradient + fit$penalty * sign(fit$weights), pmax(abs(gradient) - fit$penalty, 0))
  expect_lt(max(abs(residual)), 1e-7)
  common <- optimize(function(eta) {
    e <- c(eta, eta)
    moments$N_all * (moments$constant + 2 * sum(moments$l * e) +
      drop(t(e) %*% moments$Q %*% e)) + sum(fit$penalty * abs(e))
  }, c(-5, 5))
  expect_lte(fit$objective, common$objective + 1e-8)
})

test_that("joint empirical-mass derivatives include both arm coordinates at each site", {
  records <- .joint_test_records()
  moments <- .joint_inner_moments(records)
  increment <- c(.3, -.2)
  mass <- list(rep(1, 154), rep(1, 126))
  step <- 1e-5
  for (rule in c("soft_penalty", "hard_threshold", "quadratic_bias")) {
    lambda <- if (rule == "hard_threshold") .03 else .2
    fit <- .solve_joint_weights(moments, lambda, rule)
    derivative <- .joint_weight_derivative(moments, fit, rule, increment)
    for (site in 1:2) {
      expect_lt(abs(sum(derivative$sites[[site]])), 1e-12)
      for (row in c(1L, 49L, 80L)) {
        upper <- lower <- mass
        upper[[site]][row] <- 1 + step
        lower[[site]][row] <- 1 - step
        objective <- function(mass) {
          sum(increment * .solve_joint_weights(.joint_inner_moments(records, mass), lambda, rule)$weights)
        }
        numerical <- (objective(upper) - objective(lower)) / (2 * step)
        expect_lt(abs(derivative$sites[[site]][row] - numerical), 1e-8)
      }
    }
  }
})

test_that("all TATE aggregation modes reuse the same fitted nuisance models", {
  skip_on_cran()
  set.seed(1004)
  data <- split_data_by_site(generate_simulation_data(
    n_total = 2000L, n_target = 1000L, n_source_sizes = 1000L,
    K = 1L, p = 4L, config = "C1", dgp_type = "face",
    outcome_type = "continuous", warn_ignored = FALSE))
  fit <- suppressWarnings(run_tate_crossfit(data, n_folds = 4L,
    communication_mode = "one_round", nlambda_init = 6L, family = "gaussian",
    n_cores = 1L, verbose = FALSE, target_nuisance_method = "hou_calibrated",
    source_validation_method = "calibrated", calibration_layout = "compact"))
  estimates <- list(common_tate = fit)
  expect_identical(fit$crossfit_levels, 3L)
  for (mode in c("common_tate", "separate_arms", "joint_tate")) {
    reaggregated <- suppressWarnings(reaggregate_tate_crossfit(
      data, fit, M_tau_inference = fit$M_tau_inference, aggregation_mode = mode))
    expect_identical(reaggregated$aggregation_mode, mode)
    expect_identical(reaggregated$crossfit_levels, 3L)
    expect_true(is.finite(reaggregated$estimate))
    expect_true(is.finite(reaggregated$se))
    expect_identical(reaggregated$variance_method, "leading")
    expect_equal(reaggregated$variance,
      .multisite_pseudovalue_variance(reaggregated$all_phi_tau,
        as.integer(vapply(data, function(x) x$n, numeric(1L)))), tolerance = 1e-12)
    expect_equal(reaggregated$se, reaggregated$se_fixed_weights, tolerance = 1e-12)
    expect_equal(reaggregated$variance_weight_adjusted,
      reaggregated$variance + reaggregated$weight_layer$indirect_variance +
        reaggregated$weight_layer$cross_term, tolerance = 1e-12)
    expect_equal(reaggregated$se_weight_adjusted^2, reaggregated$variance_weight_adjusted,
      tolerance = 1e-12)
    expect_equal(c(reaggregated$ci_lower, reaggregated$ci_upper),
      reaggregated$estimate + c(-1, 1) * qnorm(.975) * reaggregated$se_fixed_weights,
      tolerance = 1e-12)
    expect_identical(reaggregated$source_validation_method, "calibrated")
    expect_identical(reaggregated$arm_results$mu1$fold_results[[1L]]$source_results$s1$alpha_ts,
                     fit$arm_results$mu1$fold_results[[1L]]$source_results$s1$alpha_ts)
    expect_identical(reaggregated$arm_results$mu0$fold_results[[1L]]$target_only$calibration$outcome,
                     fit$arm_results$mu0$fold_results[[1L]]$target_only$calibration$outcome)
    estimates[[mode]] <- reaggregated
  }
  # A cutoff-only sensitivity scan must not silently switch the estimator's
  # aggregation objective on rows whose lambda differs from the original.
  for (mode in names(estimates)) {
    original <- estimates[[mode]]
    lambdas <- original$aggregation_lambda_selection * c(1, .5, 2)
    sensitivity <- reaggregate_tate_sensitivity_grid(data, original,
      data.frame(cutoff = 1 / lambdas, M_tau_inference = original$M_tau_inference))
    expect_identical(sensitivity$results[[1L]], original)
    expect_identical(sensitivity$n_inference_refreshes, 0L)
    expect_identical(sensitivity$n_nuisance_refits, 0L)
    legacy <- original
    legacy$variance <- original$variance_weight_adjusted
    legacy$se <- original$se_weight_adjusted
    legacy$variance_method <- legacy$variance_weight_adjusted <- legacy$se_weight_adjusted <- NULL
    legacy$ci_lower <- legacy$estimate - qnorm(.975) * legacy$se
    legacy$ci_upper <- legacy$estimate + qnorm(.975) * legacy$se
    refreshed <- reaggregate_tate_sensitivity_grid(data, legacy,
      data.frame(cutoff = 1/lambdas[1L], M_tau_inference = original$M_tau_inference))
    expect_identical(refreshed$results[[1L]]$variance_method, "leading")
    expect_equal(refreshed$results[[1L]]$estimate, original$estimate, tolerance = 1e-12)
    expect_equal(refreshed$results[[1L]]$fold_weights, original$fold_weights, tolerance = 1e-12)
    expect_equal(refreshed$results[[1L]]$se, original$se_fixed_weights, tolerance = 1e-12)
    expect_identical(refreshed$n_nuisance_refits, 0L)
    for (i in seq_along(lambdas)) {
      candidate <- sensitivity$results[[i]]
      direct <- calculate_tate_crossfit_aggregation(data,
        original$arm_results$mu1, original$arm_results$mu0,
        lambda_selection = lambdas[i], aggregation_mode = mode,
        screening_rule = original$aggregation_screening_rule, verbose = FALSE)
      expect_identical(candidate$aggregation_mode, mode)
      expect_equal(candidate$estimate, direct$estimate, tolerance = 1e-12)
      expect_equal(candidate$se, direct$se, tolerance = 1e-12)
      expect_equal(candidate$fold_weights, direct$fold_weights, tolerance = 1e-12)
    }
  }
  expect_equal(estimates$common_tate$estimate, fit$estimate, tolerance = 1e-12)
  expect_equal(estimates$common_tate$se, fit$se, tolerance = 1e-12)
  separate <- estimates$separate_arms
  expect_equal(separate$estimate,
               separate$arm_results$mu1$estimate - separate$arm_results$mu0$estimate, tolerance = 1e-12)
  expect_identical(dim(estimates$joint_tate$fold_weights), c(4L, 2L))
  for (outer in seq_len(4L)) {
    expect_length(fit$arm_results$mu1$fold_results[[outer]]$source_results$s1$inner_calibrated, 3L)
    expect_equal(estimates$joint_tate$intermediates$joint_moments[[outer]]$N_all,
      2000 - length(fit$arm_results$mu1$intermediates$fold_info[[outer]]$target_idx) -
      length(fit$arm_results$mu1$intermediates$fold_info[[outer]]$source_idx[[1L]]))
  }
})
