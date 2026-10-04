.score_calibration_fixture <- function(counts, probability, outcome, fractional = FALSE) {
  support <- c(-3, -1, 0, 1, 3)
  rows <- list()
  for (index in seq_along(support)) for (arm in 0:1) {
    size <- counts[index] * if (arm == 1L) probability[index] else 1 - probability[index]
    stopifnot(abs(size - round(size)) < 1e-10)
    size <- as.integer(round(size))
    response <- if (fractional) rep(outcome[index], size) else {
      successes <- size * outcome[index]
      stopifnot(abs(successes - round(successes)) < 1e-10)
      c(rep(1, round(successes)), rep(0, size - round(successes)))
    }
    rows[[length(rows) + 1L]] <- data.frame(x = rep(support[index], size), A = arm, Y = response)
  }
  rows <- do.call(rbind, rows)
  data <- list(W_outcome = matrix(rows$x, ncol = 1L), Z_site = matrix(rows$x, ncol = 1L),
               A = rows$A, Y = rows$Y, n = nrow(rows))
  # Each block contains both treatment groups; initial plug-ins are fixed
  # population coefficients for this independent loss/score check.
  indices <- list(fold_3 = seq.int(1L, data$n, 2L), fold_4 = seq.int(2L, data$n, 2L))
  blocks <- lapply(indices, function(index) list(W_outcome = data$W_outcome[index, , drop = FALSE],
    Z_site = data$Z_site[index, , drop = FALSE], A = data$A[index], Y = data$Y[index], n = length(index)))
  list(data = data, blocks = blocks)
}

.score_calibration_sensitivity <- function(fit, source, target_design, arm, radius, target = FALSE) {
  source_design <- cbind(1, source$W_outcome)
  weight_design <- cbind(1, source$Z_site)
  mean <- plogis(drop(source_design %*% fit$outcome))
  predictor <- drop(weight_design %*% fit$weight)
  tilt <- exp(-pmax(-radius, pmin(radius, predictor)))
  selected <- source$A == arm
  derivative_weight <- -colMeans(weight_design * selected * tilt *
    (abs(predictor) < radius) * (source$Y - mean))
  derivative_outcome <- if (target) {
    colMeans(source_design * (1 - selected * (1 + tilt)) * mean * (1 - mean))
  } else {
    target_mean <- plogis(drop(target_design %*% fit$outcome))
    colMeans(target_design * target_mean * (1 - target_mean)) -
      colMeans(source_design * selected * tilt * mean * (1 - mean))
  }
  max(abs(c(derivative_weight, derivative_outcome)))
}

test_that("calibration matches the clipped score in both model-correctness branches", {
  support <- c(-3, -1, 0, 1, 3)
  for (site in c("target", "source")) for (branch in c("outcome", "weight")) {
    if (branch == "outcome") {
      counts <- c(100, 150, 150, 200, 400)
      probability <- c(.2, .4, .6, .7, .8)
      beta <- c(1.2, 2)
      truth <- plogis(beta[1] + beta[2] * support)
      gamma <- c(0, 0)
      radius <- 5
      target_counts <- if (site == "target") counts else c(200, 200, 300, 200, 100)
      # Fractional responses encode conditional expected Bernoulli losses;
      # this is a population derivative check, not an outcome simulation.
      fixture <- .score_calibration_fixture(counts, probability, truth, fractional = TRUE)
    } else {
      radius <- log(4)
      beta <- c(.4, .2)
      truth <- c(.8, .6, .3, .5, .9)
      if (site == "target") {
        counts <- rep(200, 5)
        probability <- c(.2, .25, .5, .75, .8)
        gamma <- c(0, log(3))
        target_counts <- counts
      } else {
        counts <- c(320, 200, 200, 200, 80)
        probability <- rep(.5, 5)
        gamma <- c(0, log(2))
        target_counts <- c(640, 200, 100, 50, 10)
      }
      fixture <- .score_calibration_fixture(counts, probability, truth)
    }
    data <- fixture$data
    target_design <- cbind(1, rep(support, target_counts))
    fits <- lapply(c("legacy", "score_derivative"), function(recipe) {
      derivative_radius <- if (recipe == "legacy") radius else Inf
      if (site == "target") {
        selected <- data$A != 1L
        moment <- mean(selected) * .mean_glm_gradient_site_basis(
          data$W_outcome[selected, , drop = FALSE], data$Z_site[selected, , drop = FALSE],
          beta, 1L, 1L, derivative_radius)
      } else {
        moment <- .mean_glm_gradient_site_basis(target_design[, -1L, drop = FALSE],
          target_design[, -1L, drop = FALSE], beta, 1L, 1L, derivative_radius)
      }
      .fit_fold_summed_calibration(fixture$blocks,
        setNames(list(beta, beta), names(fixture$blocks)),
        setNames(list(gamma, gamma), names(fixture$blocks)), moment,
        site = site, training_folds = 3:4, A_val = 1L, family_int = 1L, link_int = 1L,
        M_tau = radius, nlambda = 5L, max_iter = 5000L, lambda_rule = "min",
        weight_lambda = 0, outcome_lambda = 0, tol = 1e-10,
        calibration_recipe = recipe)
    })
    sensitivities <- vapply(fits, .score_calibration_sensitivity, numeric(1L),
      source = data, target_design = target_design, arm = 1L, radius = radius, target = site == "target")
    expect_true(sensitivities[1L] > 1e-6, info = paste(site, branch))
    expect_true(sensitivities[2L] < 1e-8, info = paste(site, branch))
    expect_identical(fits[[2L]]$calibration_recipe, "score_derivative")
  }
})

test_that("calibrated target propensity initialization uses the opposite-arm moment", {
  fixture <- .score_calibration_fixture(rep(200, 5), c(.2, .25, .5, .75, .8), c(.8, .6, .3, .5, .9))
  data <- fixture$data
  for (arm in 0:1) {
    fit <- .fit_initial_target_propensity(data$Z_site, data$A, arm, 8L, "min",
      initialization = "calibrated", M_tau = log(4), tol = 1e-10)
    design <- cbind(1, data$Z_site)
    tilt <- exp(-pmax(-log(4), pmin(log(4), drop(design %*% fit))))
    gradient <- colMeans(design * ((data$A != arm) - (data$A == arm) * tilt))
    penalty <- attr(fit, "lambda_used")
    residual <- c(gradient[1L], ifelse(fit[-1L] != 0,
      gradient[-1L] + penalty * sign(fit[-1L]), pmax(abs(gradient[-1L]) - penalty, 0)))
    expect_lt(max(abs(residual)), 1e-8)
  }
})

test_that("calibration settings reject ambiguous or inactive options", {
  expect_error(.validate_calibration_control(list(unknown = TRUE)), "named list", fixed = TRUE)
  expect_error(.validate_calibration_control(list(target_radius = 2), "lasso"), "require", fixed = TRUE)
  expect_error(.validate_calibration_control(list(recipe = c("legacy", "score_derivative"))),
               "one character value", fixed = TRUE)
  expect_error(validate_truncation_parameters(5, -Inf), "positive", fixed = TRUE)
})

test_that("the score-derivative recipe agrees with legacy fitting when truncation is inactive", {
  fixture <- .score_calibration_fixture(rep(200, 5), c(.2, .25, .5, .75, .8), c(.8, .6, .3, .5, .9))
  beta <- c(.1, .05)
  gamma <- c(-.7, .01)
  data <- fixture$data
  opposite <- data$A == 0L
  moment <- mean(opposite) * .mean_glm_gradient_site_basis(
    data$W_outcome[opposite, , drop = FALSE], data$Z_site[opposite, , drop = FALSE], beta, 1L, 1L, 5)
  fits <- lapply(c("legacy", "score_derivative"), function(recipe) {
    .fit_fold_summed_calibration(fixture$blocks,
      setNames(list(beta, beta), names(fixture$blocks)),
      setNames(list(gamma, gamma), names(fixture$blocks)), moment,
      "t", 3:4, 1L, 1L, 1L, 5, 8L, 10000L, "min",
      tol = 1e-10, calibration_recipe = recipe)
  })
  for (field in c("weight", "outcome")) {
    expect_equal(as.numeric(fits[[1L]][[field]]), as.numeric(fits[[2L]][[field]]), tolerance = 1e-12)
    expect_identical(attr(fits[[1L]][[field]], "lambda_used"), attr(fits[[2L]][[field]], "lambda_used"))
  }
})

test_that("three-level fitting and reaggregation preserve the requested target model", {
  skip_on_cran()
  set.seed(1103)
  data <- split_data_by_site(generate_simulation_data(n_total = 2000L, K = 1L, p = 4L,
    n_target = 1000L, n_source_sizes = 1000L, config = "C2", dgp_type = "face",
    outcome_type = "binary", warn_ignored = FALSE))
  control <- list(recipe = "score_derivative", target_propensity_initialization = "calibrated",
                  target_radius = log(9))
  expected_control <- c(control, list(source_nuisance_method = "calibrated"))
  fit <- run_tate_crossfit(data, n_folds = 4L, communication_mode = "one_round",
    target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
    calibration_control = control, calibration_layout = "compact", nuisance_solver = "proximal_newton",
    nlambda_init = 8L, n_cores = 1L, verbose = FALSE)
  expect_identical(fit$crossfit_levels, 3L)
  expect_identical(fit$calibration_control, expected_control)
  expect_true(all(is.finite(c(fit$estimate, fit$se))))
  changed_radius <- reaggregate_tate_crossfit(data, fit, M_tau_inference = 4,
                                             aggregation_mode = "joint_tate")
  expect_identical(changed_radius$calibration_control, expected_control)
  for (arm in c("mu1", "mu0")) for (outer in seq_len(fit$n_folds)) {
    before <- fit$arm_results[[arm]]$fold_results[[outer]]
    after <- changed_radius$arm_results[[arm]]$fold_results[[outer]]
    expect_equal(before$target_only$arm_propensity, after$target_only$arm_propensity, tolerance = 1e-14)
    expect_gte(min(before$target_only$arm_propensity), .1 - 1e-14)
    expect_lte(max(before$target_only$arm_propensity), .9 + 1e-14)
    expect_identical(before$target_only$calibration$calibration_recipe, "score_derivative")
    for (inner in names(before$target_only_inner)) {
      expect_identical(before$target_only_inner[[inner]]$calibration$propensity_initialization, "calibrated")
      expect_identical(before$source_results$s1$inner_calibrated[[inner]]$calibration_recipe, "score_derivative")
    }
  }
})
