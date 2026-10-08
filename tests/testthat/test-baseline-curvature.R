test_that("full-rank curvature preserves the original solve and scaled fits agree", {
  set.seed(1161)
  x <- matrix(rnorm(600), 200, 3)
  design <- cbind(1, x)
  curvature <- plogis(x[, 1]) * plogis(-x[, 1])
  score <- design * sin(x[, 2])
  score <- sweep(score, 2, colMeans(score), "-")
  sensitivity <- colMeans(design * cos(x[, 3]))
  expected <- drop(score %*% solve(crossprod(design, design * curvature) / nrow(x), sensitivity))
  rng <- .Random.seed
  solved <- .solve_baseline_curvature(design, curvature, sensitivity, list(score = score))
  expect_identical(solved$corrections$score, expected)
  expect_identical(solved$diagnostics$solver, "direct")
  expect_identical(.Random.seed, rng)

  change <- diag(c(1, 1e8, 1, 1))
  change[1, 2] <- 3e8
  scaled <- .solve_baseline_curvature(design %*% change, curvature,
    drop(crossprod(change, sensitivity)), list(score = score %*% change))
  expect_identical(scaled$diagnostics$solver, "scaled")
  expect_equal(scaled$corrections$score, expected, tolerance = 1e-10)
})

test_that("redundant standardized lasso fits match independent case-mass refits", {
  set.seed(1163)
  n <- 1000L
  x <- matrix(rnorm(n * 3), n, 3)
  a <- rbinom(n, 1, plogis(.1 + x[, 1]))
  outcomes <- list(binomial = rbinom(n, 1, plogis(.3 * a + x[, 1] + x[, 2])),
                   gaussian = .3 * a + x[, 1] + x[, 2] + rnorm(n))
  transport <- exp(.2 * x[, 3])
  augmented <- cbind(x, 2 + 3 * x[, 1])
  # Split an identified fitted coefficient between two affine copies. This
  # leaves both the predictor and glmnet's standardized L1 penalty unchanged.
  expand_model <- function(model) {
    beta <- model$coefficients
    expect_gt(abs(beta[2]), .05)
    extra <- beta[2] / 6
    model$coefficients <- c(beta[1] - 2 * extra, beta[2] / 2, beta[-c(1, 2)], extra)
    model
  }
  for (family in names(outcomes)) for (arm in 0:1) {
    y <- outcomes[[family]]
    fit <- .baseline_case_fit(x, y, a, arm, family, rep(1, n), .02, transport)
    outcome <- expand_model(fit$outcome)
    propensity <- expand_model(fit$propensity)
    result <- .baseline_aipw_influence(y, a, augmented, outcome$raw_prediction,
      propensity$raw_prediction, transport, arm, family, outcome, propensity)
    original <- .baseline_aipw_influence(y, a, x, fit$outcome$raw_prediction,
      fit$propensity$raw_prediction, transport, arm, family, fit$outcome, fit$propensity)
    expect_identical(result$estimate, original$estimate)
    expect_equal(result$influence, original$influence, tolerance = 1e-10)
    expect_equal(result$variance, original$variance, tolerance = 1e-10)
    expect_identical(result$curvature_diagnostics$PS$solver, "identifiable_subspace")
    expect_equal(result$curvature_diagnostics$PS$rank,
      1 + sum(fit$propensity$coefficients[-1] != 0))
    expect_lt(max(result$curvature_diagnostics$PS$null_compatibility), 1e-10)
    direction <- rep(c(-1, 1), length.out = n)
    for (step in c(1e-3, 1e-4)) {
      plus <- .baseline_case_fit(augmented, y, a, arm, family, 1 + step * direction, .02, transport)
      minus <- .baseline_case_fit(augmented, y, a, arm, family, 1 - step * direction, .02, transport)
      numerical <- (plus$estimate - minus$estimate) / (2 * step)
      expect_equal(mean(result$influence * direction), numerical, tolerance = 2e-5,
        info = paste(family, arm, step))
    }
  }
})

test_that("unsupported directions and near-collinearity are not silently discarded", {
  set.seed(1165)
  x <- matrix(rnorm(600), 200, 3)
  design <- cbind(1, x[, 1:2], x[, 1])
  score <- design * sin(x[, 2])
  sensitivity <- colMeans(design * cos(x[, 2]))
  expect_error(.solve_baseline_curvature(design, rep(1, 200), sensitivity + c(0, 0, 0, 1),
    list(score = score)), "non-identifiable.*sensitivity")
  broken_score <- score
  broken_score[, 4] <- broken_score[, 4] + x[, 3]
  expect_error(.solve_baseline_curvature(design, rep(1, 200), sensitivity,
    list(score = broken_score)), "non-identifiable.*score")
  design[, 4] <- design[, 4] + 1e-9 * x[, 3]
  expect_error(.solve_baseline_curvature(design, rep(1, 200), colMeans(design),
    list(score = design)), "near-collinear")
  expect_error(.solve_baseline_curvature(cbind(1, x), c(1, rep(0, 199)), rep(1, 4),
    list(score = cbind(1, x))), "unstable curvature")

  # A training-arm dependency need not hold on the evaluation population.
  training <- cbind(x[, 1], x[, 1])
  evaluation <- cbind(1, x[, 1], x[, 1] + x[, 2])
  prediction <- plogis(.2 + .4 * x[, 1])
  expect_error(.baseline_glm_adjustment(training, rep(0:1, 100), prediction,
    colMeans(cbind(1, training)), list(type = "mle", family = "binomial"),
    evaluation_design = evaluation, context = "OR arm=1"), "OR arm=1.*evaluation_design")

  # A general dependency need not preserve the standardized L1 penalty.
  features <- cbind(x[, 1:2], x[, 1] + x[, 2])
  prediction <- plogis(drop(cbind(1, features) %*% rep(.1, 4)))
  expect_error(.baseline_glm_adjustment(features, rep(0:1, 100), prediction,
    colMeans(cbind(1, features)), list(type = "glmnet", family = "binomial",
      lambda = .2, coefficients = rep(.1, 4))), "non-identifiable.*penalty")
})

test_that("density derivatives require the dependency to hold in both populations", {
  set.seed(1167)
  source <- matrix(rnorm(1200), 400, 3)
  target <- matrix(rnorm(1200), 400, 3)
  target[, 1] <- target[, 1] + .3
  phi <- sin(source[, 1]) + source[, 2]^2
  for (lambda in c(0, .02)) {
    weights <- suppressWarnings(calculate_dr_weights(source, target, lambda = lambda))
    estimate <- weighted.mean(phi, weights)
    original <- .baseline_density_influence(source, target, weights, phi, estimate, mean(weights))
    model <- attr(weights, "density_model")
    beta <- model$coefficients
    expect_gt(abs(beta[2]), .02)
    model$coefficients <- c(beta[1], beta[2] / 2, beta[-c(1, 2)], beta[2] / 2)
    augmented_weights <- weights
    attr(augmented_weights, "density_model") <- model
    augmented_source <- cbind(source, source[, 1])
    augmented_target <- cbind(target, target[, 1])
    result <- .baseline_density_influence(augmented_source, augmented_target, augmented_weights,
      phi, estimate, mean(weights))
    expect_equal(result$source, original$source, tolerance = 1e-10)
    expect_equal(result$target, original$target, tolerance = 1e-10)
    expect_identical(result$curvature_diagnostics$solver, "identifiable_subspace")

    # Even a constant shift vanishes from centered target scores, but still
    # violates the target moment and must be rejected by the design check.
    augmented_target[, 4] <- augmented_target[, 4] + .2
    expect_error(.baseline_density_influence(augmented_source, augmented_target, augmented_weights,
      phi, estimate, mean(weights)), "non-identifiable.*target_design")
  }
})

test_that("DR TATE with duplicate features retains paired-arm variance and parallel parity", {
  set.seed(1171)
  n <- 300L
  make_site <- function(shift) {
    x <- matrix(rnorm(n * 4), n, 4)
    x[, 1] <- x[, 1] + shift
    a <- rbinom(n, 1, plogis(.2 + .5 * x[, 1]))
    y <- rbinom(n, 1, plogis(.3 * a + .7 * x[, 1] + .5 * x[, 2]))
    list(W_outcome = x, Z_site = x, A = a, Y = y, n = n)
  }
  data <- list(t = make_site(.5), s1 = make_site(0), s2 = make_site(-.1))
  weights <- lapply(data[-1], function(site) {
    w <- suppressWarnings(calculate_dr_weights(site$Z_site, data$t$Z_site, lambda = .02))
    model <- attr(w, "density_model")
    beta <- model$coefficients
    expect_gt(abs(beta[2]), .02)
    model$coefficients <- c(beta[1], beta[2] / 2, beta[-c(1, 2)], beta[2] / 2)
    attr(w, "density_model") <- model
    w
  })
  data <- lapply(data, function(site) {
    site$W_outcome <- site$Z_site <- cbind(site$Z_site, site$Z_site[, 1])
    site
  })
  evaluate <- function(cores) with_seed(1172, run_all_comparisons_tate(data,
    methods = c("federated_dr", "pooled_dr"), variance_method = "analytic",
    family = "binomial", dr_weights_by_site = weights, n_cores = cores))
  serial <- evaluate(1L)
  for (method in names(serial)) {
    fit <- serial[[method]]
    expect_true(is.finite(fit$estimate) && is.finite(fit$se) && fit$se > 0)
    components <- fit$components
    expect_equal(components$variance_analytic,
      components$mu1_influence_se^2 + components$mu0_influence_se^2 -
        2 * components$cross_arm_covariance, tolerance = 1e-10)
    expect_identical(components$curvature_diagnostics$mu1$s1$density$solver,
      "identifiable_subspace")
    expect_identical(components$curvature_diagnostics$mu0$s1$density$solver,
      "identifiable_subspace")
  }
  if (.Platform$OS.type != "windows" && isTRUE(parallel::detectCores() >= 2L)) {
    parallel_result <- evaluate(2L)
    for (method in names(serial)) {
      expect_equal(parallel_result[[method]]$estimate, serial[[method]]$estimate, tolerance = 1e-12)
      expect_equal(parallel_result[[method]]$se, serial[[method]]$se, tolerance = 1e-12)
    }
  }
})
