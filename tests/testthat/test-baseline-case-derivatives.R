test_that("baseline nuisance derivatives match independent case-weight refits", {
  set.seed(1139)
  n <- 1000L
  x <- matrix(stats::rnorm(n * 5L), n, 5L)
  x[, 3L] <- .3 * x[, 3L]
  x[, 4L] <- 2 * x[, 4L]
  a <- stats::rbinom(n, 1L, stats::plogis(-.2 + 3 * x[, 1L]))
  outcomes <- list(
    binomial = stats::rbinom(n, 1L, stats::plogis(.2 + 3.5 * x[, 2L] + .5 * x[, 1L]^2)),
    gaussian = .2 + .3 * a + x[, 2L] + .5 * x[, 1L]^2 + stats::rnorm(n))
  transport <- exp(.3 * x[, 4L])
  for (family in names(outcomes)) for (arm in 0:1) for (lambda in c(0, .02)) {
    y <- outcomes[[family]]
    fit <- .baseline_case_fit(x, y, a, arm, family, rep(1, n), lambda, transport)
    result <- calculate_aipw_influence(y, a, x, fit$outcome$raw_prediction,
      fit$propensity$raw_prediction, transport, arm, family,
      outcome_model = fit$outcome, propensity_model = fit$propensity)
    expect_equal(result$estimate, fit$estimate, tolerance = 1e-12)
    expect_lt(abs(mean(result$influence)), 1e-10)
    for (index in c(1L, 117L, 396L, 810L)) {
      upper <- lower <- rep(1, n)
      step <- 1e-3
      upper[index] <- 1 + step
      lower[index] <- 1 - step
      plus <- .baseline_case_fit(x, y, a, arm, family, upper, lambda, transport)
      minus <- .baseline_case_fit(x, y, a, arm, family, lower, lambda, transport)
      for (model in c("outcome", "propensity")) {
        expect_identical(sign(plus[[model]]$coefficients), sign(minus[[model]]$coefficients))
      }
      numerical <- n * (plus$estimate - minus$estimate) / (2 * step)
      expect_equal(result$influence[index], numerical, tolerance = 2e-5,
                   info = paste(family, arm, lambda, index))
    }
  }
})

test_that("retaining baseline metadata preserves predictions and the RNG stream", {
  set.seed(1141)
  x <- matrix(stats::rnorm(1000 * 4L), 1000, 4L)
  y <- stats::rbinom(1000, 1L, stats::plogis(x[, 1L]))
  set.seed(1142)
  original <- fit_glmnet_cv(x, y, x, clip_fn = clip_propensity)
  original_rng <- .Random.seed
  set.seed(1142)
  retained <- fit_glmnet_cv(x, y, x, clip_fn = clip_propensity, retain_model = TRUE)
  expect_identical(as.numeric(original), as.numeric(retained))
  expect_identical(original_rng, .Random.seed)
  model <- attr(retained, "nuisance_model")
  expect_identical(model$type, "glmnet")
  expect_equal(model$raw_prediction, stats::plogis(drop(cbind(1, x) %*% model$coefficients)),
               tolerance = 1e-12)
})

test_that("constant outcome fallback uses its empirical mean derivative", {
  set.seed(1143)
  n <- 1000L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  a <- rep(0:1, length.out = n)
  y <- numeric(n)
  y[which(a == 1L)[1:5]] <- 1
  model <- attr(fit_glmnet_cv(x[a == 1L, ], y[a == 1L], x,
    on_degenerate_response = "constant", retain_model = TRUE), "nuisance_model")
  result <- calculate_aipw_influence(y, a, x, model$raw_prediction, rep(.5, n),
    w = exp(.2 * x[, 1L]), outcome_model = model)
  expect_identical(model$type, "constant")
  expect_true(all(is.finite(result$influence)))
  expect_lt(abs(mean(result$influence)), 1e-12)
})

.refit_density_case <- function(source_design, target_design, source_mass, target_mass,
                                 model, source_phi, target_phi = NULL) {
  active <- if (model$lambda == 0) seq_along(model$coefficients) else
    c(1L, which(model$coefficients[-1L] != 0) + 1L)
  source <- source_design[, active, drop = FALSE]
  target <- target_design[, active, drop = FALSE]
  beta <- model$coefficients[active]
  mean_target <- colSums(target * target_mass) / sum(target_mass)
  penalty <- c(0, model$lambda * sign(beta[-1L]))
  objective <- function(beta) sum(mean_target * beta) +
    sum(source_mass * exp(-drop(source %*% beta))) / sum(source_mass) + sum(penalty * beta)
  for (iteration in 1:100) {
    tilt <- exp(-drop(source %*% beta))
    gradient <- mean_target - colSums(source * source_mass * tilt) / sum(source_mass) + penalty
    if (max(abs(gradient)) < 1e-12) break
    hessian <- crossprod(source, source * source_mass * tilt) / sum(source_mass)
    direction <- drop(solve(hessian, gradient))
    step <- 1
    while (objective(beta - step * direction) > objective(beta) -
           1e-4 * step * sum(gradient * direction) + 1e-14) step <- step / 2
    beta <- beta - step * direction
  }
  stopifnot(max(abs(gradient)) < 1e-10)
  tilt <- exp(-drop(source %*% beta))
  normalized <- tilt / weighted.mean(tilt, source_mass)
  bounded <- pmin(10, pmax(.1, normalized))
  numerator <- weighted.mean(bounded * source_phi, source_mass)
  denominator <- weighted.mean(bounded, source_mass)
  if (!is.null(target_phi)) {
    numerator <- numerator + weighted.mean(target_phi, target_mass)
    denominator <- denominator + 1
  }
  list(estimate = numerator / denominator, weights = bounded)
}

test_that("density derivatives include normalization, clipping and pooled sensitivities", {
  set.seed(1147)
  n <- 1000L
  source <- matrix(stats::rnorm(n * 3L), n, 3L)
  target <- matrix(stats::rnorm(n * 3L), n, 3L)
  target[, 1L] <- target[, 1L] + 1.3
  source_phi <- sin(source[, 1L]) + source[, 2L]^2
  target_phi <- .7 + .2 * target[, 1L] + sin(target[, 2L])
  for (lambda in c(0, .05)) {
    weights <- suppressWarnings(calculate_dr_weights(source, target, lambda = lambda))
    model <- attr(weights, "density_model")
    expect_gt(sum(weights == .1 | weights == 10), 0)
    for (pooled in c(FALSE, TRUE)) {
      source_fraction <- if (pooled) .5 else 1
      normalizer <- source_fraction * mean(weights) + (1 - source_fraction)
      estimate <- (source_fraction * mean(weights * source_phi) +
        (1 - source_fraction) * mean(target_phi)) / normalizer
      correction <- .baseline_density_influence(source, target, weights, source_phi,
        estimate, normalizer, source_fraction)
      source_direct <- source_fraction / normalizer * weights * (source_phi - estimate)
      source_influence <- source_direct - mean(source_direct) + correction$source
      target_influence <- (1 - source_fraction) / normalizer *
        (target_phi - mean(target_phi)) + correction$target
      for (site in c("source", "target")) for (index in c(1L, 214L, 531L, 971L)) {
        plus <- minus <- list(source = rep(1, n), target = rep(1, n))
        step <- 1e-3
        plus[[site]][index] <- 1 + step
        minus[[site]][index] <- 1 - step
        evaluated <- vapply(list(plus, minus), function(mass) {
          .refit_density_case(cbind(1, source), cbind(1, target), mass$source, mass$target,
            model, source_phi, if (pooled) target_phi else NULL)$estimate
        }, numeric(1L))
        numerical <- n * diff(rev(evaluated)) / (2 * step)
        actual <- if (site == "source") source_influence[index] else target_influence[index]
        expect_equal(actual, numerical, tolerance = 2e-5,
                     info = paste(lambda, pooled, site, index))
      }
    }
  }
})

test_that("pooled DR influence matches a complete fixed-lambda case-weight refit", {
  set.seed(1151)
  n <- 1000L
  make_site <- function(shift) {
    x <- matrix(stats::rnorm(n * 3L), n, 3L)
    x[, 1L] <- x[, 1L] + shift
    a <- stats::rbinom(n, 1L, stats::plogis(.2 + x[, 1L]))
    y <- stats::rbinom(n, 1L, stats::plogis(.3 * a + x[, 2L] + .5 * x[, 1L]^2 + shift))
    list(W_outcome = x, Z_site = x, A = a, Y = y, n = n)
  }
  data <- list(t = make_site(1), s1 = make_site(0))
  weights <- suppressWarnings(calculate_dr_weights(data$s1$Z_site, data$t$Z_site, lambda = .05))
  x <- rbind(data$t$W_outcome, data$s1$W_outcome)
  y <- c(data$t$Y, data$s1$Y)
  a <- c(data$t$A, data$s1$A)
  for (arm in 0:1) {
    pooled <- estimate_pooled_dr(data, dr_weights_by_site = list(s1 = weights),
      A_val = arm, variance_method = "analytic")
    base <- with_seed(pooled$components$nuisance_seed, calculate_weighted_site_aipw(
      y, a, x, weights = c(rep(1, n), weights), A_val = arm))
    expect_equal(pooled$estimate, base$estimate, tolerance = 1e-12)
    lambdas <- c(base$propensity_model$lambda, base$outcome_model$lambda)
    for (site in c("t", "s1")) for (index in c(93L, 687L)) {
      upper <- lower <- list(t = rep(1, n), s1 = rep(1, n))
      step <- 1e-3
      upper[[site]][index] <- 1 + step
      lower[[site]][index] <- 1 - step
      estimates <- vapply(list(upper, lower), function(mass) {
        mass <- lapply(mass, function(w) n * w / sum(w))
        transported <- .refit_density_case(cbind(1, data$s1$Z_site), cbind(1, data$t$Z_site),
          mass$s1, mass$t, attr(weights, "density_model"), rep(0, n))$weights
        .baseline_case_fit(x, y, a, arm, "binomial", c(mass$t, mass$s1), lambdas,
          c(rep(1, n), transported))$estimate
      }, numeric(1L))
      numerical <- n * (estimates[1L] - estimates[2L]) / (2 * step)
      block <- pooled$influence_blocks[[site]]
      actual <- block$weight * (block$influence[index] - mean(block$influence))
      expect_equal(actual, numerical, tolerance = 5e-5, info = paste(arm, site, index))
    }
  }
})
