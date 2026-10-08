# Independent fixed-lambda case-mass refits used by derivative regression tests.
.baseline_case_fit <- function(x, y, a, arm, family, mass, lambda, transport) {
  selected <- a == arm
  design <- cbind(1, x)
  fit_model <- function(rows, response, family, lambda) {
    if (lambda == 0) {
      fit <- stats::glm.fit(design[rows, , drop = FALSE], response,
        family = if (family == "binomial") stats::binomial() else stats::gaussian(),
        weights = mass[rows], control = stats::glm.control(epsilon = 1e-12, maxit = 100L))
      beta <- fit$coefficients
      type <- "mle"
    } else {
      fit <- glmnet::glmnet(x[rows, , drop = FALSE], response, family = family,
        weights = mass[rows], lambda = lambda, alpha = 1, standardize = TRUE,
        thresh = 1e-14, maxit = 100000L)
      beta <- as.numeric(stats::coef(fit, s = lambda))
      type <- "glmnet"
    }
    predictor <- drop(design %*% beta)
    raw <- if (family == "binomial") stats::plogis(predictor) else predictor
    list(type = type, family = family, coefficients = beta, lambda = lambda, raw_prediction = raw)
  }
  lambda <- rep(lambda, length.out = 2L)
  propensity <- fit_model(rep(TRUE, length(a)), a, "binomial", lambda[1L])
  outcome <- fit_model(selected, y[selected], family, lambda[2L])
  clipped_ps <- pmin(.99, pmax(.01, propensity$raw_prediction))
  clipped_or <- if (family == "binomial") pmin(.999, pmax(.001, outcome$raw_prediction)) else
    outcome$raw_prediction
  probability <- if (arm == 1L) clipped_ps else 1 - clipped_ps
  phi <- clipped_or + selected * (y - clipped_or) / probability
  list(estimate = sum(mass * transport * phi) / sum(mass * transport),
       outcome = outcome, propensity = propensity)
}
