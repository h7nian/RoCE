# Conditional empirical derivatives for the fitted comparison estimators.
# Lambda and the lasso active set are held fixed. These calculations do not
# establish uniform inference through variable selection or CV transitions.

BASELINE_INFERENCE_POLICY <- "conditional_active_set_v1"

.baseline_glm_adjustment <- function(x, y, raw_prediction, sensitivity, model = NULL) {
  x <- as.matrix(x)
  n <- nrow(x)
  design <- cbind(1, x)
  family <- model$family %||% "binomial"
  if (identical(model$type, "constant")) {
    return(as.numeric(sensitivity) * (y - mean(y)))
  }
  active <- if (is.null(model$coefficients)) seq_len(ncol(design)) else
    c(1L, which(model$coefficients[-1L] != 0) + 1L)
  selected <- design[, active, drop = FALSE]
  derivative <- if (family == "gaussian") rep(1, n) else raw_prediction * (1 - raw_prediction)
  hessian <- crossprod(selected, selected * derivative) / n
  score <- selected * (y - raw_prediction)
  score <- sweep(score, 2L, colMeans(score), "-")
  if (identical(model$type, "glmnet") && length(active) > 1L) {
    # glmnet standardizes x using its empirical population variance. Changing
    # case masses therefore also changes lambda * sd(x_j) * |beta_j| on the
    # original coefficient scale, even with the selected lambda held fixed.
    features <- x[, active[-1L] - 1L, drop = FALSE]
    centered <- sweep(features, 2L, colMeans(features), "-")
    variance <- colMeans(centered^2)
    if (any(variance <= 0)) stop("An active baseline feature has zero variance.", call. = FALSE)
    scale_derivative <- sweep(sweep(centered^2, 2L, variance, "-"), 2L, 2 * sqrt(variance), "/")
    score[, -1L] <- score[, -1L, drop = FALSE] - sweep(scale_derivative, 2L,
      model$lambda * sign(model$coefficients[active[-1L]]), "*")
  }
  direction <- tryCatch(solve(hessian, sensitivity[active]), error = function(error) {
    stop("Baseline nuisance derivative has singular active-set curvature: ",
         conditionMessage(error), call. = FALSE)
  })
  drop(score %*% direction)
}

.baseline_aipw_influence <- function(y, a, x, m_hat, pi_hat, w, A_val, family,
                                      outcome_model = NULL, propensity_model = NULL) {
  n <- length(y)
  design <- cbind(1, as.matrix(x))
  raw_outcome <- outcome_model$raw_prediction %||% as.numeric(m_hat)
  raw_propensity <- propensity_model$raw_prediction %||% as.numeric(pi_hat)
  outcome_model <- outcome_model %||% list(type = "mle", family = family)
  propensity_model <- propensity_model %||% list(type = "mle", family = "binomial")
  m_hat <- clip_outcome_pred(raw_outcome, family)
  pi_hat <- clip_propensity(raw_propensity)
  selected <- a == A_val
  probability <- if (A_val == 1L) pi_hat else 1 - pi_hat
  phi <- m_hat + selected * (y - m_hat) / probability
  normalizer <- mean(w)
  estimate <- mean(w * phi) / normalizer
  outcome_inside <- if (family == "binomial") {
    raw_outcome > OUTCOME_PRED_LOWER & raw_outcome < OUTCOME_PRED_UPPER
  } else rep(TRUE, n)
  outcome_slope <- w * (1 - selected / probability) * outcome_inside / normalizer
  outcome_correction <- numeric(n)
  if (identical(outcome_model$type, "constant")) {
    sensitivity <- mean(outcome_slope)
  } else {
    derivative <- if (family == "gaussian") rep(1, n) else raw_outcome * (1 - raw_outcome)
    sensitivity <- colMeans(design * outcome_slope * derivative)
  }
  outcome_correction[selected] <- n / sum(selected) * .baseline_glm_adjustment(
    x[selected, , drop = FALSE], y[selected], raw_outcome[selected], sensitivity, outcome_model)
  propensity_inside <- raw_propensity > PROP_SCORE_LOWER & raw_propensity < PROP_SCORE_UPPER
  propensity_slope <- (if (A_val == 1L) -1 else 1) * selected * (y - m_hat) /
    probability^2 * raw_propensity * (1 - raw_propensity) * propensity_inside
  propensity_sensitivity <- colMeans(design * w * propensity_slope) / normalizer
  propensity_correction <- .baseline_glm_adjustment(x, a, raw_propensity,
    propensity_sensitivity, propensity_model)
  direct <- w * (phi - estimate) / normalizer
  influence <- direct + outcome_correction + propensity_correction
  list(estimate = estimate, influence = influence, variance = mean(influence^2) / n,
       phi = phi, X_int = design, m_hat = m_hat, pi_hat = pi_hat, d = normalizer,
       outcome_correction = outcome_correction, propensity_correction = propensity_correction,
       inference_scope = "conditional_nuisance_lambda_and_active_sets")
}

.baseline_density_influence <- function(Z_source, Z_target, weights, phi,
                                         estimate, normalizer, source_fraction = 1) {
  model <- attr(weights, "density_model")
  if (is.null(model)) {
    stop("Density-weight derivatives require weights returned by calculate_dr_weights().",
         call. = FALSE)
  }
  source_design <- cbind(1, as.matrix(Z_source))
  target_design <- cbind(1, as.matrix(Z_target))
  beta <- model$coefficients
  if (length(beta) != ncol(source_design) || length(phi) != nrow(source_design)) {
    stop("Density nuisance metadata and evaluation rows are not conformable.", call. = FALSE)
  }
  tilt <- exp(-drop(source_design %*% beta))
  if (any(!is.finite(tilt)) || any(tilt <= WEIGHT_MIN | tilt >= WEIGHT_MAX)) {
    stop("Density nuisance reached a numerical tilt guard; its derivative is not supported.", call. = FALSE)
  }
  normalized <- tilt / mean(tilt)
  bounded <- pmax(DR_WEIGHT_LOWER, pmin(DR_WEIGHT_UPPER, normalized))
  if (!isTRUE(all.equal(as.numeric(weights), bounded, tolerance = 1e-10))) {
    stop("Density nuisance metadata do not reproduce the supplied weights.", call. = FALSE)
  }
  active <- if (model$lambda == 0) seq_along(beta) else c(1L, which(beta[-1L] != 0) + 1L)
  source_active <- source_design[, active, drop = FALSE]
  target_active <- target_design[, active, drop = FALSE]
  inside <- normalized > DR_WEIGHT_LOWER & normalized < DR_WEIGHT_UPPER
  tilted_mean <- colMeans(source_active * normalized)
  normalized_design <- sweep(source_active, 2L, tilted_mean, "-")
  sensitivity <- -source_fraction / normalizer *
    colMeans(normalized_design * inside * normalized * (phi - estimate))
  hessian <- crossprod(source_active, source_active * tilt) / nrow(source_design)
  direction <- tryCatch(solve(hessian, sensitivity), error = function(error) {
    stop("Density nuisance derivative has singular active-set curvature: ",
         conditionMessage(error), call. = FALSE)
  })
  source_score <- source_active * tilt
  source_score <- sweep(source_score, 2L, colMeans(source_score), "-")
  target_score <- sweep(target_active, 2L, colMeans(target_active), "-")
  # Differentiating the empirical normalization of exp(-Z beta) contributes
  # directly through source masses as well as through the fitted beta.
  normalization_correction <- -source_fraction / normalizer *
    mean(inside * normalized * (phi - estimate)) * (normalized - 1)
  list(source = drop(source_score %*% direction) + normalization_correction,
       target = -drop(target_score %*% direction), active = active)
}
