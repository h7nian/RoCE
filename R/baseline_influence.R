# Conditional empirical derivatives for the fitted comparison estimators.
# Lambda and the lasso active set are held fixed. These calculations do not
# establish uniform inference through variable selection or CV transitions.

BASELINE_INFERENCE_POLICY <- "conditional_active_set_v2"

# Solve only in directions identified by the training design. The coordinate
# change below is for differentiation, not a refit of the penalized model.
# A score or evaluation direction outside this space is an unsupported
# derivative, not a reason to silently truncate a singular value.
.solve_baseline_curvature <- function(design, curvature, sensitivity, scores,
                                      checks = list(), lambda = NA_real_,
                                      context = "Baseline nuisance derivative") {
  design <- as.matrix(design)
  n <- nrow(design)
  p <- ncol(design)
  rank <- NA_integer_
  fail <- function(reason) {
    stop(sprintf("%s: %s [active_columns=%d, rank=%s, lambda=%s].",
      context, reason, p, as.character(rank), format(lambda, digits = 12)), call. = FALSE)
  }
  if (!n || !p || length(curvature) != n || length(sensitivity) != p ||
      any(!is.finite(design)) || any(design[, 1L] != 1) ||
      any(!is.finite(curvature)) || any(curvature < 0) ||
      any(!is.finite(sensitivity))) fail("invalid curvature inputs")
  if (!length(scores) || is.null(names(scores)) || anyDuplicated(names(scores))) {
    fail("score blocks must be named")
  }
  centers <- colMeans(design)
  centers[1L] <- 0
  centered <- sweep(design, 2L, centers, "-")
  scales <- sqrt(colMeans(centered^2))
  scales[scales == 0] <- 1
  if (any(!is.finite(scales))) fail("non-finite feature scales")
  scaled_design <- sweep(centered, 2L, scales, "/")
  # Z = X A; vectors transform by A' and score rows by A. In particular,
  # intercept centering must also be applied to the penalty derivatives.
  transform_rows <- function(rows) {
    rows <- as.matrix(rows)
    if (ncol(rows) != p || any(!is.finite(rows))) fail("invalid score/check block")
    sweep(rows - outer(rows[, 1L], centers), 2L, scales, "/")
  }
  scaled_scores <- lapply(scores, transform_rows)
  scaled_checks <- lapply(checks, transform_rows)
  scaled_sensitivity <- drop(transform_rows(matrix(sensitivity, nrow = 1L)))
  decomposition <- svd(scaled_design / sqrt(n), nu = 0L, nv = p)
  relative_values <- c(decomposition$d, rep(0, p - length(decomposition$d))) /
    max(decomposition$d)
  rank_tolerance <- 100 * .Machine$double.eps * max(1L, p)
  curvature_tolerance <- 1e-7
  compatibility_tolerance <- 1e-8
  retained <- relative_values > rank_tolerance
  rank <- sum(retained)
  if (any(relative_values[retained] < curvature_tolerance)) {
    fail("near-collinear design is not numerical redundancy")
  }
  basis <- decomposition$v[, retained, drop = FALSE]
  compatibility <- numeric()
  if (rank < p) {
    null_basis <- decomposition$v[, !retained, drop = FALSE]
    null_fraction <- function(rows) {
      magnitude <- max(abs(rows))
      if (magnitude == 0) return(0)
      rows <- rows / magnitude
      sqrt(sum((rows %*% null_basis)^2) / sum(rows^2))
    }
    compatibility <- c(sensitivity = null_fraction(matrix(scaled_sensitivity, nrow = 1L)),
      vapply(scaled_scores, null_fraction, numeric(1L)),
      vapply(scaled_checks, null_fraction, numeric(1L)))
    incompatible <- which(compatibility > compatibility_tolerance)
    if (length(incompatible)) {
      fail(paste0("non-identifiable derivative; null-space incompatibility in ",
        paste(sprintf("%s=%.3g", names(compatibility)[incompatible],
          compatibility[incompatible]), collapse = ", ")))
    }
  }
  weighted_design <- scaled_design * sqrt(curvature / n)
  weighted <- svd(weighted_design %*% basis, nu = 0L, nv = rank)
  curvature_ratio <- min(weighted$d) / max(weighted$d)
  if (!is.finite(curvature_ratio) || curvature_ratio < curvature_tolerance ||
      min(weighted$d) <= 0) fail("unstable curvature in the identifiable space")

  # Retain the original arithmetic for a well-conditioned, full-rank fit.
  hessian <- crossprod(design, design * curvature) / n
  raw_rcond <- if (all(is.finite(hessian))) rcond(hessian) else 0
  direction <- if (rank == p && is.finite(raw_rcond) && raw_rcond >= 1e-12) {
    tryCatch(solve(hessian, sensitivity), error = function(error) NULL)
  } else NULL
  if (!is.null(direction)) {
    corrections <- lapply(scores, function(score) drop(score %*% direction))
    solver <- "direct"
  } else {
    weighted_basis <- basis %*% weighted$v
    scaled_direction <- drop(weighted_basis %*%
      (drop(crossprod(weighted_basis, scaled_sensitivity)) / weighted$d^2))
    corrections <- lapply(scaled_scores, function(score) drop(score %*% scaled_direction))
    solver <- if (rank == p) "scaled" else "identifiable_subspace"
  }
  if (any(vapply(corrections, function(x) any(!is.finite(x)), logical(1L)))) {
    fail("non-finite derivative correction")
  }
  list(corrections = corrections, diagnostics = list(
    solver = solver, active_columns = p, rank = rank, lambda = lambda,
    raw_rcond = raw_rcond, curvature_ratio = curvature_ratio,
    rank_tolerance = rank_tolerance, compatibility_tolerance = compatibility_tolerance,
    null_compatibility = compatibility))
}

.baseline_glm_adjustment <- function(x, y, raw_prediction, sensitivity, model = NULL,
                                     evaluation_design = NULL, context = "GLM") {
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
  score <- selected * (y - raw_prediction)
  score <- sweep(score, 2L, colMeans(score), "-")
  checks <- list()
  if (!is.null(evaluation_design)) {
    checks$evaluation_design <- evaluation_design[, active, drop = FALSE]
  }
  if (identical(model$type, "glmnet") && length(active) > 1L) {
    # glmnet standardizes x using its empirical population variance. Changing
    # case masses therefore also changes lambda * sd(x_j) * |beta_j| on the
    # original coefficient scale, even with the selected lambda held fixed.
    features <- x[, active[-1L] - 1L, drop = FALSE]
    centered <- sweep(features, 2L, colMeans(features), "-")
    variance <- colMeans(centered^2)
    if (any(variance <= 0)) stop("An active baseline feature has zero variance.", call. = FALSE)
    scale_derivative <- sweep(sweep(centered^2, 2L, variance, "-"), 2L, 2 * sqrt(variance), "/")
    penalty <- model$lambda * sign(model$coefficients[active[-1L]])
    penalty_score <- cbind(0, sweep(scale_derivative, 2L, penalty, "*"))
    checks$data_score <- score
    checks$penalty_score <- penalty_score
    checks$penalty_gradient <- matrix(c(0, sqrt(variance) * penalty), nrow = 1L)
    score <- score - penalty_score
  }
  solved <- .solve_baseline_curvature(selected, derivative, sensitivity[active],
    scores = list(score = score), checks = checks, lambda = model$lambda %||% NA_real_,
    context = paste("Baseline nuisance derivative", context))
  result <- solved$corrections$score
  attr(result, "curvature_diagnostics") <- solved$diagnostics
  result
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
  outcome_adjustment <- .baseline_glm_adjustment(
    x[selected, , drop = FALSE], y[selected], raw_outcome[selected], sensitivity, outcome_model,
    evaluation_design = design, context = paste0("OR arm=", A_val))
  outcome_correction[selected] <- n / sum(selected) * as.numeric(outcome_adjustment)
  propensity_inside <- raw_propensity > PROP_SCORE_LOWER & raw_propensity < PROP_SCORE_UPPER
  propensity_slope <- (if (A_val == 1L) -1 else 1) * selected * (y - m_hat) /
    probability^2 * raw_propensity * (1 - raw_propensity) * propensity_inside
  propensity_sensitivity <- colMeans(design * w * propensity_slope) / normalizer
  propensity_adjustment <- .baseline_glm_adjustment(x, a, raw_propensity,
    propensity_sensitivity, propensity_model, context = paste0("PS arm=", A_val))
  propensity_correction <- as.numeric(propensity_adjustment)
  direct <- w * (phi - estimate) / normalizer
  influence <- direct + outcome_correction + propensity_correction
  list(estimate = estimate, influence = influence, variance = mean(influence^2) / n,
       phi = phi, X_int = design, m_hat = m_hat, pi_hat = pi_hat, d = normalizer,
       outcome_correction = outcome_correction, propensity_correction = propensity_correction,
       curvature_diagnostics = list(OR = attr(outcome_adjustment, "curvature_diagnostics"),
         PS = attr(propensity_adjustment, "curvature_diagnostics")),
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
  source_score <- source_active * tilt
  source_score <- sweep(source_score, 2L, colMeans(source_score), "-")
  target_score <- sweep(target_active, 2L, colMeans(target_active), "-")
  penalty <- c(0, model$lambda * sign(beta[active[-1L]]))
  solved <- .solve_baseline_curvature(source_active, tilt, sensitivity,
    scores = list(source_score = source_score, target_score = target_score),
    checks = list(target_design = target_active, penalty_gradient = matrix(penalty, nrow = 1L)),
    lambda = model$lambda, context = "Density nuisance derivative")
  # Differentiating the empirical normalization of exp(-Z beta) contributes
  # directly through source masses as well as through the fitted beta.
  normalization_correction <- -source_fraction / normalizer *
    mean(inside * normalized * (phi - estimate)) * (normalized - 1)
  list(source = solved$corrections$source_score + normalization_correction,
       target = -solved$corrections$target_score, active = active,
       curvature_diagnostics = solved$diagnostics)
}
