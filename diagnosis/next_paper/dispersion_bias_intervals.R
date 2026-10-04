# Known-Gaussian reference with explicit bounds on pooled-candidate bias.
# This module does not assume that fitted invalid RoCE candidates are regular.

noncentrality_upper_bound <- function(statistic, degrees_freedom, failure_probability, tolerance = 1e-9) {
  if (!is.numeric(statistic) || !length(statistic) || any(!is.finite(statistic)) || any(statistic < 0) ||
      length(degrees_freedom) != 1L || !is.finite(degrees_freedom) || degrees_freedom < 1 ||
      degrees_freedom != floor(degrees_freedom) || length(failure_probability) != 1L ||
      !is.finite(failure_probability) || failure_probability <= 0 || failure_probability >= 1 ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Noncentrality inversion requires nonnegative statistics, positive integer df and a failure probability in (0,1)")
  }
  log_probability <- log(failure_probability)
  log_cdf <- function(values, noncentrality) {
    result <- pchisq(values, degrees_freedom, ncp = noncentrality, log.p = TRUE)
    if (anyNA(result) || any(result > 0)) stop("Noncentral chi-square CDF evaluation failed")
    result
  }
  upper <- numeric(length(statistic))
  active <- which(log_cdf(statistic, 0) > log_probability)
  if (!length(active)) return(upper)
  values <- statistic[active]
  lower <- numeric(length(active))
  high <- pmax(1, values - degrees_freedom + 10 * sqrt(degrees_freedom + values) + 20)
  for (iteration in seq_len(100L)) {
    enlarge <- log_cdf(values, high) > log_probability
    if (!any(enlarge)) break
    high[enlarge] <- high[enlarge] * 2
    if (any(!is.finite(high))) stop("Noncentrality upper bound exceeded finite arithmetic")
  }
  if (any(log_cdf(values, high) > log_probability)) stop("Could not bracket the noncentrality upper bound")
  for (iteration in seq_len(100L)) {
    unfinished <- which(high - lower > tolerance * (1 + high))
    if (!length(unfinished)) break
    midpoint <- lower[unfinished] + (high[unfinished] - lower[unfinished]) / 2
    above <- log_cdf(values[unfinished], midpoint) > log_probability
    lower[unfinished[above]] <- midpoint[above]
    high[unfinished[!above]] <- midpoint[!above]
  }
  if (any(high - lower > tolerance * (1 + high))) stop("Noncentrality inversion did not converge")
  # Retain the conservative endpoint of the root bracket, not its midpoint.
  upper[active] <- high
  upper
}

checked_dispersion_covariance <- function(covariance) {
  covariance <- as.matrix(covariance)
  if (!is.numeric(covariance) || nrow(covariance) < 1L || nrow(covariance) != ncol(covariance) ||
      any(!is.finite(covariance)) || max(abs(covariance - t(covariance))) > 1e-12 * max(abs(covariance))) {
    stop("A finite symmetric covariance matrix is required")
  }
  covariance <- (covariance + t(covariance)) / 2
  factor <- tryCatch(chol(covariance), error = function(error) NULL)
  if (is.null(factor)) stop("The dispersion pivot requires a positive-definite covariance")
  list(covariance = covariance, inverse = chol2inv(factor))
}

make_dispersion_calibration <- function(covariance, valid_minimum, known_valid = integer(),
                                        bound_method = c("subset_exact", "uniform_upper"), max_subsets = 50000L) {
  checked <- checked_dispersion_covariance(covariance)
  covariance <- checked$covariance
  count <- nrow(covariance)
  if (length(valid_minimum) != 1L || !is.finite(valid_minimum) || valid_minimum < 1 ||
      valid_minimum > count || valid_minimum != floor(valid_minimum) || !is.numeric(known_valid) ||
      any(!is.finite(known_valid)) || any(known_valid < 1 | known_valid > count | known_valid != floor(known_valid)) ||
      anyDuplicated(known_valid) || length(known_valid) > valid_minimum) stop("Invalid guaranteed or known valid-candidate counts")
  bound_method <- match.arg(bound_method)
  precision <- rowSums(checked$inverse)
  variance <- 1 / sum(precision)
  weights <- variance * precision
  subset_variance <- function(indices) {
    factor <- chol(covariance[indices, indices, drop = FALSE])
    1 / sum(chol2inv(factor))
  }
  examined <- 0L
  if (valid_minimum == count) {
    bias_factor <- 0
  } else if (bound_method == "subset_exact") {
    if (length(max_subsets) != 1L || !is.finite(max_subsets) || max_subsets < 1 || max_subsets != floor(max_subsets)) {
      stop("max_subsets must be a positive integer")
    }
    additional <- valid_minimum - length(known_valid)
    optional <- setdiff(seq_len(count), known_valid)
    if (lchoose(length(optional), additional) > log(max_subsets) + 1e-12) {
      stop("Exact subset enumeration exceeds max_subsets; select uniform_upper explicitly")
    }
    subsets <- if (additional == 0L) list(integer()) else combn(optional, additional, simplify = FALSE)
    values <- vapply(subsets, function(subset) subset_variance(c(known_valid, subset)), numeric(1L))
    bias_factor <- max(values) - variance
    examined <- length(subsets)
  } else {
    bias_factor <- max(diag(covariance)) / valid_minimum + (valid_minimum - 1) / valid_minimum *
      max(covariance[row(covariance) != col(covariance)]) - variance
    if (length(known_valid)) bias_factor <- min(bias_factor, subset_variance(known_valid) - variance)
  }
  if (!is.finite(bias_factor) || bias_factor < -1e-10 * max(diag(covariance))) stop("Invalid quadratic bias factor")
  list(covariance = covariance, inverse = checked$inverse, weights = weights, variance = variance,
       bias_factor = max(0, bias_factor), degrees_freedom = count - 1L,
       valid_minimum = valid_minimum, known_valid = as.integer(known_valid),
       bound_method = bound_method, subsets_examined = examined)
}

dispersion_arm_statistics <- function(estimates, calibration, failure_probability) {
  estimates <- as.matrix(estimates)
  if (!is.numeric(estimates) || ncol(estimates) != length(calibration$weights) ||
      any(!is.finite(estimates))) stop("Candidate draws must match the calibrated covariance dimension")
  means <- drop(estimates %*% calibration$weights)
  if (calibration$bias_factor == 0) {
    return(list(mean = means, statistic = rep(0, nrow(estimates)),
      noncentrality_upper = rep(0, nrow(estimates)), bias_allowance = rep(0, nrow(estimates))))
  }
  centered <- sweep(estimates, 1L, means, "-")
  statistic <- rowSums((centered %*% calibration$inverse) * centered)
  if (any(statistic < -1e-10)) stop("The dispersion quadratic form is negative")
  upper <- noncentrality_upper_bound(pmax(statistic, 0), calibration$degrees_freedom, failure_probability)
  list(mean = means, statistic = statistic, noncentrality_upper = upper,
       bias_allowance = sqrt(calibration$bias_factor * upper))
}

make_arm_dispersion_calibration <- function(covariance, valid_mu1, valid_mu0,
                                            bound_method = c("subset_exact", "uniform_upper")) {
  checked <- checked_dispersion_covariance(covariance)
  covariance <- checked$covariance
  width <- nrow(covariance) / 2
  if (width != floor(width) || width < 1L ||
      length(valid_mu1) != 1L || length(valid_mu0) != 1L ||
      any(!is.finite(c(valid_mu1, valid_mu0))) || any(c(valid_mu1, valid_mu0) < 0) ||
      any(c(valid_mu1, valid_mu0) > width - 1L) || any(c(valid_mu1, valid_mu0) != floor(c(valid_mu1, valid_mu0)))) {
    stop("Arm validity counts must lie between zero and the source count")
  }
  first <- seq_len(width)
  second <- width + first
  bound_method <- match.arg(bound_method)
  mu1 <- make_dispersion_calibration(covariance[first, first], valid_mu1 + 1L, 1L, bound_method)
  mu0 <- make_dispersion_calibration(covariance[second, second], valid_mu0 + 1L, 1L, bound_method)
  contrast <- c(mu1$weights, -mu0$weights)
  anchor_contrast <- numeric(2L * width); anchor_contrast[c(1L, width + 1L)] <- c(1, -1)
  list(mu1 = mu1, mu0 = mu0, width = width,
       variance = drop(crossprod(contrast, covariance %*% contrast)),
       anchor_variance = drop(crossprod(anchor_contrast, covariance %*% anchor_contrast)),
       anchor_only = valid_mu1 == 0 && valid_mu0 == 0, bound_method = bound_method)
}

arm_dispersion_intervals <- function(estimates, calibration, alpha = .05,
                                      dispersion_fraction = .5, anchor_fraction = NULL) {
  estimates <- as.matrix(estimates)
  if (!is.numeric(estimates) || ncol(estimates) != 2L * calibration$width || any(!is.finite(estimates))) {
    stop("Draw columns must contain all mu1 candidates followed by all mu0 candidates")
  }
  probabilities <- c(alpha, dispersion_fraction, anchor_fraction)
  if (length(alpha) != 1L || length(dispersion_fraction) != 1L ||
      (!is.null(anchor_fraction) && length(anchor_fraction) != 1L) ||
      any(!is.finite(probabilities)) || any(probabilities <= 0 | probabilities >= 1)) stop("Error allocations must lie in (0,1)")
  anchor_estimate <- estimates[, 1L] - estimates[, calibration$width + 1L]
  if (calibration$anchor_only) {
    radius <- qnorm(alpha / 2, lower.tail = FALSE) * sqrt(calibration$anchor_variance)
    return(data.frame(pooled_estimate = anchor_estimate, lower = anchor_estimate - radius, upper = anchor_estimate + radius,
      bias_allowance = 0, used_anchor_fallback = FALSE))
  }
  interval_alpha <- alpha * if (is.null(anchor_fraction)) 1 else (1 - anchor_fraction)
  active <- sum(c(calibration$mu1$bias_factor, calibration$mu0$bias_factor) > 0)
  dispersion_alpha <- if (active) interval_alpha * dispersion_fraction else 0
  noise_alpha <- interval_alpha - dispersion_alpha
  first <- seq_len(calibration$width)
  mu1 <- dispersion_arm_statistics(estimates[, first, drop = FALSE], calibration$mu1,
                                    if (active) dispersion_alpha / active else NULL)
  mu0 <- dispersion_arm_statistics(estimates[, calibration$width + first, drop = FALSE], calibration$mu0,
                                    if (active) dispersion_alpha / active else NULL)
  point <- mu1$mean - mu0$mean
  bias <- mu1$bias_allowance + mu0$bias_allowance
  radius <- qnorm(noise_alpha / 2, lower.tail = FALSE) * sqrt(calibration$variance) + bias
  lower <- point - radius
  upper <- point + radius
  fallback <- rep(FALSE, length(point))
  if (!is.null(anchor_fraction)) {
    anchor_radius <- qnorm(alpha * anchor_fraction / 2, lower.tail = FALSE) * sqrt(calibration$anchor_variance)
    anchor_lower <- anchor_estimate - anchor_radius
    anchor_upper <- anchor_estimate + anchor_radius
    lower <- pmax(lower, anchor_lower); upper <- pmin(upper, anchor_upper)
    fallback <- lower > upper
    lower[fallback] <- anchor_lower[fallback]; upper[fallback] <- anchor_upper[fallback]
  }
  data.frame(pooled_estimate = point, lower = lower, upper = upper, bias_allowance = bias,
             used_anchor_fallback = fallback)
}
