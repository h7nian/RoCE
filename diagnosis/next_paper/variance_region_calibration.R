# Gaussian reference only: source source_confidence_sets.R before this file.
# Coverage requires a simultaneous variance region with the stated failure
# probability; nuisance-fitted score variances do not automatically supply it.

check_variance_region <- function(shared_bounds, private_bounds) {
  if (!is.numeric(shared_bounds) || length(shared_bounds) != 2L ||
      any(!is.finite(shared_bounds)) || shared_bounds[1L] < 0 ||
      shared_bounds[1L] > shared_bounds[2L]) {
    stop("shared_bounds must contain finite ordered nonnegative bounds")
  }
  if (!is.matrix(private_bounds) || !is.numeric(private_bounds) ||
      ncol(private_bounds) != 2L || nrow(private_bounds) == 0L ||
      any(!is.finite(private_bounds)) || any(private_bounds[, 1L] < 0) ||
      any(private_bounds[, 2L] <= 0) ||
      any(private_bounds[, 1L] > private_bounds[, 2L]) ||
      any(!is.finite(shared_bounds[2L] + private_bounds[, 2L]))) {
    stop("private_bounds must have ordered nonnegative lower and positive upper bounds")
  }
  invisible(NULL)
}

normal_interval_probability <- function(critical, location, scale) {
  if (length(critical) != 1L || !is.finite(critical) || critical < 0 ||
      length(location) != length(scale) || any(!is.finite(c(location, scale))) ||
      any(scale < 0)) stop("Invalid normal interval parameters")
  probability <- as.numeric(abs(location) <= critical)
  positive <- scale > 0
  probability[positive] <- pnorm((critical - location[positive]) / scale[positive]) -
    pnorm((-critical - location[positive]) / scale[positive])
  pmin(1, pmax(0, probability))
}

loading_region_probabilities <- function(critical, common_noise, loading_lower, loading_upper) {
  # Relax the dependence of the conditional mean and scale on the loading.
  # At a fixed mean, symmetric-normal interval probability has no interior
  # minimum as a function of its nonnegative standard deviation.
  location <- abs(common_noise) * loading_upper
  pmin(normal_interval_probability(critical, location, sqrt(1 - loading_upper^2)),
       normal_interval_probability(critical, location, sqrt(1 - loading_lower^2)))
}

variance_region_vote_cdf <- function(critical, valid_minimum, votes_required,
                                    shared_bounds, private_bounds) {
  check_variance_region(shared_bounds, private_bounds)
  check_integer(valid_minimum, "valid_minimum", 1, nrow(private_bounds))
  check_integer(votes_required, "votes_required", 1, valid_minimum)
  if (length(critical) != 1L || !is.finite(critical) || critical < 0) {
    stop("critical must be finite and nonnegative")
  }
  if (critical == 0) return(0)
  loading_lower <- sqrt(shared_bounds[1L] / (shared_bounds[1L] + private_bounds[, 2L]))
  loading_upper <- if (shared_bounds[2L] == 0) rep(0, nrow(private_bounds)) else
    sqrt(shared_bounds[2L] / (shared_bounds[2L] + private_bounds[, 1L]))
  if (all(loading_lower == loading_upper) && all(loading_lower == loading_lower[1L])) {
    return(equicorrelated_vote_cdf(critical, valid_minimum, votes_required, loading_lower[1L]^2))
  }
  conditional_tail <- function(common_noise) {
    probability <- loading_region_probabilities(critical, common_noise, loading_lower, loading_upper)
    poisson_binomial_tail(sort(probability)[seq_len(valid_minimum)], votes_required)
  }
  integrand <- function(values) vapply(values, conditional_tail, numeric(1L)) * dnorm(values)
  # The endpoint minimum and changing worst subset create kinks. Integrating
  # separate finite intervals avoids roundoff failures in the infinite-domain
  # transformation. Omitting |U| > 8 only decreases the lower bound (by < 2e-15).
  knots <- sort(unique(c(0:8, if (any(loading_upper == 1)) min(8, critical))))
  integrate_piece <- function(lower, upper, absolute_tolerance = 1e-10, depth = 0L) {
    result <- integrate(integrand, lower, upper, rel.tol = 1e-8,
                        abs.tol = absolute_tolerance, subdivisions = 1000L, stop.on.error = FALSE)
    if (identical(result$message, "OK")) return(result$value)
    if (depth >= 6L) {
      stop(structure(list(message = paste("Variance-region integration failed:", result$message),
                          call = NULL), class = c("variance_region_integration_error", "error", "condition")))
    }
    # Retry smaller intervals with the SAME total absolute-error budget.
    midpoint <- (lower + upper) / 2
    integrate_piece(lower, midpoint, absolute_tolerance / 2, depth + 1L) +
      integrate_piece(midpoint, upper, absolute_tolerance / 2, depth + 1L)
  }
  pieces <- vapply(seq_len(length(knots) - 1L), function(index) {
    integrate_piece(knots[index], knots[index + 1L])
  }, numeric(1L))
  2 * sum(pieces)
}

make_source_region_calibration <- function(valid_minimum, votes_required,
                                           shared_bounds, private_bounds,
                                           anchor_variance_upper,
                                           alpha = .05, variance_alpha = .005,
                                           anchor_fraction = .5) {
  check_variance_region(shared_bounds, private_bounds)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(variance_alpha) != 1L || !is.finite(variance_alpha) ||
      variance_alpha < 0 || variance_alpha >= alpha ||
      length(anchor_variance_upper) != 1L || !is.finite(anchor_variance_upper) ||
      anchor_variance_upper <= 0) stop("Invalid variance uncertainty budget or anchor bound")
  check_integer(valid_minimum, "valid_minimum", 1, nrow(private_bounds))
  calibration <- make_source_calibration(nrow(private_bounds), valid_minimum, votes_required,
    alpha = alpha - variance_alpha, anchor_fraction = anchor_fraction)
  marginal_critical <- calibration$source_critical
  objective <- function(value) variance_region_vote_cdf(value, valid_minimum, votes_required,
    shared_bounds, private_bounds) - (1 - calibration$source_alpha)
  # If the envelope cannot improve on the universal marginal bound, retain it.
  # This also avoids an unnecessarily large root bracket for wide regions.
  candidate <- tryCatch({
    if (objective(marginal_critical) >= 0) {
      uniroot(objective, c(0, marginal_critical), tol = 1e-8)$root + 1e-7
    } else NULL
  }, variance_region_integration_error = function(error) error)
  envelope_critical <- NULL
  if (inherits(candidate, "variance_region_integration_error")) {
    # The marginal rule has its own coverage proof; failed quadrature is never
    # accepted as a numerical value. Record this distinct fallback explicitly.
    calibration$integration_fallback_reason <- conditionMessage(candidate)
  } else if (!is.null(candidate)) {
    envelope_critical <- candidate
    calibration$source_critical <- min(marginal_critical, envelope_critical)
  }
  calibration$method <- "shared_target_variance_region"
  calibration$calibration_bound <- if (!is.null(envelope_critical) &&
    envelope_critical < marginal_critical) "conditional_variance_region" else "marginal_bound"
  calibration$envelope_critical <- envelope_critical
  calibration$alpha <- alpha
  calibration$variance_alpha <- variance_alpha
  calibration$shared_bounds <- shared_bounds
  calibration$private_bounds <- private_bounds
  calibration$anchor_variance_upper <- anchor_variance_upper
  calibration$source_total_variances <- shared_bounds[2L] + private_bounds[, 2L]
  calibration
}

source_region_interval <- function(anchor_estimate, source_estimates, calibration) {
  if (!identical(calibration$method, "shared_target_variance_region")) {
    stop("A variance-region calibration is required")
  }
  source_search_interval(anchor_estimate, sqrt(calibration$anchor_variance_upper),
    source_estimates, sqrt(calibration$source_total_variances), calibration)
}

gaussian_variance_bounds <- function(sample_variances, degrees_freedom, failure_probability) {
  # Exact for variances of independent normal observations (or their means
  # when sample_variances have already been divided by the sample sizes).
  # Bonferroni requires no independence among the component variance estimates.
  if (!is.numeric(sample_variances) || !length(sample_variances) ||
      any(!is.finite(sample_variances)) || any(sample_variances < 0) ||
      !(length(degrees_freedom) %in% c(1L, length(sample_variances))) ||
      any(!is.finite(degrees_freedom)) || any(degrees_freedom <= 0) ||
      length(failure_probability) != 1L || !is.finite(failure_probability) ||
      failure_probability <= 0 || failure_probability >= 1) stop("Invalid Gaussian variance inputs")
  tail_probability <- failure_probability / (2 * length(sample_variances))
  bounds <- cbind(lower = degrees_freedom * sample_variances /
                   qchisq(tail_probability, degrees_freedom, lower.tail = FALSE),
                 upper = degrees_freedom * sample_variances /
                   qchisq(tail_probability, degrees_freedom))
  if (any(!is.finite(bounds))) stop("Variance bounds overflowed")
  bounds
}
