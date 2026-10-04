# Research reference constructions, separate from the frozen RoCE estimator.
# Coverage requires a valid target anchor and at least valid_minimum unbiased
# source candidates with the marginal calibration stated in theory.md.

check_integer <- function(value, name, lower, upper) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper) {
    stop(name, " must be an integer in [", lower, ", ", upper, "]", call. = FALSE)
  }
}

quorum_intervals <- function(lower, upper, min_votes) {
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != length(upper) ||
      !length(lower) || any(!is.finite(c(lower, upper))) || any(lower > upper)) {
    stop("Finite, ordered interval endpoints of equal positive length are required", call. = FALSE)
  }
  check_integer(min_votes, "min_votes", 1, length(lower))
  endpoints <- sort(unique(c(lower, upper)))
  starts <- tabulate(match(lower, endpoints), nbins = length(endpoints))
  ends <- tabulate(match(upper, endpoints), nbins = length(endpoints))
  count_left <- 0L
  open_interval <- FALSE
  result <- list()
  for (index in seq_along(endpoints)) {
    # Closed intervals include both endpoints, including isolated quorum points.
    count_at <- count_left + starts[index]
    count_right <- count_at - ends[index]
    if (!open_interval && count_at >= min_votes) {
      interval_start <- endpoints[index]
      open_interval <- TRUE
    }
    if (open_interval && count_right < min_votes) {
      result[[length(result) + 1L]] <- c(interval_start, endpoints[index])
      open_interval <- FALSE
    }
    count_left <- count_right
  }
  if (open_interval || count_left != 0L) stop("Interval sweep did not close")
  intervals <- if (length(result)) do.call(rbind, result) else matrix(numeric(), 0L, 2L)
  colnames(intervals) <- c("lower", "upper")
  intervals
}

equicorrelated_vote_cdf <- function(critical, valid_minimum, votes_required, correlation) {
  if (critical <= 0) return(0)
  if (correlation == 1 || valid_minimum == 1L) return(2 * pnorm(critical) - 1)
  if (correlation == 0) {
    return(pbinom(votes_required - 1L, valid_minimum, 2 * pnorm(critical) - 1,
                  lower.tail = FALSE))
  }
  integrand <- function(common_noise) {
    location <- sqrt(correlation) * common_noise
    scale <- sqrt(1 - correlation)
    probability <- pnorm((critical - location) / scale) - pnorm((-critical - location) / scale)
    probability <- pmin(1, pmax(0, probability))
    pbinom(votes_required - 1L, valid_minimum, probability, lower.tail = FALSE) * dnorm(common_noise)
  }
  integrate(integrand, -Inf, Inf, rel.tol = 1e-9, abs.tol = 1e-11,
            subdivisions = 1000L, stop.on.error = TRUE)$value
}

poisson_binomial_tail <- function(probabilities, votes_required) {
  if (!is.numeric(probabilities) || !length(probabilities) ||
      any(!is.finite(probabilities)) || any(probabilities < 0 | probabilities > 1)) {
    stop("Bernoulli probabilities must be a nonempty numeric vector in [0, 1]")
  }
  check_integer(votes_required, "votes_required", 0, length(probabilities))
  if (votes_required == 0L) return(1)
  # Keep masses below the threshold and an absorbing mass at/above it.
  below <- numeric(votes_required)
  below[1L] <- 1
  tail <- 0
  for (probability in probabilities) {
    previous <- below
    tail <- tail + previous[votes_required] * probability
    below <- previous * (1 - probability)
    if (votes_required > 1L) {
      below[2:votes_required] <- below[2:votes_required] +
        previous[1:(votes_required - 1L)] * probability
    }
  }
  min(1, max(0, tail))
}

shared_target_vote_cdf <- function(critical, valid_minimum, votes_required,
                                   shared_variance, source_private_variances) {
  if (length(shared_variance) != 1L || !is.finite(shared_variance) || shared_variance < 0 ||
      !is.numeric(source_private_variances) || !length(source_private_variances) ||
      any(!is.finite(source_private_variances)) || any(source_private_variances <= 0) ||
      any(!is.finite(shared_variance + source_private_variances))) {
    stop("A finite nonnegative shared variance and positive private variances are required")
  }
  check_integer(valid_minimum, "valid_minimum", 1, length(source_private_variances))
  check_integer(votes_required, "votes_required", 1, valid_minimum)
  if (critical <= 0) return(0)
  if (shared_variance == 0) return(equicorrelated_vote_cdf(critical, valid_minimum, votes_required, 0))
  if (all(source_private_variances == source_private_variances[1L])) {
    correlation <- shared_variance / (shared_variance + source_private_variances[1L])
    return(equicorrelated_vote_cdf(critical, valid_minimum, votes_required, correlation))
  }
  marginal_se <- sqrt(shared_variance + source_private_variances)
  private_se <- sqrt(source_private_variances)
  conditional_tail <- function(common_noise) {
    location <- sqrt(shared_variance) * common_noise
    probability <- pnorm((critical * marginal_se - location) / private_se) -
      pnorm((-critical * marginal_se - location) / private_se)
    # This adversarial subset gives a lower bound for any fixed valid subset.
    probability <- sort(pmin(1, pmax(0, probability)))[seq_len(valid_minimum)]
    poisson_binomial_tail(probability, votes_required)
  }
  integrand <- function(values) vapply(values, conditional_tail, numeric(1L)) * dnorm(values)
  2 * integrate(integrand, 0, Inf, rel.tol = 1e-8, abs.tol = 1e-10,
                subdivisions = 1000L, stop.on.error = TRUE)$value
}

make_source_calibration <- function(num_sources, valid_minimum, votes_required,
                                    alpha = .05, anchor_fraction = .5,
                                    method = c("marginal_bound", "equicorrelated_gaussian", "shared_target_gaussian"),
                                    source_correlation = NULL, shared_variance = NULL,
                                    source_private_variances = NULL) {
  check_integer(num_sources, "num_sources", 1, .Machine$integer.max)
  check_integer(valid_minimum, "valid_minimum", 0, num_sources)
  check_integer(votes_required, "votes_required", as.integer(valid_minimum > 0), valid_minimum)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(anchor_fraction) != 1L || !is.finite(anchor_fraction) ||
      anchor_fraction <= 0 || anchor_fraction >= 1) stop("Invalid probability allocation")
  method <- match.arg(method)
  if (valid_minimum == 0L) {
    if (!is.null(source_correlation) || !is.null(shared_variance) || !is.null(source_private_variances)) {
      stop("Source covariance is unused when no valid-source count is assumed")
    }
    return(list(num_sources = as.integer(num_sources), valid_minimum = 0L, votes_required = 0L,
                alpha = alpha, anchor_alpha = alpha, source_alpha = 0,
                method = "target_only", anchor_critical = qnorm(1 - alpha / 2)))
  }
  if (method != "shared_target_gaussian" &&
      (!is.null(shared_variance) || !is.null(source_private_variances))) {
    stop("Variance components are only used by shared_target_gaussian")
  }
  anchor_alpha <- alpha * anchor_fraction
  source_alpha <- alpha - anchor_alpha
  marginal_tail <- source_alpha * (valid_minimum - votes_required + 1) / valid_minimum
  calibration_bound <- method
  factor_critical <- source_total_variances <- NULL
  if (method == "marginal_bound") {
    if (!is.null(source_correlation)) stop("Correlation is not used by the marginal bound")
    source_critical <- qnorm(1 - marginal_tail / 2)
  } else if (method == "equicorrelated_gaussian") {
    if (length(source_correlation) != 1L || !is.finite(source_correlation) ||
        source_correlation < 0 || source_correlation > 1) {
      stop("Known standardized source-error correlation must be in [0, 1]")
    }
    if (source_correlation == 0) {
      probability <- qbeta(1 - source_alpha, votes_required, valid_minimum - votes_required + 1)
      source_critical <- qnorm((1 + probability) / 2)
    } else if (source_correlation == 1 || valid_minimum == 1L) {
      source_critical <- qnorm(1 - source_alpha / 2)
    } else {
      upper <- qnorm(1 - source_alpha / (2 * valid_minimum))
      objective <- function(value) equicorrelated_vote_cdf(
        value, valid_minimum, votes_required, source_correlation) - (1 - source_alpha)
      source_critical <- uniroot(objective, c(0, upper), tol = 1e-10)$root + 1e-8
    }
  } else {
    if (!is.null(source_correlation)) stop("Use variance components for shared_target_gaussian")
    if (length(source_private_variances) != num_sources) stop("Provide one private variance per source")
    # Validate components even before using the zero-CDF shortcut.
    shared_target_vote_cdf(0, valid_minimum, votes_required, shared_variance, source_private_variances)
    source_total_variances <- shared_variance + source_private_variances
    upper <- qnorm(1 - source_alpha / (2 * num_sources))
    objective <- function(value) shared_target_vote_cdf(value, valid_minimum, votes_required,
      shared_variance, source_private_variances) - (1 - source_alpha)
    factor_critical <- uniroot(objective, c(0, upper), tol = 1e-8)$root + 1e-7
    marginal_critical <- qnorm(1 - marginal_tail / 2)
    source_critical <- min(factor_critical, marginal_critical)
    calibration_bound <- if (factor_critical <= marginal_critical) "conditional_worst_subset" else "marginal_bound"
  }
  list(num_sources = as.integer(num_sources), valid_minimum = as.integer(valid_minimum),
       votes_required = as.integer(votes_required), alpha = alpha,
       anchor_alpha = anchor_alpha, source_alpha = source_alpha,
       marginal_tail = marginal_tail, method = method, source_correlation = source_correlation,
       calibration_bound = calibration_bound, factor_critical = factor_critical,
       shared_variance = shared_variance, source_private_variances = source_private_variances,
       source_total_variances = source_total_variances,
       anchor_critical = qnorm(1 - anchor_alpha / 2), source_critical = source_critical)
}

source_search_interval <- function(anchor_estimate, anchor_se, source_estimates,
                                   source_se, calibration) {
  if (length(anchor_estimate) != 1L || !is.finite(anchor_estimate) ||
      length(anchor_se) != 1L || !is.finite(anchor_se) || anchor_se <= 0) {
    stop("A finite anchor estimate and positive standard error are required")
  }
  anchor_interval <- anchor_estimate + c(-1, 1) * calibration$anchor_critical * anchor_se
  if (calibration$valid_minimum == 0L) {
    return(list(lower = anchor_interval[1L], upper = anchor_interval[2L],
                intervals = matrix(anchor_interval, 1L, dimnames = list(NULL, c("lower", "upper"))),
                used_anchor_fallback = FALSE, source_components = 0L, interval_components = 1L))
  }
  if (
      length(source_estimates) != calibration$num_sources ||
      length(source_se) != length(source_estimates) ||
      any(!is.finite(c(source_estimates, source_se))) || any(source_se <= 0)) {
    stop("Finite estimates and positive standard errors matching the calibration are required")
  }
  if (!is.null(calibration$source_total_variances) &&
      any(abs(source_se^2 - calibration$source_total_variances) > 1e-10 * calibration$source_total_variances)) {
    stop("Source standard errors do not match the calibrated variance components")
  }
  source_intervals <- quorum_intervals(
    source_estimates - calibration$source_critical * source_se,
    source_estimates + calibration$source_critical * source_se,
    calibration$votes_required)
  intervals <- source_intervals
  if (nrow(intervals)) {
    intervals[, "lower"] <- pmax(intervals[, "lower"], anchor_interval[1L])
    intervals[, "upper"] <- pmin(intervals[, "upper"], anchor_interval[2L])
    intervals <- intervals[intervals[, "lower"] <= intervals[, "upper"], , drop = FALSE]
  }
  fallback <- nrow(intervals) == 0L
  if (fallback) intervals <- matrix(anchor_interval, 1L, dimnames = list(NULL, c("lower", "upper")))
  # The hull is an ordinary reported CI; retain the actual confidence set too.
  list(lower = min(intervals[, "lower"]), upper = max(intervals[, "upper"]),
       intervals = intervals, used_anchor_fallback = fallback,
       source_components = nrow(source_intervals), interval_components = nrow(intervals))
}
