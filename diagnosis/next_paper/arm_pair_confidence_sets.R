# Dependence-robust Gaussian-marginal reference for different valid arm sets.
# Source source_confidence_sets.R first for integer checks and endpoint sweep.

arm_pair_moments <- function(means, covariance, tolerance = 1e-10) {
  means <- as.matrix(means)
  covariance <- as.matrix(covariance)
  if (!is.numeric(means) || ncol(means) != 2L || nrow(means) < 1L || any(!is.finite(means)) ||
      !identical(colnames(means), c("mu1", "mu0")) || is.null(rownames(means)) ||
      anyDuplicated(rownames(means)) || any(!nzchar(rownames(means))) ||
      rownames(means)[1L] != "target_anchor") {
    stop("Means require mu1/mu0 columns and distinct candidate row names, with the target first")
  }
  width <- nrow(means)
  labels <- c(paste0("mu1:", rownames(means)), paste0("mu0:", rownames(means)))
  if (!is.numeric(covariance) || !identical(dim(covariance), c(2L * width, 2L * width)) ||
      any(!is.finite(covariance)) || !identical(rownames(covariance), labels) ||
      !identical(colnames(covariance), labels)) stop("Covariance labels must follow all mu1 candidates then all mu0 candidates")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) stop("tolerance must be positive")
  scale <- max(abs(covariance))
  if (max(abs(covariance - t(covariance))) > tolerance * scale) stop("Covariance must be symmetric")
  covariance <- (covariance + t(covariance)) / 2
  if (min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values) < -tolerance * scale) {
    stop("Covariance is not positive semidefinite")
  }
  first <- seq_len(width)
  second <- width + first
  variances <- outer(diag(covariance)[first], diag(covariance)[second], "+") -
    2 * covariance[first, second, drop = FALSE]
  if (min(variances) < -4 * tolerance * scale) stop("Negative contrast variance exceeds rounding tolerance")
  pairs <- expand.grid(treated = rownames(means), control = rownames(means), stringsAsFactors = FALSE)
  pairs$estimate <- as.vector(outer(means[, "mu1"], means[, "mu0"], "-"))
  pairs$variance <- as.vector(pmax(variances, 0))
  pairs$se <- sqrt(pairs$variance)
  pairs
}

arm_pair_search_interval <- function(means, covariance, valid_mu1, valid_mu0,
                                     votes_required = NULL, alpha = .05, anchor_fraction = .5) {
  pairs <- arm_pair_moments(means, covariance)
  count <- nrow(means) - 1L
  check_integer(valid_mu1, "valid_mu1", 0, count)
  check_integer(valid_mu0, "valid_mu0", 0, count)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      length(anchor_fraction) != 1L || !is.finite(anchor_fraction) ||
      anchor_fraction <= 0 || anchor_fraction >= 1) stop("alpha and anchor_fraction must be in (0,1)")
  valid_pairs <- (valid_mu1 + 1) * (valid_mu0 + 1)
  anchor_only <- valid_mu1 == 0L && valid_mu0 == 0L
  anchor_alpha <- if (anchor_only) alpha else alpha * anchor_fraction
  anchor <- pairs$estimate[1L] + c(-1, 1) * qnorm(anchor_alpha / 2, lower.tail = FALSE) * pairs$se[1L]
  if (is.null(votes_required)) votes_required <- max(1, floor(valid_pairs / 2))
  check_integer(votes_required, "votes_required", 1, valid_pairs)
  pair_tail <- if (anchor_only) NA_real_ else (alpha - anchor_alpha) * (valid_pairs - votes_required + 1) / valid_pairs
  intervals <- matrix(anchor, 1L, dimnames = list(NULL, c("lower", "upper")))
  fallback <- FALSE
  source_components <- 0L
  if (!anchor_only) {
    critical <- qnorm(pair_tail / 2, lower.tail = FALSE)
    intervals <- quorum_intervals(pairs$estimate - critical * pairs$se,
                                  pairs$estimate + critical * pairs$se, votes_required)
    source_components <- nrow(intervals)
    if (nrow(intervals)) {
      intervals[, "lower"] <- pmax(intervals[, "lower"], anchor[1L])
      intervals[, "upper"] <- pmin(intervals[, "upper"], anchor[2L])
      intervals <- intervals[intervals[, "lower"] <= intervals[, "upper"], , drop = FALSE]
    }
    fallback <- nrow(intervals) == 0L
    if (fallback) intervals <- matrix(anchor, 1L, dimnames = list(NULL, c("lower", "upper")))
  }
  list(lower = min(intervals[, "lower"]), upper = max(intervals[, "upper"]),
       intervals = intervals, used_anchor_fallback = fallback, source_components = source_components,
       anchor_only = anchor_only, valid_pairs = valid_pairs, votes_required = votes_required,
       pair_tail = pair_tail, anchor_alpha = anchor_alpha, pairs = pairs)
}
