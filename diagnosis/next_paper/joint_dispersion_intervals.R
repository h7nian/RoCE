# Joint Gaussian bias-dispersion reference. Source the checked scalar helpers
# in dispersion_bias_intervals.R first. This is not a fitted-RoCE validity claim.

joint_mean_gls <- function(covariance, design, indices = seq_len(nrow(covariance))) {
  precision <- chol2inv(chol(covariance[indices, indices, drop = FALSE]))
  selected_design <- design[indices, , drop = FALSE]
  precision_design <- precision %*% selected_design
  mean_covariance <- chol2inv(chol(crossprod(selected_design, precision_design)))
  mean_weights <- mean_covariance %*% t(precision_design)
  coefficients <- numeric(nrow(covariance))
  coefficients[indices] <- drop(c(1, -1) %*% mean_weights)
  list(coefficients = coefficients, variance = drop(crossprod(c(1, -1), mean_covariance %*% c(1, -1))),
       mean_weights = mean_weights, mean_covariance = mean_covariance, precision = precision)
}

shared_prediction_variances <- function(covariance, tolerance = 1e-12,
                                       allow_private_correlation = FALSE) {
  width <- nrow(covariance) / 2L
  count <- width - 1L
  if (count < 2L) stop("Shared-prediction verification needs at least two sources")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) stop("A positive structure tolerance is required")
  if(!is.logical(allow_private_correlation) || length(allow_private_correlation)!=1L || is.na(allow_private_correlation)) stop("allow_private_correlation must be TRUE or FALSE")
  sources <- list(2:width, width + (2:width))
  anchors <- c(1L, width + 1L)
  prediction <- matrix(0, 2L, 2L)
  cross <- matrix(0, 2L, 2L)
  private <- matrix(0, count, 2L)
  error <- 0
  for (arm in 1:2) {
    block <- covariance[sources[[arm]], sources[[arm]], drop = FALSE]
    shared <- block[1L, 2L]
    error <- max(error, abs(block[row(block) != col(block)] - shared))
    prediction[arm, arm] <- shared
    private[, arm] <- diag(block) - shared
    for (target_arm in 1:2) {
      values <- covariance[anchors[target_arm], sources[[arm]]]
      cross[target_arm, arm] <- values[1L]
      error <- max(error, abs(values - values[1L]))
    }
  }
  block <- covariance[sources[[1L]], sources[[2L]], drop = FALSE]
  shared_cross <- if(allow_private_correlation) block[1L,2L] else block[1L,1L]
  prediction[1L, 2L] <- prediction[2L, 1L] <- shared_cross
  checked_cross <- if(allow_private_correlation) block[row(block)!=col(block)] else block
  error <- max(error, abs(checked_cross-shared_cross))
  private_cross <- diag(block)-shared_cross
  scale <- max(abs(covariance))
  common <- rbind(cbind(covariance[anchors, anchors], cross), cbind(t(cross), prediction))
  if (error > tolerance * scale || min(private) < -tolerance * scale ||
      min(eigen(common, symmetric = TRUE, only.values = TRUE)$values) < -tolerance * scale) {
    stop(if(allow_private_correlation)
      "Covariance does not have shared predictions and independent source-private blocks" else
      "Covariance does not have shared predictions and zero cross-arm private covariance; no structural approximation is made")
  }
  result <- list(private_variances = pmax(private, 0), maximum_structure_error = error)
  if(allow_private_correlation) {
    deviation_scale <- sqrt(result$private_variances[,1L])*sqrt(result$private_variances[,2L])
    if(any(abs(private_cross)>deviation_scale+tolerance*scale) ||
       any(deviation_scale==0 & private_cross!=0)) stop("Invalid private cross-arm covariance block")
    result$private_cross_covariances <- private_cross
  }
  result
}

inflate_shared_private_covariance <- function(covariance) {
  structure <- shared_prediction_variances(covariance,allow_private_correlation=TRUE)
  private <- structure$private_variances
  denominator <- sqrt(private[,1L])*sqrt(private[,2L])
  correlation <- numeric(nrow(private))
  positive <- denominator>0
  correlation[positive] <- structure$private_cross_covariances[positive]/denominator[positive]
  inflation <- 1+abs(correlation)
  inflated <- covariance
  width <- nrow(private)+1L
  for(site in seq_len(nrow(private))) {
    positions <- c(site+1L,width+site+1L)
    diagonal <- cbind(positions,positions)
    inflated[diagonal] <- covariance[diagonal]+(inflation[site]-1)*private[site,]
    inflated[positions[1L],positions[2L]] <- inflated[positions[2L],positions[1L]] <-
      covariance[positions[1L],positions[2L]]-structure$private_cross_covariances[site]
  }
  list(covariance=inflated,inflation=inflation,maximum_structure_error=structure$maximum_structure_error)
}

validate_joint_score_bias <- function(value,width) {
  if(is.null(value)) return(NULL)
  if(!is.numeric(value) || !length(value) || any(!is.finite(value)) || any(value<0)) {
    stop("valid_score_bias must contain finite nonnegative mean-error bounds")
  }
  if(is.matrix(value) && (any(dim(value)!=c(width,2L)) ||
      (!is.null(colnames(value)) && !identical(colnames(value),c("mu1","mu0"))))) {
    stop("A valid_score_bias matrix must have target/source rows and mu1/mu0 columns")
  }
  if(length(value)==1L) value <- rep(value,2L*width)
  if(length(value)!=2L*width) stop("valid_score_bias must be a scalar or match both arm blocks")
  if(all(value==0)) return(NULL)
  as.numeric(value)
}

shared_prediction_bias_envelope <- function(covariance,bound_covariance,design,known,reference,counts,budgets) {
  width <- nrow(covariance)/2L
  private <- shared_prediction_variances(bound_covariance)$private_variances
  active <- which(counts>0L)
  if(length(active) && any(private[,active,drop=FALSE]<=0)) {
    stop("Shared mean-error bounds require positive private variances; use known_valid_upper for a reference-only bound")
  }
  precision_design <- reference$precision%*%design[known,,drop=FALSE]
  residual <- matrix(0,2L*width,2L*width)
  residual[known,known] <- reference$precision-
    precision_design%*%reference$mean_covariance%*%t(precision_design)
  projection <- diag(2L*width)-residual%*%covariance
  transformed <- drop(crossprod(abs(projection),budgets))
  anchors <- c(1L,width+1L)
  anchor <- bound_covariance[anchors,anchors]
  cross <- bound_covariance[anchors,c(2L,width+2L)]
  prediction <- matrix(c(bound_covariance[2L,3L],bound_covariance[2L,width+3L],
    bound_covariance[2L,width+3L],bound_covariance[width+2L,width+3L]),2L)
  curvature <- anchor+prediction-cross-t(cross)
  contrast <- c(1,-1)
  right <- drop((anchor-t(cross))%*%contrast)
  source_bounds <- c(max(transformed[2:width]),max(transformed[width+(2:width)]))
  limits <- lapply(1:2,function(arm) {
    if(counts[arm]==0L) return(0)
    ordered <- sort(private[,arm])
    unique(c(1/sum(1/head(ordered,counts[arm])),1/sum(1/tail(ordered,counts[arm]))))
  })
  corners <- expand.grid(limits)
  values <- apply(corners,1L,function(pool) {
    beta <- numeric(2L)
    if(length(active)) {
      matrix <- curvature[active,active,drop=FALSE]+diag(pool[active],nrow=length(active))
      beta[active] <- chol2inv(chol(matrix))%*%right[active]
    }
    sum(abs(contrast-beta)*transformed[anchors])+sum(abs(beta)*source_bounds)
  })
  if(any(!is.finite(values))) stop("Nonfinite shared mean-error envelope")
  max(values)
}

joint_score_bias_allowance <- function(calibration,noncentrality_upper) {
  if(is.null(calibration$valid_score_bias)) {
    return(sqrt(calibration$bias_factor*noncentrality_upper))
  }
  reference_factor <- max(0,calibration$reference_variance-calibration$variance)
  reference <- calibration$reference_bias_allowance+sqrt(reference_factor*noncentrality_upper)
  pairs <- calibration$valid_subset_bias_bounds
  if(is.null(pairs)) return(reference)
  # Each row is a paired (mean-error envelope, variance difference). Unknown
  # validity patterns require a maximum across subset rows before the minimum
  # with the separate known-reference bound.
  subset <- rep(0,length(noncentrality_upper))
  for(index in seq_len(nrow(pairs))) {
    subset <- pmax(subset,pairs[index,"bias"]+sqrt(pairs[index,"variance"]*noncentrality_upper))
  }
  pmin(reference,subset)
}

make_shared_prediction_arm_calibration <- function(covariance, valid_mu1, valid_mu0) {
  # Exact marginal comparator for large K. Its point/interval program remains
  # the earlier armwise program; only combinatorial subset search is replaced.
  calibration <- make_arm_dispersion_calibration(covariance, valid_mu1, valid_mu0, "uniform_upper")
  structure <- shared_prediction_variances(covariance)
  validity <- c(valid_mu1, valid_mu0)
  for (arm in 1:2) {
    name <- c("mu1", "mu0")[arm]
    selected <- c(1L, 1L + head(order(structure$private_variances[, arm], decreasing = TRUE), validity[arm]))
    block <- calibration[[name]]$covariance[selected, selected, drop = FALSE]
    worst_variance <- 1 / sum(chol2inv(chol(block)))
    calibration[[name]]$bias_factor <- if (validity[arm] == calibration$width - 1L) 0 else
      max(0, worst_variance - calibration[[name]]$variance)
    calibration[[name]]$subsets_examined <- 1L
    calibration[[name]]$bound_method <- "shared_prediction_exact"
  }
  calibration$bound_method <- "shared_prediction_exact"
  calibration
}

make_joint_dispersion_calibration <- function(covariance, valid_mu1, valid_mu0,
    bound_method = c("subset_exact", "shared_prediction_exact", "shared_prediction_bound", "known_valid_upper"),
    max_subsets = 50000L,valid_score_bias=NULL) {
  checked <- checked_dispersion_covariance(covariance)
  covariance <- checked$covariance
  width <- nrow(covariance) / 2L
  if (width != floor(width) || width < 2L) stop("The covariance must contain a target and sources in both arms")
  count <- width - 1L
  score_bias <- validate_joint_score_bias(valid_score_bias,width)
  validity <- c(valid_mu1, valid_mu0)
  if (length(valid_mu1) != 1L || length(valid_mu0) != 1L || any(!is.finite(validity)) ||
      any(validity < 0 | validity > count | validity != floor(validity))) stop("Invalid arm-specific guaranteed source counts")
  bound_method <- match.arg(bound_method)
  design <- cbind(mu1 = rep(c(1, 0), each = width), mu0 = rep(c(0, 1), each = width))
  full <- joint_mean_gls(covariance, design)
  known <- sort(unique(c(1L, width + 1L, if (valid_mu1 == count) seq_len(width),
                        if (valid_mu0 == count) width + seq_len(width))))
  reference <- joint_mean_gls(covariance, design, known)
  uncertain <- setdiff(seq_len(2L * width), known)
  # Estimate the unrestricted biases of uncertain coordinates after using the
  # known-valid observations to estimate both means and predict correlated noise.
  bias_map <- matrix(0, length(uncertain), 2L * width)
  bias_precision <- matrix(numeric(), 0L, 0L)
  if (length(uncertain)) {
    noise_prediction <- covariance[uncertain, known, drop = FALSE] %*% reference$precision
    mean_design <- design[uncertain, , drop = FALSE] - noise_prediction %*% design[known, , drop = FALSE]
    bias_map[, uncertain] <- diag(length(uncertain))
    bias_map[, known] <- -noise_prediction - mean_design %*% reference$mean_weights
    bias_covariance <- bias_map %*% covariance %*% t(bias_map)
    bias_precision <- chol2inv(chol(bias_covariance))
  }
  examined <- 0L
  structure_error <- NA_real_
  inflation <- NULL
  subset_bias_bounds <- NULL
  uncapped_subset_variance <- NULL
  worst <- known
  if (!length(uncertain)) {
    worst_variance <- full$variance
  } else if (bound_method == "known_valid_upper") {
    worst_variance <- reference$variance
  } else if (bound_method %in% c("shared_prediction_exact","shared_prediction_bound")) {
    bound_covariance <- covariance
    if(bound_method=="shared_prediction_bound") {
      inflated <- inflate_shared_private_covariance(covariance)
      bound_covariance <- inflated$covariance
      inflation <- inflated$inflation
    }
    structure <- shared_prediction_variances(bound_covariance)
    structure_error <- structure$maximum_structure_error
    if(!is.null(inflation)) structure_error <- max(structure_error,inflated$maximum_structure_error)
    selected <- lapply(1:2, function(arm) head(order(structure$private_variances[, arm], decreasing = TRUE), validity[arm]))
    worst <- sort(c(1L, width + 1L, 1L + selected[[1L]], width + 1L + selected[[2L]]))
    worst_variance <- joint_mean_gls(bound_covariance, design, worst)$variance
    uncapped_subset_variance <- worst_variance
    if(bound_method=="shared_prediction_bound") {
      # Every permitted subset contains the known-valid reference. Its
      # original variance is another valid upper bound, without inflating
      # known-valid observations that the reference already uses.
      worst_variance <- min(worst_variance,reference$variance)
    }
    examined <- 1L
  } else {
    if (length(max_subsets) != 1L || !is.finite(max_subsets) || max_subsets < 1 || max_subsets != floor(max_subsets)) {
      stop("max_subsets must be a positive integer")
    }
    if (sum(lchoose(count, validity)) > log(max_subsets) + 1e-12) {
      stop("Joint subset enumeration exceeds max_subsets; explicitly select a justified structural or upper bound")
    }
    subsets <- lapply(validity, function(number)
      if (number == 0L) list(integer()) else combn(count, number, simplify = FALSE))
    if(!is.null(score_bias)) {
      subset_bias_bounds <- matrix(0,nrow=prod(lengths(subsets)),ncol=2L,
        dimnames=list(NULL,c("bias","variance")))
    }
    worst_variance <- -Inf
    for (first in subsets[[1L]]) for (second in subsets[[2L]]) {
      indices <- c(1L, width + 1L, first + 1L, second + width + 1L)
      candidate_fit <- joint_mean_gls(covariance, design, indices)
      candidate <- candidate_fit$variance
      if(!is.null(score_bias)) {
        subset_bias_bounds[examined+1L,] <-
          c(sum(abs(candidate_fit$coefficients)*score_bias),max(0,candidate-full$variance))
      }
      if (candidate > worst_variance) { worst_variance <- candidate; worst <- sort(indices) }
      examined <- examined + 1L
    }
  }
  factor <- worst_variance - full$variance
  if (!is.finite(factor) || factor < -1e-10 * max(diag(covariance))) stop("Invalid joint bias factor")
  result <- list(covariance = covariance, width = width, design = design, coefficients = full$coefficients,
       variance = full$variance, reference_coefficients = reference$coefficients,
       reference_variance = reference$variance, known_valid = known, uncertain = uncertain,
       bias_map = bias_map, bias_precision = bias_precision, degrees_freedom = length(uncertain),
       bias_factor = max(factor, 0), bound_subset = worst, subsets_examined = examined,
       bound_method = bound_method, maximum_structure_error = structure_error,
       anchor_only = all(validity == 0))
  if(bound_method=="shared_prediction_bound") {
    result$private_variance_inflation <- inflation
    result$subset_variance_upper <- worst_variance
  }
  if(!is.null(score_bias)) {
    if(length(uncertain) && !result$anchor_only &&
       bound_method%in%c("shared_prediction_exact","shared_prediction_bound")) {
      envelope <- shared_prediction_bias_envelope(covariance,bound_covariance,design,known,reference,validity,score_bias)
      subset_bias_bounds <- matrix(c(envelope,max(0,uncapped_subset_variance-full$variance)),
        nrow=1L,dimnames=list(NULL,c("bias","variance")))
    }
    result$valid_score_bias <- score_bias
    result$reference_bias_allowance <- sum(abs(reference$coefficients)*score_bias)
    result$valid_subset_bias_bounds <- subset_bias_bounds
  }
  result
}

intersect_reference_interval <- function(lower, upper, reference_lower, reference_upper,
                                         empty_intersection = c("shorter", "reference")) {
  policy <- match.arg(empty_intersection)
  inputs <- list(lower, upper, reference_lower, reference_upper)
  if (!length(lower) || any(vapply(inputs,length,integer(1L)) != length(lower)) ||
      any(!is.finite(unlist(inputs))) || any(lower>upper | reference_lower>reference_upper)) {
    stop("Finite, nonempty component intervals of matching lengths are required")
  }
  combined_lower <- pmax(lower, reference_lower)
  combined_upper <- pmin(upper, reference_upper)
  empty <- combined_lower > combined_upper
  use_reference <- empty & (policy == "reference" | reference_upper-reference_lower <= upper-lower)
  use_original <- empty & !use_reference
  combined_lower[use_reference] <- reference_lower[use_reference]
  combined_upper[use_reference] <- reference_upper[use_reference]
  combined_lower[use_original] <- lower[use_original]
  combined_upper[use_original] <- upper[use_original]
  list(lower=combined_lower, upper=combined_upper, empty_intersection=empty, reference_fallback=use_reference)
}

joint_dispersion_intervals <- function(estimates, calibration, alpha = .05,
    dispersion_fraction = .5, reference_fraction = .5, empty_intersection = c("shorter", "reference")) {
  empty_intersection <- match.arg(empty_intersection)
  estimates <- as.matrix(estimates)
  if (!is.numeric(estimates) || ncol(estimates) != length(calibration$coefficients) || any(!is.finite(estimates))) {
    stop("Finite draws must match the joint calibration dimension")
  }
  probabilities <- c(alpha, dispersion_fraction, reference_fraction)
  if (length(alpha) != 1L || length(dispersion_fraction) != 1L ||
      (!is.null(reference_fraction) && length(reference_fraction) != 1L) ||
      any(!is.finite(probabilities)) || any(probabilities <= 0 | probabilities >= 1)) stop("Error allocations must lie in (0,1)")
  reference <- drop(estimates %*% calibration$reference_coefficients)
  reference_bias <- if(is.null(calibration$reference_bias_allowance)) 0 else calibration$reference_bias_allowance
  if (calibration$anchor_only || calibration$degrees_freedom == 0L) {
    radius <- qnorm(alpha / 2, lower.tail = FALSE) * sqrt(calibration$reference_variance)+reference_bias
    return(data.frame(estimate = reference, lower = reference - radius, upper = reference + radius,
      statistic = 0, noncentrality_upper = 0, bias_allowance = reference_bias, reference_fallback = FALSE, empty_intersection = FALSE))
  }
  interval_alpha <- alpha * if (is.null(reference_fraction)) 1 else 1 - reference_fraction
  dispersion_alpha <- interval_alpha * dispersion_fraction
  noise_alpha <- interval_alpha - dispersion_alpha
  contrasts <- estimates %*% t(calibration$bias_map)
  statistic <- rowSums((contrasts %*% calibration$bias_precision) * contrasts)
  if (any(statistic < -1e-10)) stop("Negative joint dispersion statistic")
  upper_noncentrality <- noncentrality_upper_bound(pmax(statistic, 0), calibration$degrees_freedom, dispersion_alpha)
  bias <- joint_score_bias_allowance(calibration,upper_noncentrality)
  point <- drop(estimates %*% calibration$coefficients)
  radius <- qnorm(noise_alpha / 2, lower.tail = FALSE) * sqrt(calibration$variance) + bias
  lower <- point - radius; upper <- point + radius
  fallback <- rep(FALSE, nrow(estimates))
  empty <- fallback
  if (!is.null(reference_fraction)) {
    reference_radius <- qnorm(alpha * reference_fraction / 2, lower.tail = FALSE) * sqrt(calibration$reference_variance)+reference_bias
    reference_lower <- reference - reference_radius; reference_upper <- reference + reference_radius
    combined <- intersect_reference_interval(lower,upper,reference_lower,reference_upper,empty_intersection)
    lower <- combined$lower; upper <- combined$upper
    fallback <- combined$reference_fallback; empty <- combined$empty_intersection
  }
  data.frame(estimate = point, lower = lower, upper = upper, statistic = statistic,
             noncentrality_upper = upper_noncentrality, bias_allowance = bias,
             reference_fallback = fallback, empty_intersection = empty)
}
