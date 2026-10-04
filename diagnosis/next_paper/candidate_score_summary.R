# Extract source-candidate moments from saved complete cross-fitted TATE fits.
# This is a diagnostic interface; it does not change fitting or aggregation.

score_mean_covariance <- function(values) {
  values <- as.matrix(values)
  if (!is.numeric(values) || nrow(values) < 2L || any(!is.finite(values))) {
    stop("Score values must be a finite numeric matrix with at least two rows")
  }
  centered <- sweep(values, 2L, colMeans(values), "-")
  crossprod(centered) / nrow(values)^2
}

check_score_indices <- function(indices, size, label) {
  if (!is.numeric(indices) || !length(indices) || any(!is.finite(indices)) ||
      any(indices != floor(indices)) || any(indices < 1L | indices > size) ||
      anyDuplicated(indices)) stop(label, ": invalid or duplicated observation indices")
  as.integer(indices)
}

pool_source_fold_moments <- function(fold_means, fold_variances, fold_sizes) {
  fold_means <- as.matrix(fold_means)
  fold_variances <- as.matrix(fold_variances)
  fold_sizes <- as.matrix(fold_sizes)
  if (!identical(dim(fold_means), dim(fold_variances)) ||
      !identical(dim(fold_means), dim(fold_sizes)) || !length(fold_means) ||
      any(!is.finite(c(fold_means, fold_variances, fold_sizes))) ||
      any(fold_variances < 0) || any(fold_sizes < 1 | fold_sizes != floor(fold_sizes))) {
    stop("Source fold moments require matching finite matrices and positive sample counts")
  }
  sizes <- colSums(fold_sizes)
  means <- colSums(fold_means * fold_sizes) / sizes
  # Fold variances use divisor n_fold, matching the saved score moments.
  # Include between-fold mean differences when centering at the site mean.
  differences <- sweep(fold_means, 2L, means, "-")
  variances <- colSums(fold_sizes * (fold_variances + differences^2)) / sizes^2
  list(means = means, variances = variances, sizes = sizes)
}

common_target_diagnostics <- function(predictions, comparison_sources) {
  if (anyDuplicated(comparison_sources) || any(!comparison_sources %in% colnames(predictions))) {
    stop("Comparison sources must be distinct prediction columns")
  }
  count <- length(comparison_sources)
  result <- list(comparison_sources = comparison_sources, comparison_count = count,
    common_variance = NA_real_, relative_covariance_error = NA_real_,
    leading_eigenvalue_fraction = NA_real_, minimum_correlation = NA_real_,
    prediction_disagreement_rmse = NA_real_)
  if (count < 2L) return(result) # A single valid source cannot test a common-factor restriction.
  values <- predictions[, comparison_sources, drop = FALSE]
  covariance <- score_mean_covariance(values)
  common <- rowMeans(values)
  common_variance <- as.numeric(score_mean_covariance(common))
  denominator <- sqrt(sum(covariance^2))
  eigenvalues <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  marginal_variances <- diag(covariance)
  result$common_variance <- common_variance
  result$relative_covariance_error <- if (denominator > 0)
    sqrt(sum((covariance - common_variance)^2)) / denominator else 0
  result$leading_eigenvalue_fraction <- if (sum(eigenvalues) > 0)
    max(eigenvalues) / sum(eigenvalues) else NA_real_
  if (all(marginal_variances > 0)) {
    correlation <- covariance / sqrt(outer(marginal_variances, marginal_variances))
    result$minimum_correlation <- min(correlation[upper.tri(correlation)])
  }
  result$prediction_disagreement_rmse <- sqrt(mean((values - common)^2))
  result
}

extract_candidate_score_summary <- function(fitted_tate, data_split,
                                             comparison_sources = NULL, keep_scores = FALSE) {
  sites <- setdiff(names(data_split), "t")
  if (!length(sites) || !identical(names(fitted_tate$source_estimates), sites)) {
    stop("Source order in the fitted estimates must match the data sites")
  }
  if (length(keep_scores) != 1L || is.na(keep_scores) || !is.logical(keep_scores)) {
    stop("keep_scores must be TRUE or FALSE")
  }
  if (is.null(comparison_sources)) comparison_sources <- sites
  sizes <- vapply(data_split[c("t", sites)], function(site) as.integer(site$n), integer(1L))
  if (any(sizes < 2L)) stop("Each site needs at least two observations")
  target <- target_centered_by_fold <- matrix(NA_real_, sizes[1L], length(sites) + 1L,
    dimnames = list(NULL, c("target_anchor", sites)))
  residuals <- residuals_centered_by_fold <- lapply(sizes[-1L], function(size) rep(NA_real_, size))
  target_assignments <- integer(sizes[1L])
  source_assignments <- lapply(sizes[-1L], integer)
  folds <- fitted_tate$intermediates$fold_info
  if (length(folds) != fitted_tate$n_folds) stop("Saved outer-fold count is inconsistent")
  fold_means <- fold_variances <- fold_sizes <- matrix(NA_real_, length(folds), length(sites),
    dimnames = list(NULL, sites))
  fold_number <- 0L
  for (fold in folds) {
    fold_number <- fold_number + 1L
    target_indices <- check_score_indices(fold$target_idx, sizes[1L], "Target fold")
    if (length(fold$varphi_ot) != length(target_indices) ||
        any(vapply(fold$zeta_components, length, integer(1L)) != length(target_indices)) ||
        length(fold$zeta_components) != length(sites) ||
        length(fold$xi_components) != length(sites)) stop("Fold score dimensions are inconsistent")
    target[target_indices, 1L] <- fold$varphi_ot + fold$fold_target_estimate
    target_centered_by_fold[target_indices, 1L] <- fold$varphi_ot
    target_assignments <- target_assignments + tabulate(target_indices, nbins = sizes[1L])
    for (index in seq_along(sites)) {
      source_indices <- check_score_indices(fold$source_idx[[index]], sizes[index + 1L], "Source fold")
      if (length(fold$xi_components[[index]]) != length(source_indices)) {
        stop("Source score and observation-index lengths differ")
      }
      target[target_indices, index + 1L] <- fold$zeta_components[[index]] + fold$mu_pred_ts[index]
      target_centered_by_fold[target_indices, index + 1L] <- fold$zeta_components[[index]]
      residuals[[index]][source_indices] <- fold$xi_components[[index]] + fold$delta_ts[index]
      residuals_centered_by_fold[[index]][source_indices] <- fold$xi_components[[index]]
      source_assignments[[index]] <- source_assignments[[index]] +
        tabulate(source_indices, nbins = sizes[index + 1L])
      fold_means[fold_number, index] <- fold$delta_ts[index]
      fold_variances[fold_number, index] <- fold$V_s[index]
      fold_sizes[fold_number, index] <- length(source_indices)
    }
  }
  if (any(target_assignments != 1L) || any(unlist(source_assignments) != 1L)) {
    stop("Outer folds must cover every site observation exactly once")
  }
  target_covariance <- score_mean_covariance(target)
  private_variances <- vapply(residuals, function(values) as.numeric(score_mean_covariance(values)), numeric(1L))
  pooled <- pool_source_fold_moments(fold_means, fold_variances, fold_sizes)
  message_error <- max(abs(pooled$means - vapply(residuals, mean, numeric(1L))),
                       abs(pooled$variances - private_variances))
  if (message_error > 1e-12) stop("Source summary messages disagree with individual score components")
  covariance <- target_covariance
  diag(covariance)[-1L] <- diag(covariance)[-1L] + private_variances
  estimates <- colMeans(target) + c(0, vapply(residuals, mean, numeric(1L)))
  reference <- c(target_anchor = fitted_tate$target_only$estimate, fitted_tate$source_estimates)
  estimate_error <- max(abs(estimates - reference))
  anchor_variance_error <- abs(covariance[1L, 1L] - fitted_tate$target_only$variance)
  if (!is.finite(estimate_error) || estimate_error > 1e-12 ||
      !is.finite(anchor_variance_error) || anchor_variance_error > 1e-12) {
    stop("Reconstructed candidate means or target variance disagree with the saved fit")
  }
  result <- list(estimates = estimates, covariance = covariance, site_sizes = sizes,
    target_covariance = target_covariance, source_private_variances = private_variances,
    target_covariance_within_folds = score_mean_covariance(target_centered_by_fold),
    private_variances_within_folds = vapply(residuals_centered_by_fold,
      function(values) as.numeric(score_mean_covariance(values)), numeric(1L)),
    common_target = common_target_diagnostics(target[, sites, drop = FALSE], comparison_sources),
    common_target_within_folds = common_target_diagnostics(
      target_centered_by_fold[, sites, drop = FALSE], comparison_sources),
    estimate_identity_error = estimate_error, anchor_variance_identity_error = anchor_variance_error,
    source_message_identity_error = message_error,
    source_message_moments = list(means = fold_means, variances = fold_variances, sizes = fold_sizes),
    inference_scope = "Empirical candidate-score moments; nuisance remainder and Gaussian inference require separate validation")
  if (keep_scores) result$score_records <- list(target = target, source = residuals)
  result
}
