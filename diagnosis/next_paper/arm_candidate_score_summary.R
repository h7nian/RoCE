# Joint arm-score moments for the next-paper candidate interface.
# Source candidate_score_summary.R first. No nuisance models are refitted.

pool_source_fold_cross_covariance <- function(means_mu1, means_mu0, covariances, sizes) {
  matrices <- lapply(list(means_mu1, means_mu0, covariances, sizes), as.matrix)
  if (!length(matrices[[1L]]) ||
      any(!vapply(matrices, function(value) identical(dim(value), dim(matrices[[1L]])), logical(1L))) ||
      any(!is.finite(unlist(matrices))) || any(matrices[[4L]] < 1 | matrices[[4L]] != floor(matrices[[4L]]))) {
    stop("Cross-arm messages require matching finite matrices and positive integer counts")
  }
  means_mu1 <- matrices[[1L]]
  means_mu0 <- matrices[[2L]]
  covariances <- matrices[[3L]]
  sizes <- matrices[[4L]]
  total_sizes <- colSums(sizes)
  centered_mu1 <- sweep(means_mu1, 2L, colSums(sizes * means_mu1) / total_sizes, "-")
  centered_mu0 <- sweep(means_mu0, 2L, colSums(sizes * means_mu0) / total_sizes, "-")
  colSums(sizes * (covariances + centered_mu1 * centered_mu0)) / total_sizes^2
}

extract_joint_arm_score_summary <- function(fitted_tate, data_split, keep_scores = FALSE) {
  if (!is.logical(keep_scores) || length(keep_scores) != 1L || is.na(keep_scores)) {
    stop("keep_scores must be TRUE or FALSE")
  }
  if (!all(c("mu1", "mu0") %in% names(fitted_tate$arm_results))) {
    stop("Saved TATE must contain both arm fits")
  }
  arms <- lapply(fitted_tate$arm_results[c("mu1", "mu0")], function(fit)
    extract_candidate_score_summary(fit, data_split, keep_scores = TRUE))
  candidates <- names(arms$mu1$estimates)
  if (!identical(candidates, names(arms$mu0$estimates))) stop("Arm candidate order differs")
  width <- length(candidates)
  sources <- candidates[-1L]
  target_scores <- cbind(arms$mu1$score_records$target, arms$mu0$score_records$target)
  labels <- c(paste0("mu1:", candidates), paste0("mu0:", candidates))
  colnames(target_scores) <- labels
  target_covariance <- score_mean_covariance(target_scores)
  covariance <- target_covariance
  source_covariances <- vector("list", length(sources))
  names(source_covariances) <- sources
  fold_covariances <- matrix(NA_real_, fitted_tate$n_folds, length(sources), dimnames = list(NULL, sources))
  for (fold_index in seq_len(fitted_tate$n_folds)) {
    mu1 <- fitted_tate$arm_results$mu1$intermediates$fold_info[[fold_index]]
    mu0 <- fitted_tate$arm_results$mu0$intermediates$fold_info[[fold_index]]
    tate <- fitted_tate$intermediates$fold_info[[fold_index]]
    if (!identical(mu1$target_idx, mu0$target_idx) || !identical(mu1$source_idx, mu0$source_idx) ||
        !identical(mu1$source_idx, tate$source_idx)) stop("Two-arm fold observation order differs")
    # The saved difference variance recovers the paired local covariance;
    # this identity is also checked against the retained residual vectors.
    fold_covariances[fold_index, ] <- (mu1$V_s + mu0$V_s - tate$V_s) / 2
    direct <- vapply(seq_along(sources), function(index)
      mean(mu1$xi_components[[index]] * mu0$xi_components[[index]]), numeric(1L))
    if (max(abs(direct - fold_covariances[fold_index, ])) > 1e-12) {
      stop("Saved local TATE variance disagrees with the cross-arm residual covariance")
    }
  }
  pooled_cross <- pool_source_fold_cross_covariance(
    arms$mu1$source_message_moments$means, arms$mu0$source_message_moments$means,
    fold_covariances, arms$mu1$source_message_moments$sizes)
  message_error <- 0
  for (index in seq_along(sources)) {
    scores <- cbind(mu1 = arms$mu1$score_records$source[[index]],
                    mu0 = arms$mu0$score_records$source[[index]])
    private <- score_mean_covariance(scores)
    positions <- c(index + 1L, width + index + 1L)
    covariance[positions, positions] <- covariance[positions, positions] + private
    source_covariances[[index]] <- private
    message_error <- max(message_error, abs(private[1L, 2L] - pooled_cross[index]))
  }
  if (message_error > 1e-12) stop("Pooled source cross-arm messages disagree with individual scores")
  means <- cbind(mu1 = arms$mu1$estimates, mu0 = arms$mu0$estimates)
  contrast <- cbind(diag(width), -diag(width))
  difference_covariance <- contrast %*% covariance %*% t(contrast)
  dimnames(difference_covariance) <- list(candidates, candidates)
  reference <- extract_candidate_score_summary(fitted_tate, data_split)
  estimate_error <- max(abs(means[, "mu1"] - means[, "mu0"] - reference$estimates))
  covariance_error <- max(abs(difference_covariance - reference$covariance))
  if (estimate_error > 1e-12 || covariance_error > 1e-12) {
    stop("Joint arm summaries do not reproduce the saved TATE candidate moments")
  }
  result <- list(means = means, covariance = covariance,
    target_covariance = target_covariance, source_private_covariances = source_covariances,
    source_fold_cross_covariances = fold_covariances, site_sizes = arms$mu1$site_sizes,
    estimate_identity_error = estimate_error, covariance_identity_error = covariance_error,
    source_message_identity_error = message_error,
    inference_scope = "Empirical joint arm score moments; normal approximation and nuisance remainders are not established by this extraction")
  if (keep_scores) result$score_records <- list(target = target_scores,
    source = Map(function(mu1, mu0) cbind(mu1 = mu1, mu0 = mu0),
      arms$mu1$score_records$source, arms$mu0$score_records$source))
  result
}
