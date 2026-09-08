# Experimental common-TATE weights. Not the current RoCE default.
# Objective: n_target * Var_hat(eta) + n_target^power * sum(delta_hat^2 * eta^2).
.quadratic_bias_weights <- function(moments, power = .75) {
  required <- c("V_ot", "V_t", "V_s", "C_ot", "C_cross", "n_t", "n_s",
                "avg_target_est", "avg_source_est")
  if (!is.list(moments) || !all(required %in% names(moments))) stop("incomplete weight moments")
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= .5 || power >= 1) stop("power must lie strictly between 1/2 and 1")
  K <- length(moments$V_t)
  for (field in c("V_ot", "n_t", "avg_target_est")) {
    x <- moments[[field]]
    if (!is.numeric(x) || length(x) != 1L || !is.null(dim(x)) || !is.finite(x)) {
      stop("invalid scalar weight moments: ", field)
    }
  }
  if (moments$V_ot < 0 || moments$n_t < 2 || moments$n_t != floor(moments$n_t) || K < 1L) {
    stop("invalid scalar weight moments")
  }
  for (field in c("V_t", "V_s", "C_ot", "n_s", "avg_source_est")) {
    x <- moments[[field]]
    if (!is.numeric(x) || length(x) != K || !is.null(dim(x)) || any(!is.finite(x))) {
      stop("invalid vector weight moments: ", field)
    }
  }
  if (any(moments$V_t < 0) || any(moments$V_s < 0) || any(moments$n_s < 2) ||
      any(moments$n_s != floor(moments$n_s))) stop("invalid variance or site sample size")
  C <- moments$C_cross
  if (!is.matrix(C) || !is.numeric(C) || !identical(dim(C), c(K, K)) || any(!is.finite(C))) {
    stop("invalid shared-target covariance matrix")
  }
  scale <- max(abs(c(C, moments$V_ot, moments$V_t, moments$V_s)))
  if (max(abs(C-t(C))) > 100*.Machine$double.eps*K*max(scale, .Machine$double.xmin)) {
    stop("shared-target covariance must be symmetric")
  }
  C <- (C+t(C))/2
  diag(C) <- moments$V_t
  ones <- rep(1, K)
  A <- moments$V_ot*tcrossprod(ones)-outer(moments$C_ot, ones)-
    outer(ones, moments$C_ot)+C+diag(moments$n_t*moments$V_s/moments$n_s, K)
  A <- (A+t(A))/2
  # A is the target-n-scaled variance of source-minus-target contrasts.
  # Require positive curvature; never add an unreported statistical ridge.
  if (any(diag(A) <= 0)) stop("contrast covariance lacks positive diagonal curvature")
  A_scaled <- A/outer(sqrt(diag(A)), sqrt(diag(A)))
  if (inherits(try(chol(A_scaled), silent = TRUE), "try-error")) {
    stop("contrast covariance must be positive definite")
  }
  discrepancy <- moments$avg_source_est-moments$avg_target_est
  penalty <- moments$n_t^power*discrepancy^2
  H <- A+diag(penalty, K)
  linear <- moments$C_ot-moments$V_ot
  diagonal_scale <- sqrt(diag(H))
  normalized <- H/outer(diagonal_scale, diagonal_scale)
  cholesky <- chol(normalized)
  normalized_weights <- backsolve(cholesky,
    forwardsolve(t(cholesky), -linear/diagonal_scale))
  weights <- drop(normalized_weights/diagonal_scale)
  residual <- max(abs(drop(H %*% weights)+linear))
  relative_residual <- residual/max(max(abs(linear)), max(abs(H))*max(abs(weights)),
                                    .Machine$double.xmin)
  if (any(!is.finite(weights)) || !is.finite(relative_residual) || relative_residual > 1e-10) {
    stop("quadratic bias-weight normal equations failed")
  }
  list(weights = weights, discrepancy = discrepancy, penalty = penalty, power = power,
       target_training_n = moments$n_t, normal_equation_error = residual,
       relative_normal_equation_error = relative_residual,
       variance_quadratic = A, linear = linear)
}

.common_tate_from_weights <- function(fit, weights) {
  info <- fit$intermediates$fold_info
  sizes <- fit$intermediates$sample_sizes
  source_sizes <- sizes$n_source
  K <- length(source_sizes)
  if (!is.matrix(weights) || !identical(dim(weights), c(length(info), K)) ||
      any(!is.finite(weights))) stop("invalid common TATE fold weights")
  site_sizes <- c(sizes$n_t, source_sizes)
  N <- sum(site_sizes)
  values <- counts <- lapply(site_sizes, function(n) numeric(n))
  for (k in seq_along(info)) {
    entry <- info[[k]]; eta <- weights[k, ]
    target <- (1-sum(eta))*(entry$varphi_ot+entry$fold_target_estimate)
    for (j in seq_len(K)) {
      target <- target+eta[j]*(entry$zeta_components[[j]]+entry$mu_pred_ts[j])
      index <- entry$source_idx[[j]]
      values[[j+1L]][index] <- N/source_sizes[j]*eta[j]*(entry$xi_components[[j]]+entry$delta_ts[j])
      counts[[j+1L]] <- counts[[j+1L]]+tabulate(index, nbins = source_sizes[j])
    }
    values[[1L]][entry$target_idx] <- N/sizes$n_t*target
    counts[[1L]] <- counts[[1L]]+tabulate(entry$target_idx, nbins = sizes$n_t)
  }
  if (!all(unlist(counts) == 1) || any(!is.finite(unlist(values)))) stop("invalid observation-level TATE reconstruction")
  list(estimate = mean(unlist(values)),
       variance = sum(vapply(values, function(x) sum((x-mean(x))^2), numeric(1)))/N^2)
}
