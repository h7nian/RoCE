# Analytic empirical-measure derivative of the quadratic weight candidate.
# Nuisance predictions and fold membership are FIXED in this module. This is
# not a claim of a complete fitted-estimator IF or a resampling algorithm.

.quadratic_inner_records <- function(inner) {
  components <- inner$components
  K <- length(inner$V_s)
  target <- lapply(components, function(x) {
    cbind(x$varphi_ot+x$avg_target_est,
          sweep(do.call(cbind, x$zeta_components), 2L, x$mu_pred_ts, "+"))
  })
  source <- lapply(seq_len(K), function(j) lapply(components, function(x)
    x$xi_components[[j]]+x$delta_ts[j]))
  list(target = target, source = source,
       target_ids = lapply(components, `[[`, "target_idx"),
       source_ids = lapply(seq_len(K), function(j) lapply(components, function(x) x$source_idx[[j]])),
       n_t = inner$n_t, n_s = inner$n_s)
}

.quadratic_inner_moments <- function(records, mass) {
  target_mass <- lapply(records$target_ids, function(ids) mass[[1L]][ids])
  weighted_center <- function(x, w) {
    if (is.matrix(x)) sweep(x, 2L, colSums(x*w)/sum(w), "-") else x-sum(w*x)/sum(w)
  }
  centered_target <- do.call(rbind, Map(weighted_center, records$target, target_mass))
  raw_target <- do.call(rbind, records$target)
  tw <- unlist(target_mass, use.names = FALSE)
  target_mean <- colSums(raw_target*tw)/sum(tw)
  covariance <- crossprod(centered_target, centered_target*tw)/sum(tw)
  K <- length(records$source)
  source_mean <- source_variance <- numeric(K)
  centered_source <- source_raw <- source_w <- vector("list", K)
  for (j in seq_len(K)) {
    sw <- lapply(records$source_ids[[j]], function(ids) mass[[j+1L]][ids])
    centered_source[[j]] <- unlist(Map(weighted_center, records$source[[j]], sw), use.names = FALSE)
    source_raw[[j]] <- unlist(records$source[[j]], use.names = FALSE)
    source_w[[j]] <- unlist(sw, use.names = FALSE)
    source_mean[j] <- sum(source_w[[j]]*source_raw[[j]])/sum(source_w[[j]])
    source_variance[j] <- sum(source_w[[j]]*centered_source[[j]]^2)/sum(source_w[[j]])
  }
  C <- covariance[-1L, -1L, drop = FALSE]
  vt <- diag(C); diag(C) <- 0
  list(moments = list(V_ot = covariance[1L, 1L], V_t = vt, V_s = source_variance,
    C_ot = covariance[1L, -1L], C_cross = C, n_t = records$n_t, n_s = records$n_s,
    avg_target_est = target_mean[1L], avg_source_est = target_mean[-1L]+source_mean,
    mu_pred_ts = target_mean[-1L], delta_ts = source_mean),
    target_centered = centered_target, target_raw = raw_target,
    source_centered = centered_source, source_raw = source_raw)
}

.quadratic_outer_values <- function(fit, weights) {
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  values <- counts <- lapply(sizes, function(n) numeric(n))
  increments <- matrix(0, nrow(weights), ncol(weights))
  for (k in seq_len(nrow(weights))) {
    info <- fit$intermediates$fold_info[[k]]
    raw_target <- info$varphi_ot+info$fold_target_estimate
    eta <- weights[k, ]
    target <- (1-sum(eta))*raw_target
    for (j in seq_len(ncol(weights))) {
      raw_prediction <- info$zeta_components[[j]]+info$mu_pred_ts[j]
      raw_source <- info$xi_components[[j]]+info$delta_ts[j]
      target <- target+eta[j]*raw_prediction
      values[[j+1L]][info$source_idx[[j]]] <- eta[j]*raw_source
      counts[[j+1L]] <- counts[[j+1L]]+tabulate(info$source_idx[[j]], nbins = sizes[j+1L])
      increments[k, j] <- sum(raw_prediction-raw_target)/sizes[1L]+sum(raw_source)/sizes[j+1L]
    }
    values[[1L]][info$target_idx] <- target
    counts[[1L]] <- counts[[1L]]+tabulate(info$target_idx, nbins = sizes[1L])
  }
  if (!all(unlist(counts) == 1) || any(!is.finite(unlist(values)))) stop("outer observations do not form a finite partition")
  list(values = values, increments = increments, sizes = sizes)
}

.quadratic_weight_functional <- function(fit, mass = NULL) {
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  if (is.null(mass)) mass <- lapply(sizes, function(n) rep(1, n))
  if (!is.list(mass) || length(mass) != length(sizes) ||
      any(lengths(mass) != sizes) || any(!is.finite(unlist(mass))) || any(unlist(mass) <= 0)) {
    stop("empirical masses must be positive, finite and match site sizes")
  }
  inner <- fit$intermediates$inner_fold_info
  K <- length(sizes)-1L
  answers <- lapply(inner, function(x) {
    records <- .quadratic_inner_records(x)
    .quadratic_bias_weights(.quadratic_inner_moments(records, mass)$moments, .75)
  })
  weights <- do.call(rbind, lapply(answers, `[[`, "weights"))
  stopifnot(identical(dim(weights), c(length(inner), K)))
  outer <- .quadratic_outer_values(fit, weights)
  sum(vapply(seq_along(sizes), function(s) sum(mass[[s]]*outer$values[[s]])/sum(mass[[s]]), numeric(1)))
}

.quadratic_weight_influence <- function(fit) {
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  mass <- lapply(sizes, function(n) rep(1, n))
  inner <- fit$intermediates$inner_fold_info
  records <- lapply(inner, .quadratic_inner_records)
  centered <- lapply(records, .quadratic_inner_moments, mass = mass)
  for (k in seq_along(inner)) {
    for (field in c("V_ot", "V_t", "V_s", "C_ot", "C_cross", "avg_target_est", "avg_source_est")) {
      stopifnot(max(abs(inner[[k]][[field]]-centered[[k]]$moments[[field]])) < 1e-10)
    }
  }
  answers <- lapply(centered, function(x) .quadratic_bias_weights(x$moments, .75))
  weights <- do.call(rbind, lapply(answers, `[[`, "weights"))
  outer <- .quadratic_outer_values(fit, weights)
  direct <- lapply(seq_along(sizes), function(s) (outer$values[[s]]-mean(outer$values[[s]]))/sizes[s])
  indirect <- lapply(sizes, function(n) numeric(n))
  for (k in seq_along(inner)) {
    data <- centered[[k]]; m <- data$moments; answer <- answers[[k]]
    eta <- answer$weights; delta <- answer$discrepancy
    H <- answer$variance_quadratic+diag(answer$penalty, length(eta))
    # L is the exact global estimator derivative with respect to this fold's
    # common source weights; it includes different target/source fold fractions.
    h <- drop(solve(H, outer$increments[k, ]))
    A_target <- answer$variance_quadratic-diag(m$n_t*m$V_s/m$n_s, length(eta))
    q <- data$target_centered[, -1L, drop = FALSE]-data$target_centered[, 1L]
    qh <- drop(q %*% h); qe <- drop(q %*% eta)
    raw_difference <- data$target_raw[, -1L, drop = FALSE]-data$target_raw[, 1L]
    mean_difference <- m$mu_pred_ts-m$avg_target_est
    delta_derivative <- sweep(raw_difference, 2L, mean_difference, "-")
    nt <- m$n_t; lambda <- nt^.75
    target_term <- -(data$target_centered[, 1L]*qh-sum(h*answer$linear)+
      qh*qe-sum(h*drop(A_target %*% eta))+
      2*lambda*drop(delta_derivative %*% (h*delta*eta)))/nt
    ids <- unlist(records[[k]]$target_ids, use.names = FALSE)
    stopifnot(!anyDuplicated(ids), !any(ids %in% fit$intermediates$fold_info[[k]]$target_idx))
    indirect[[1L]][ids] <- indirect[[1L]][ids]+target_term
    for (j in seq_along(eta)) {
      ids <- unlist(records[[k]]$source_ids[[j]], use.names = FALSE)
      stopifnot(!anyDuplicated(ids), !any(ids %in% fit$intermediates$fold_info[[k]]$source_idx[[j]]))
      ns <- m$n_s[j]
      source_term <- -h[j]*eta[j]*(nt/ns^2*(data$source_centered[[j]]^2-m$V_s[j])+
        2*lambda*delta[j]/ns*(data$source_raw[[j]]-m$delta_ts[j]))
      indirect[[j+1L]][ids] <- indirect[[j+1L]][ids]+source_term
    }
  }
  full <- Map(`+`, direct, indirect)
  stopifnot(max(abs(vapply(full, sum, numeric(1)))) < 1e-10,
            max(abs(vapply(indirect, sum, numeric(1)))) < 1e-10)
  fixed_variance <- sum(unlist(direct)^2)
  indirect_variance <- sum(unlist(indirect)^2)
  cross_term <- 2*sum(unlist(direct)*unlist(indirect))
  variance <- sum(unlist(full)^2)
  baseline <- .common_tate_from_weights(fit, weights)
  estimate <- sum(vapply(outer$values, mean, numeric(1)))
  stopifnot(abs(baseline$variance-fixed_variance) < 1e-12,
            abs(baseline$estimate-estimate) < 1e-12,
            abs(variance-fixed_variance-indirect_variance-cross_term) < 1e-12)
  list(estimate = estimate, weights = weights, direct = direct, indirect = indirect,
       gradient = full, fixed_variance = fixed_variance,
       weight_linearized_variance = variance, indirect_variance = indirect_variance,
       direct_indirect_cross_term = cross_term, nuisance_fits_differentiated = FALSE)
}
