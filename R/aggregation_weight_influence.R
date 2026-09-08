# Delta-method contribution of the learned common source weights to the TATE
# variance: main.tex sec:adaptive_aggregation ("Variance of the learned
# weights") and supplemental.tex supp:weight_layer. Nuisance fits and fold
# membership are held fixed; only the inner-fold moments that determine the
# weights are differentiated with respect to the empirical mass of each
# observation.

# Raw influence records of the inner folds of one outer fold, keyed so that a
# target row holds (varphi_ot + tau_ot, zeta_j + mu_pred_j) and a source-j
# entry holds xi_j + delta_j, exactly the quantities averaged by
# .average_aggregation_components().
.inner_fold_records <- function(inner_info) {
  components <- inner_info$components
  K <- length(inner_info$V_s)
  target <- lapply(components, function(x) {
    cbind(x$varphi_ot + x$avg_target_est,
          sweep(do.call(cbind, x$zeta_components), 2L, x$mu_pred_ts, "+"))
  })
  source <- lapply(seq_len(K), function(j) {
    lapply(components, function(x) x$xi_components[[j]] + x$delta_ts[j])
  })
  list(
    target = target,
    source = source,
    target_ids = lapply(components, `[[`, "target_idx"),
    source_ids = lapply(seq_len(K), function(j) {
      lapply(components, function(x) x$source_idx[[j]])
    }),
    n_t = inner_info$n_t,
    n_s = inner_info$n_s
  )
}

# Pooled inner-fold moments as functionals of observation masses (unit mass
# reproduces the stored Phase-1b components). `mass` holds one vector per site
# indexed by the original site observation, so its lengths are the full site
# sizes, not the inner training sizes. Centering is within inner fold;
# covariances are pooled across inner folds with the total mass.
.inner_fold_moments <- function(records, mass) {
  K <- length(records$source)
  weighted_center <- function(x, w) {
    if (is.matrix(x)) sweep(x, 2L, colSums(x * w) / sum(w), "-") else x - sum(w * x) / sum(w)
  }
  target_mass <- lapply(records$target_ids, function(ids) mass[[1L]][ids])
  target_centered <- do.call(rbind, Map(weighted_center, records$target, target_mass))
  target_raw <- do.call(rbind, records$target)
  target_weight <- unlist(target_mass, use.names = FALSE)
  target_mean <- colSums(target_raw * target_weight) / sum(target_weight)
  covariance <- crossprod(target_centered, target_centered * target_weight) / sum(target_weight)
  source_mean <- source_variance <- numeric(K)
  source_centered <- source_raw <- vector("list", K)
  for (j in seq_len(K)) {
    source_mass <- lapply(records$source_ids[[j]], function(ids) mass[[j + 1L]][ids])
    source_centered[[j]] <- unlist(
      Map(weighted_center, records$source[[j]], source_mass), use.names = FALSE
    )
    source_raw[[j]] <- unlist(records$source[[j]], use.names = FALSE)
    weight <- unlist(source_mass, use.names = FALSE)
    source_mean[j] <- sum(weight * source_raw[[j]]) / sum(weight)
    source_variance[j] <- sum(weight * source_centered[[j]]^2) / sum(weight)
  }
  cross <- covariance[-1L, -1L, drop = FALSE]
  V_t <- diag(cross)
  diag(cross) <- 0
  list(
    moments = list(
      V_ot = covariance[1L, 1L], V_t = V_t, V_s = source_variance,
      C_ot = covariance[1L, -1L], C_cross = cross,
      n_t = records$n_t, n_s = records$n_s,
      avg_target_est = target_mean[1L],
      avg_source_est = target_mean[-1L] + source_mean,
      mu_pred_ts = target_mean[-1L], delta_ts = source_mean
    ),
    target_centered = target_centered, target_raw = target_raw,
    source_centered = source_centered, source_raw = source_raw
  )
}

# Derivative of the aggregated estimate with respect to each outer fold's
# weights: the fold-k1 discrepancy tau_j - tau_ot with the actual target and
# source fold fractions (n_folds x K).
.outer_fold_increments <- function(fold_info, n_t, n_source_full) {
  K <- length(n_source_full)
  do.call(rbind, lapply(fold_info, function(info) {
    raw_target <- info$varphi_ot + info$fold_target_estimate
    vapply(seq_len(K), function(j) {
      sum(info$zeta_components[[j]] + info$mu_pred_ts[j] - raw_target) / n_t +
        sum(info$xi_components[[j]] + info$delta_ts[j]) / n_source_full[j]
    }, numeric(1L))
  }))
}

# Quadratic form of N_all * Var(eta) = N_all * (c + 2 l'eta + eta' Q eta); the
# optimizer's PSD ridge enters the curvature as 2 Q + ridge I.
.variance_quadratic_form <- function(m, psd_ridge) {
  K <- length(m$V_t)
  Q <- matrix(m$V_ot, K, K) - outer(m$C_ot, rep(1, K)) - outer(rep(1, K), m$C_ot)
  cross <- m$C_cross
  diag(cross) <- 0
  Q <- (Q + cross) / m$n_t
  diag(Q) <- diag(Q) + m$V_t / m$n_t + m$V_s / m$n_s + psd_ridge / 2
  list(Q = Q, l = (m$C_ot - m$V_ot) / m$n_t, N_all = m$n_t + sum(m$n_s))
}

# Wald ingredients of the penalty p_j = (lambda t_j - 1)_+, t_j = |d_j| / s_j,
# mirroring .compute_phase2_weights(): the discrepancy variance is floored at
# VARIANCE_MIN exactly as the optimizer floors it, and a floored coordinate has
# a constant standard error, hence a zero derivative.
.wald_penalty_ingredients <- function(m, lambda) {
  d <- m$avg_target_est - m$avg_source_est
  s2 <- m$V_ot / m$n_t + m$V_t / m$n_t + m$V_s / m$n_s - 2 * m$C_ot / m$n_t
  floored <- s2 < VARIANCE_MIN
  s <- sqrt(pmax(s2, VARIANCE_MIN))
  t <- abs(d) / s
  list(d = d, s = s, t = t, floored = floored,
       penalty_active = lambda * t > 1,
       kink = abs(lambda * t - 1) < WEIGHT_LAYER_KINK_TOLERANCE)
}

# Per-observation derivative of one outer fold's contribution to the estimate
# through its weights (supplemental.tex eq:supp_weight_derivative). `lambda` is
# the fold's Wald multiplier. Under "hard_threshold" the retained sources were
# solved without a penalty, so only the kink count uses the cutoff and the
# discrete inclusion event is not differentiated. Under "quadratic_bias" every
# coordinate is active and the smooth penalty n_t^(power-1) delta_j^2 enters
# the curvature and, through delta_j, the derivative.
.weight_rule_derivative <- function(data, eta, screening_rule, lambda, psd_ridge, increment) {
  m <- data$moments
  form <- .variance_quadratic_form(m, psd_ridge)
  K <- length(eta)
  quadratic <- identical(screening_rule, "quadratic_bias")
  penalized <- identical(screening_rule, "soft_penalty")
  wald <- .wald_penalty_ingredients(m, lambda)
  kinks <- if (quadratic) 0L else sum(wald$kink)
  discrepancy <- m$avg_source_est - m$avg_target_est
  penalty_rate <- if (quadratic) m$n_t^(AGG_QUADRATIC_BIAS_POWER - 1) else 0
  penalty_curvature <- penalty_rate * discrepancy^2
  active <- if (quadratic) seq_len(K) else which(abs(eta) > WEIGHT_LAYER_ACTIVE_TOLERANCE)
  target_gradient <- numeric(nrow(data$target_raw))
  source_gradient <- lapply(m$n_s, function(n) numeric(n))
  if (length(active) == 0L) {
    return(list(target = target_gradient, source = source_gradient, kinks = kinks))
  }
  sign_eta <- sign(eta[active])
  penalty <- if (penalized) pmax(lambda * wald$t[active] - 1, 0) else rep(0, length(active))
  Q_active <- form$Q[active, active, drop = FALSE] +
    diag(penalty_curvature[active], length(active))
  stationarity <- Q_active %*% eta[active] + form$l[active] +
    sign_eta * penalty / (2 * form$N_all)
  stationarity_scale <- max(abs(Q_active), abs(form$l[active]), penalty / (2 * form$N_all))
  if (max(abs(stationarity)) > WEIGHT_LAYER_KKT_TOLERANCE * stationarity_scale) {
    stop(sprintf(
      ".weight_rule_derivative: fold weights violate the stationarity conditions (residual %.3e, scale %.3e).",
      max(abs(stationarity)), stationarity_scale
    ), call. = FALSE)
  }
  h <- drop(solve(Q_active, increment[active]))
  W_t <- m$n_t
  c0 <- data$target_centered[, 1L]
  cs <- data$target_centered[, -1L, drop = FALSE]
  y0 <- data$target_raw[, 1L]
  ys <- data$target_raw[, -1L, drop = FALSE]
  dV_ot <- (c0^2 - m$V_ot) / W_t
  dC_ot <- sweep(c0 * cs, 2L, m$C_ot, "-") / W_t
  dV_t <- sweep(cs^2, 2L, m$V_t, "-") / W_t
  dmu_ot <- (y0 - m$avg_target_est) / W_t
  dmu_pred <- sweep(ys, 2L, m$mu_pred_ts, "-") / W_t
  penalty_scale <- if (penalized) lambda * wald$penalty_active / (2 * form$N_all) else rep(0, K)
  variance_slope <- ifelse(wald$floored, 0, 1)
  for (a in seq_along(active)) {
    j <- active[a]
    dQ_eta <- numeric(length(c0))
    for (k in active) {
      dC_cross <- if (j == k) 0 else (cs[, j] * cs[, k] - m$C_cross[j, k]) / W_t
      dQ_eta <- dQ_eta + eta[k] * (dV_ot - dC_ot[, j] - dC_ot[, k] + dC_cross) / m$n_t
    }
    dQ_eta <- dQ_eta + eta[j] * dV_t[, j] / m$n_t
    dl <- (dC_ot[, j] - dV_ot) / m$n_t
    # Target observations move the discrepancy through mu_ot and mu_pred_j.
    ddiscrepancy <- dmu_pred[, j] - dmu_ot
    penalty_term <- if (quadratic) {
      2 * penalty_rate * discrepancy[j] * ddiscrepancy * eta[j]
    } else {
      ds2 <- variance_slope[j] * (dV_ot + dV_t[, j] - 2 * dC_ot[, j]) / m$n_t
      dt <- sign(wald$d[j]) * (-ddiscrepancy) / wald$s[j] -
        abs(wald$d[j]) * ds2 / (2 * wald$s[j]^3)
      sign_eta[a] * penalty_scale[j] * dt
    }
    target_gradient <- target_gradient - h[a] * (dQ_eta + dl + penalty_term)
    # Source-j observations enter only through V_s,j and delta_j.
    dV_s <- (data$source_centered[[j]]^2 - m$V_s[j]) / m$n_s[j]
    ddelta <- (data$source_raw[[j]] - m$delta_ts[j]) / m$n_s[j]
    source_penalty_term <- if (quadratic) {
      2 * penalty_rate * discrepancy[j] * ddelta * eta[j]
    } else {
      dt_source <- -sign(wald$d[j]) * ddelta / wald$s[j] -
        abs(wald$d[j]) * variance_slope[j] * (dV_s / m$n_s[j]) / (2 * wald$s[j]^3)
      sign_eta[a] * penalty_scale[j] * dt_source
    }
    source_gradient[[j]] <- source_gradient[[j]] -
      h[a] * (eta[j] * dV_s / m$n_s[j] + source_penalty_term)
  }
  list(target = target_gradient, source = source_gradient, kinks = kinks)
}

# Weight-layer gradient of every observation, accumulated over the outer folds
# whose inner records contain it. An active PSD ridge depends on the moments
# through an eigenvalue and is not differentiated here, so it is rejected
# rather than treated as a constant.
.weight_layer_gradient <- function(fold_info, inner_fold_info, fold_weights,
                                   fold_lambdas, screening_rule, fold_weight_psd_ridge,
                                   n_t, n_source_full) {
  if (any(fold_weight_psd_ridge > 0)) {
    stop(".weight_layer_gradient: the weight-layer derivative is not available when the PSD ridge is active.",
         call. = FALSE)
  }
  increments <- .outer_fold_increments(fold_info, n_t, n_source_full)
  unit_mass <- c(list(rep(1, n_t)), lapply(n_source_full, function(n) rep(1, n)))
  target <- numeric(n_t)
  source <- lapply(n_source_full, function(n) numeric(n))
  kink_cells <- 0L
  for (k in seq_along(inner_fold_info)) {
    records <- .inner_fold_records(inner_fold_info[[k]])
    data <- .inner_fold_moments(records, unit_mass)
    for (field in c("V_ot", "V_t", "V_s", "C_ot", "C_cross", "avg_target_est", "avg_source_est")) {
      stored <- as.numeric(inner_fold_info[[k]][[field]])
      recomputed <- as.numeric(data$moments[[field]])
      if (length(stored) != length(recomputed) || !all(is.finite(c(stored, recomputed)))) {
        stop(sprintf(
          ".weight_layer_gradient: inner-fold %s of outer fold %d is missing or non-finite (stored: %s; recomputed: %s).",
          field, k, paste(signif(stored, 6), collapse = ","),
          paste(signif(recomputed, 6), collapse = ",")
        ), call. = FALSE)
      }
      difference <- max(abs(stored - recomputed))
      if (difference > WEIGHT_LAYER_MOMENT_TOLERANCE * max(1, max(abs(stored)))) {
        stop(sprintf(
          ".weight_layer_gradient: recomputed inner-fold %s of outer fold %d differs from the stored Phase-1b component by %.3e.",
          field, k, difference
        ), call. = FALSE)
      }
    }
    fold <- .weight_rule_derivative(
      data, fold_weights[k, ], screening_rule, fold_lambdas[k],
      fold_weight_psd_ridge[k], increments[k, ]
    )
    kink_cells <- kink_cells + fold$kinks
    ids <- unlist(records$target_ids, use.names = FALSE)
    if (anyDuplicated(ids) || any(ids %in% fold_info[[k]]$target_idx)) {
      stop(".weight_layer_gradient: inner target records overlap the outer evaluation fold.",
           call. = FALSE)
    }
    target[ids] <- target[ids] + fold$target
    for (j in seq_along(fold$source)) {
      ids <- unlist(records$source_ids[[j]], use.names = FALSE)
      if (anyDuplicated(ids) || any(ids %in% fold_info[[k]]$source_idx[[j]])) {
        stop(".weight_layer_gradient: inner source records overlap the outer evaluation fold.",
             call. = FALSE)
      }
      source[[j]][ids] <- source[[j]][ids] + fold$source[[j]]
    }
  }
  list(target = target, source = source, kink_cells = kink_cells,
       fold_source_cells = length(fold_weights))
}

# Variance of the mean pseudo-value with the weight layer added: the direct
# gradient is the within-site centered pseudo-value over N_all, so the
# fixed-weight part equals .multisite_pseudovalue_variance().
.weight_layer_variance <- function(all_phi_agg, n_t, n_source_full, gradient) {
  sizes <- c(n_t, n_source_full)
  N_all <- sum(sizes)
  ends <- cumsum(sizes)
  direct <- lapply(seq_along(sizes), function(s) {
    phi <- all_phi_agg[(ends[s] - sizes[s] + 1L):ends[s]]
    (phi - mean(phi)) / N_all
  })
  indirect <- c(list(gradient$target), gradient$source)
  total <- Map(`+`, direct, indirect)
  site_sums <- vapply(indirect, sum, numeric(1L))
  site_scales <- vapply(indirect, function(g) max(1, sum(abs(g))), numeric(1L))
  if (any(abs(site_sums) > WEIGHT_LAYER_MOMENT_TOLERANCE * site_scales)) {
    stop(".weight_layer_variance: weight-layer gradients must sum to zero within each site.",
         call. = FALSE)
  }
  direct_all <- unlist(direct, use.names = FALSE)
  indirect_all <- unlist(indirect, use.names = FALSE)
  list(
    variance = sum(unlist(total, use.names = FALSE)^2),
    fixed_variance = sum(direct_all^2),
    indirect_variance = sum(indirect_all^2),
    cross_term = 2 * sum(direct_all * indirect_all),
    kink_cells = gradient$kink_cells,
    fold_source_cells = gradient$fold_source_cells
  )
}
