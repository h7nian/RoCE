# cross_fitting_aggregation.R - Shared aggregation utilities for cross-fitting algorithms
# Extracts common code between run_crossfit communication modes (two_round/one_round)
# to follow DRY (Don't Repeat Yourself) principle.
#
# ============================================================================
# VERSION A: Outer k + Inner m cross-fitting with Sample-level Variance
# ============================================================================
# This implements the "Version A" aggregation scheme:
#   - Outer loop k=1..K_f: evaluation fold E_k, training set T_k
#   - Step 1: Nuisance θ̂^(-k) trained on T_k (reuses existing fold_results)
#   - Step 2: Inner M-fold on T_k estimates variance components for η
#   - Step 3: Weight optimization η̂^(-k) on T_k (independent of E_k)
#   - Step 4: Per-sample aggregated IF φ̂_{agg,i} on E_k
#   - Step 5: Group-centered sample-level variance (1/N²)Σ_g Σ_{i∈g}(φ̂_{agg,i} - φ̄_g)²
#
#' Compute a site-stratified pseudo-value variance
#'
#' Given observation-level pseudo-values from all outer folds (each
#' observation appearing exactly once), compute
#' \deqn{\operatorname{Var}(\bar\Phi)
#' =N^{-2}\sum_g\sum_{i\in g}(\Phi_i-\bar\Phi_g)^2,}
#' where groups are the target and each source site.
#'
#' Under the multi-sample (site-stratified) asymptotic regime (fixed sites,
#' independent samples per site with fixed sampling fractions), the variance of
#' the overall mean depends on within-site variation only; global centering would
#' add an ANOVA between-site term that does not correspond to sampling
#' variability.
#'
#' The caller computes the point estimate as the overall pseudo-value mean;
#' this helper supplies its matching within-site-centered variance.
#'
#' @param phi Finite numeric vector of observation-level pseudo-values, ordered
#'   contiguously by site.
#' @param group_sizes Non-negative integer vector giving the number of entries
#'   in \code{phi} contributed by each site.
#' @return The site-stratified variance of the overall pseudo-value mean.
#' @keywords internal
.multisite_pseudovalue_variance <- function(phi, group_sizes) {
  if (!is.numeric(phi) || any(!is.finite(phi))) {
    stop(".multisite_pseudovalue_variance: phi must be a finite numeric vector.")
  }
  if (!is.numeric(group_sizes) || any(!is.finite(group_sizes)) ||
      any(group_sizes < 0) || any(group_sizes != floor(group_sizes))) {
    stop(".multisite_pseudovalue_variance: group_sizes must be non-negative integers.")
  }
  group_sizes <- as.integer(group_sizes)
  N_all <- length(phi)
  if (sum(group_sizes) != N_all) {
    stop(".multisite_pseudovalue_variance: sum(group_sizes) must equal length(phi).")
  }
  if (N_all == 0L) return(0)

  wss <- 0
  start <- 1L
  for (g_size in group_sizes) {
    if (g_size > 0L) {
      end <- start + g_size - 1L
      phi_g <- phi[start:end]
      wss <- wss + stable_sum_kahan((phi_g - mean(phi_g))^2)
      start <- end + 1L
    }
  }
  max(wss / (N_all^2), 0)
}

.pseudovalue_assignment_counts <- function(indices, values, size, label) {
  if (length(size) != 1L || !is.finite(size) || size < 1L ||
      size != floor(size)) {
    stop(label, ": destination size must be one positive integer.",
         call. = FALSE)
  }
  if (!is.numeric(indices) || length(indices) != length(values) ||
      any(!is.finite(indices)) || any(indices != floor(indices)) ||
      any(indices < 1L) || any(indices > size) ||
      !is.numeric(values) || any(!is.finite(values))) {
    stop(
      label,
      ": indices and finite pseudo-values must have matching valid lengths.",
      call. = FALSE
    )
  }
  tabulate(as.integer(indices), nbins = as.integer(size))
}

.assemble_target_pseudovalues <- function(fold_info, n_t, caller) {
  target_values <- numeric(n_t)
  assignment_count <- integer(n_t)
  for (info in fold_info) {
    values <- info$varphi_ot + info$fold_target_estimate
    counts <- .pseudovalue_assignment_counts(
      info$target_idx, values, n_t, caller
    )
    target_values[info$target_idx] <- values
    assignment_count <- assignment_count + counts
  }
  if (any(assignment_count != 1L)) {
    stop(
      caller,
      ": target folds must cover every observation exactly once.",
      call. = FALSE
    )
  }
  target_values
}

aggregate_fold_estimates <- function(fold_aggregated_estimates, all_phi_agg,
                                     N_all, n_t, n_source_full,
                                     fold_weights, n_folds, verbose,
                                     fold_info, inner_fold_info,
                                     fold_lambdas, fold_weight_psd_ridge,
                                     screening_rule = c("soft_penalty", "hard_threshold", "quadratic_bias")) {
  screening_rule <- match.arg(screening_rule)
  # =========================================================================
  # UNIFIED ESTIMATING EQUATION (Pitfall 1 fix)
  # =========================================================================
  # The point estimate and variance both derive from the SAME set {Φ̂_i}:
  #
  #   μ̂_agg  = Φ̄ = (1/N) Σᵢ Φ̂_i
  #   Var(μ̂) = (1/N²) Σ_g Σ_{i∈g} (Φ̂_i - Φ̄_g)²
  #
  # where Φ̂_i are the UNCENTERED per-sample IF values. This guarantees
  #   μ̂_agg - μ = (1/N) Σᵢ Φ̂_i - μ = (1/N) Σᵢ (Φ̂_i - μ)
  # and the multi-sample variance estimator is a within-group sample variance
  # of {Φ̂_i} (groups = sites), scaled by 1/N^2.
  #
  # With balanced folds, Φ̄ equals the fold-average (1/K_f)Σ μ̂_agg,k.
  # The fold-average is retained for diagnostics.
  # =========================================================================
  final_estimate <- mean(all_phi_agg)  # = Φ̄ = μ̂_agg

  # Average weights for reporting
  K <- if (is.matrix(fold_weights)) ncol(fold_weights) else 0
  average_weights <- if (K > 0) colMeans(fold_weights) else numeric(0)

  # Sample-level multi-sample variance (site-stratified):
  #   Var = (1/N^2) * sum_g sum_{i in g} (Phi_i - mean_g(Phi))^2
  # where groups are target and each source site.
  if (length(all_phi_agg) != N_all) {
    stop("Length mismatch: all_phi_agg must have length N_all")
  }
  if (!is.numeric(n_t) || length(n_t) != 1 || n_t < 0) {
    stop("n_t must be a non-negative scalar")
  }
  if (!is.numeric(n_source_full)) {
    stop("n_source_full must be a numeric vector")
  }
  group_sizes <- c(as.integer(n_t), as.integer(n_source_full))
  if (sum(group_sizes) != N_all) {
    stop("Group sizes (n_t + sum(n_source_full)) must equal N_all")
  }

  # The reported variance adds the delta-method contribution of the learned
  # fold weights to the fixed-weight pseudo-value variance (main.tex
  # sec:adaptive_aggregation, "Variance of the learned weights").
  weight_layer <- .weight_layer_variance(
    all_phi_agg, n_t, n_source_full,
    .weight_layer_gradient(
      fold_info, inner_fold_info, fold_weights, fold_lambdas, screening_rule,
      fold_weight_psd_ridge, n_t, n_source_full
    )
  )
  final_variance <- weight_layer$variance
  variance_fixed_weights <- weight_layer$fixed_variance

  se <- sqrt(final_variance)
  ci_lower <- final_estimate - Z_ALPHA_05 * se
  ci_upper <- final_estimate + Z_ALPHA_05 * se

  # Diagnostic: fold-average estimate for comparison
  fold_avg_estimate <- mean(fold_aggregated_estimates)

  if (verbose) {
    cat(sprintf("   Unified estimating equation (N_all=%d):\n", N_all))
    cat(sprintf("      - Point estimate (Phi_bar): %.6f\n", final_estimate))
    cat(sprintf("      - Fold-average estimate:    %.6f  (diff=%.2e)\n",
                fold_avg_estimate, abs(final_estimate - fold_avg_estimate)))
    cat(sprintf("      - Variance (within-site, fixed weights): %.6f\n", variance_fixed_weights))
    cat(sprintf("      - Variance (with weight layer): %.6f (SE=%.4f)\n", final_variance, se))
    cat(paste0("   95% CI: [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n"))
  }

  return(list(
    final_estimate = final_estimate,
    final_variance = final_variance,
    variance_fixed_weights = variance_fixed_weights,
    se = se,
    se_fixed_weights = sqrt(variance_fixed_weights),
    weight_layer = weight_layer[c("indirect_variance", "cross_term", "kink_cells", "fold_source_cells")],
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    average_weights = average_weights
  ))
}

.compute_cross_covariance_matrix <- function(zeta_list, K) {
  if (K <= 1) {
    return(matrix(0, nrow = K, ncol = K))
  }

  valid <- vapply(zeta_list, function(x) !is.null(x), logical(1L))
  if (sum(valid) <= 1) {
    return(matrix(0, nrow = K, ncol = K))
  }

  zeta_matrix <- do.call(cbind, zeta_list[valid])
  n_fold <- nrow(zeta_matrix)
  if (is.null(n_fold) || n_fold <= 0) {
    C_valid <- matrix(0, nrow = sum(valid), ncol = sum(valid))
  } else {
    zeta_centered <- scale(zeta_matrix, center = TRUE, scale = FALSE)
    C_valid <- crossprod(zeta_centered) / n_fold
  }

  full_C <- matrix(0, nrow = K, ncol = K)
  valid_idx <- which(valid)
  full_C[valid_idx, valid_idx] <- C_valid
  diag(full_C) <- 0
  full_C
}

.compute_phase1_fold_info <- function(k1, source_sites, K, fold_results,
                                      target_folds, source_folds,
                                      M_tau_inference,
                                      family_int, link_int, A_val) {
  target_fold_k1 <- materialize_fold(target_folds, k1)
  target_only_k1 <- fold_results[[k1]]$target_only
  varphi_ot_k1 <- target_only_k1$varphi_ot

  fold_target_estimate <- target_only_k1$estimate
  fold_source_estimates <- numeric(K)
  mu_pred_ts_k1 <- numeric(K)
  delta_ts_k1 <- numeric(K)

  V_t_k1 <- numeric(K)
  V_s_k1 <- numeric(K)
  C_ot_k1 <- numeric(K)
  zeta_components_k1 <- vector("list", K)
  xi_components_k1 <- vector("list", K)
  source_variance_clip_diagnostics <- vector("list", K)
  source_idx_k1 <- vector("list", K)

  for (i in seq_along(source_sites)) {
    s <- source_sites[i]
    fold_params <- fold_results[[k1]]$source_results[[s]]
    source_fold_k1 <- materialize_fold(source_folds[[s]], k1)
    source_idx_k1[[i]] <- source_fold_k1$original_idx
    fold_source_estimates[i] <- fold_params$mu_ts
    mu_pred_ts_k1[i] <- fold_params$mu_pred_ts
    delta_ts_k1[i] <- fold_params$delta_ts

    source_var_result <- calculate_source_variance_cpp(
      source_fold_k1$Z_site, source_fold_k1$A, source_fold_k1$Y,
      fold_params$gamma_s, fold_params$alpha_ts, fold_params$delta_ts,
      source_fold_k1$W_outcome, M_tau_inference,
      family_int, link_int, A_val
    )
    V_s_k1[i] <- source_var_result$V_s
    xi_components_k1[[i]] <- source_var_result$xi_components
    source_variance_clip_diagnostics[[i]] <- source_var_result$clip_diagnostics

    target_var_result <- calculate_target_variance_cpp(
      target_fold_k1$W_outcome, fold_params$alpha_ts, fold_params$mu_pred_ts,
      family_int, link_int
    )
    V_t_k1[i] <- target_var_result$V_t
    zeta_components_k1[[i]] <- target_var_result$zeta_components

    C_ot_k1[i] <- calculate_covariance_term_cpp(
      varphi_ot_k1, target_var_result$zeta_components
    )
  }

  C_cross_k1 <- .compute_cross_covariance_matrix(zeta_components_k1, K)

  list(
    V_ot = target_only_k1$V_ot,
    V_t = V_t_k1,
    V_s = V_s_k1,
    C_ot = C_ot_k1,
    C_cross = C_cross_k1,
    fold_target_estimate = fold_target_estimate,
    fold_source_estimates = fold_source_estimates,
    mu_pred_ts = mu_pred_ts_k1,
    delta_ts = delta_ts_k1,
    varphi_ot = varphi_ot_k1,
    zeta_components = zeta_components_k1,
    xi_components = xi_components_k1,
    source_variance_clip_diagnostics = source_variance_clip_diagnostics,
    target_idx = target_fold_k1$original_idx,
    source_idx = source_idx_k1
  )
}

.centered_second_moment <- function(x) {
  if (length(x) == 0L) return(0)
  mean((x - mean(x))^2)
}

.pool_fold_pairwise_estimates <- function(
    mu_pred_matrix, delta_matrix, target_fold_sizes, source_fold_sizes) {
  mu_pred_matrix <- as.matrix(mu_pred_matrix)
  delta_matrix <- as.matrix(delta_matrix)
  source_fold_sizes <- as.matrix(source_fold_sizes)
  target_fold_sizes <- as.numeric(target_fold_sizes)
  if (!identical(dim(mu_pred_matrix), dim(delta_matrix)) ||
      !identical(dim(mu_pred_matrix), dim(source_fold_sizes)) ||
      length(target_fold_sizes) != nrow(mu_pred_matrix) ||
      any(!is.finite(mu_pred_matrix)) || any(!is.finite(delta_matrix)) ||
      any(!is.finite(target_fold_sizes)) ||
      any(!is.finite(source_fold_sizes)) ||
      any(target_fold_sizes <= 0) || any(source_fold_sizes <= 0)) {
    stop(paste0(
      ".pool_fold_pairwise_estimates: component matrices and positive fold ",
      "sizes must have matching dimensions."
    ))
  }
  vapply(seq_len(ncol(mu_pred_matrix)), function(j) {
    stats::weighted.mean(mu_pred_matrix[, j], target_fold_sizes) +
      stats::weighted.mean(delta_matrix[, j], source_fold_sizes[, j])
  }, numeric(1L))
}

.combine_tate_fold_info <- function(info1, info0, K) {
  if (!identical(info1$target_idx, info0$target_idx)) {
    stop(".combine_tate_fold_info: treated/control target evaluation indices differ.")
  }
  if (length(info1$source_idx) != K || length(info0$source_idx) != K) {
    stop(".combine_tate_fold_info: source index lists do not match K.")
  }
  for (j in seq_len(K)) {
    if (!identical(info1$source_idx[[j]], info0$source_idx[[j]])) {
      stop(sprintf(
        ".combine_tate_fold_info: treated/control source evaluation indices differ for source %d.",
        j
      ))
    }
  }

  varphi_tau <- info1$varphi_ot - info0$varphi_ot
  zeta_tau <- lapply(seq_len(K), function(j) {
    info1$zeta_components[[j]] - info0$zeta_components[[j]]
  })
  xi_tau <- lapply(seq_len(K), function(j) {
    info1$xi_components[[j]] - info0$xi_components[[j]]
  })

  V_t <- vapply(zeta_tau, .centered_second_moment, numeric(1L))
  V_s <- vapply(xi_tau, .centered_second_moment, numeric(1L))
  C_ot <- vapply(zeta_tau, function(z) {
    calculate_covariance_term_cpp(varphi_tau, z)
  }, numeric(1L))

  list(
    V_ot = .centered_second_moment(varphi_tau),
    V_t = V_t,
    V_s = V_s,
    C_ot = C_ot,
    C_cross = .compute_cross_covariance_matrix(zeta_tau, K),
    fold_target_estimate = info1$fold_target_estimate - info0$fold_target_estimate,
    fold_source_estimates = info1$fold_source_estimates - info0$fold_source_estimates,
    mu_pred_ts = info1$mu_pred_ts - info0$mu_pred_ts,
    delta_ts = info1$delta_ts - info0$delta_ts,
    varphi_ot = varphi_tau,
    zeta_components = zeta_tau,
    xi_components = xi_tau,
    source_variance_clip_diagnostics = list(
      treated = info1$source_variance_clip_diagnostics,
      control = info0$source_variance_clip_diagnostics
    ),
    target_idx = info1$target_idx,
    source_idx = info1$source_idx
  )
}

.compute_inner_source_arm_info <- function(
    source_fold, target_fold, gamma, alpha,
    M_tau_inference, family_int, link_int, A_val) {
  correction <- calculate_correction_term_cpp(
    source_fold$Z_site, source_fold$A, source_fold$Y,
    gamma, alpha, source_fold$W_outcome, M_tau_inference,
    family_int, link_int, A_val
  )
  delta <- correction$delta_ts
  mu_pred <- mean(predict_glm_cpp(
    target_fold$W_outcome, alpha, family_int, link_int
  ))
  source_variance <- calculate_source_variance_cpp(
    source_fold$Z_site, source_fold$A, source_fold$Y,
    gamma, alpha, delta, source_fold$W_outcome, M_tau_inference,
    family_int, link_int, A_val
  )
  target_variance <- calculate_target_variance_cpp(
    target_fold$W_outcome, alpha, mu_pred, family_int, link_int
  )

  list(
    estimate = mu_pred + delta,
    mu_pred = mu_pred,
    delta = delta,
    xi = source_variance$xi_components,
    zeta = target_variance$zeta_components,
    V_s = source_variance$V_s,
    V_t = target_variance$V_t,
    clip_diagnostics = source_variance$clip_diagnostics
  )
}

.compute_phase1b_inner_fold_info <- function(k1, n_folds, source_sites, K,
                                             fold_results, target_folds,
                                             source_folds, M_tau_inference,
                                             family_int, link_int, A_val) {
  secondary_folds <- setdiff(1:n_folds, k1)
  n_inner <- length(secondary_folds)
  components <- setNames(vector("list", n_inner), paste0("k2_", secondary_folds))

  for (k2 in secondary_folds) {
    k2_key <- paste0("k2_", k2)

    to_inner <- fold_results[[k1]]$target_only_inner[[k2_key]]
    varphi_ot_inner <- to_inner$varphi_ot
    V_ot_k2 <- to_inner$V_ot
    target_est_k2 <- to_inner$estimate

    target_fold_k2 <- materialize_fold(target_folds, k2)
    V_t_k2 <- numeric(K)
    V_s_k2 <- numeric(K)
    C_ot_k2 <- numeric(K)
    source_est_k2 <- numeric(K)
    mu_pred_k2 <- numeric(K)
    delta_k2 <- numeric(K)
    n_s_val_k2 <- numeric(K)
    zeta_list_inner <- vector("list", K)
    xi_list_inner <- vector("list", K)
    source_idx_inner <- stats::setNames(vector("list", K), source_sites)

    for (i in seq_along(source_sites)) {
      s <- source_sites[i]
      src <- fold_results[[k1]]$source_results[[s]]
      gamma_k2 <- src$per_k2_gamma[[k2_key]]
      alpha_k2 <- src$per_k2_alpha[[k2_key]]

      if (is.null(gamma_k2) || is.null(alpha_k2)) {
        missing_fields <- c(
          if (is.null(gamma_k2)) "per_k2_gamma" else character(0),
          if (is.null(alpha_k2)) "per_k2_alpha" else character(0)
        )
        stop(sprintf(
          ".aggregate_inner_training_components: missing %s for outer fold %s, inner fold %s, source '%s'.",
          paste(missing_fields, collapse = " and "), k1, k2, s
        ), call. = FALSE)
      }

      source_fold_k2 <- materialize_fold(source_folds[[s]], k2)
      n_s_val_k2[i] <- source_fold_k2$n
      source_idx_inner[[s]] <- source_fold_k2$original_idx

      arm_info <- .compute_inner_source_arm_info(
        source_fold = source_fold_k2,
        target_fold = target_fold_k2,
        gamma = gamma_k2,
        alpha = alpha_k2,
        M_tau_inference = M_tau_inference,
        family_int = family_int,
        link_int = link_int,
        A_val = A_val
      )
      source_est_k2[i] <- arm_info$estimate
      mu_pred_k2[i] <- arm_info$mu_pred
      delta_k2[i] <- arm_info$delta
      V_s_k2[i] <- arm_info$V_s
      V_t_k2[i] <- arm_info$V_t
      zeta_list_inner[[i]] <- arm_info$zeta
      xi_list_inner[[i]] <- arm_info$xi
      C_ot_k2[i] <- calculate_covariance_term_cpp(varphi_ot_inner, arm_info$zeta)
    }

    C_cross_k2 <- .compute_cross_covariance_matrix(zeta_list_inner, K)

    components[[k2_key]] <- list(
      outer_fold = as.integer(k1),
      inner_fold = as.integer(k2),
      target_idx = target_fold_k2$original_idx,
      source_idx = source_idx_inner,
      V_ot = V_ot_k2,
      V_t = V_t_k2,
      V_s = V_s_k2,
      C_ot = C_ot_k2,
      C_cross = C_cross_k2,
      avg_target_est = target_est_k2,
      avg_source_est = source_est_k2,
      mu_pred_ts = mu_pred_k2,
      delta_ts = delta_k2,
      n_t = target_fold_k2$n,
      n_s = n_s_val_k2,
      varphi_ot = varphi_ot_inner,
      zeta_components = zeta_list_inner,
      xi_components = xi_list_inner
    )
  }

  averaged <- .average_aggregation_components(components)
  averaged$components <- components
  averaged
}

.combine_tate_aggregation_component <- function(component1, component0, K) {
  required <- c(
    "avg_target_est", "avg_source_est", "n_t", "n_s",
    "mu_pred_ts", "delta_ts",
    "varphi_ot", "zeta_components", "xi_components"
  )
  for (arm_component in list(component1, component0)) {
    missing <- required[!vapply(required, function(name) {
      !is.null(arm_component[[name]])
    }, logical(1L))]
    if (length(missing) > 0L) {
      stop(sprintf(
        ".combine_tate_aggregation_component: missing field(s): %s.",
        paste(missing, collapse = ", ")
      ))
    }
  }
  if (!identical(component1$n_t, component0$n_t) ||
      !identical(component1$n_s, component0$n_s)) {
    stop(".combine_tate_aggregation_component: treated/control sample sizes differ.")
  }
  index_fields <- c("outer_fold", "inner_fold", "target_idx", "source_idx")
  missing_indices <- index_fields[!vapply(index_fields, function(field) {
    !is.null(component1[[field]]) && !is.null(component0[[field]])
  }, logical(1L))]
  if (length(missing_indices) > 0L) {
    stop(sprintf(
      ".combine_tate_aggregation_component: missing fold index field(s): %s.",
      paste(missing_indices, collapse = ", ")
    ))
  }
  if (!identical(component1$outer_fold, component0$outer_fold) ||
      !identical(component1$inner_fold, component0$inner_fold) ||
      !identical(component1$target_idx, component0$target_idx) ||
      !identical(component1$source_idx, component0$source_idx)) {
    stop(paste0(
      ".combine_tate_aggregation_component: treated/control inner-fold ",
      "observation indices differ."
    ))
  }

  varphi_tau <- component1$varphi_ot - component0$varphi_ot
  zeta_tau <- lapply(seq_len(K), function(j) {
    component1$zeta_components[[j]] - component0$zeta_components[[j]]
  })
  xi_tau <- lapply(seq_len(K), function(j) {
    component1$xi_components[[j]] - component0$xi_components[[j]]
  })

  list(
    outer_fold = component1$outer_fold,
    inner_fold = component1$inner_fold,
    target_idx = component1$target_idx,
    source_idx = component1$source_idx,
    V_ot = .centered_second_moment(varphi_tau),
    V_t = vapply(zeta_tau, .centered_second_moment, numeric(1L)),
    V_s = vapply(xi_tau, .centered_second_moment, numeric(1L)),
    C_ot = vapply(zeta_tau, function(zeta) {
      calculate_covariance_term_cpp(varphi_tau, zeta)
    }, numeric(1L)),
    C_cross = .compute_cross_covariance_matrix(zeta_tau, K),
    avg_target_est = component1$avg_target_est - component0$avg_target_est,
    avg_source_est = component1$avg_source_est - component0$avg_source_est,
    mu_pred_ts = component1$mu_pred_ts - component0$mu_pred_ts,
    delta_ts = component1$delta_ts - component0$delta_ts,
    n_t = component1$n_t,
    n_s = component1$n_s,
    varphi_ot = varphi_tau,
    zeta_components = zeta_tau,
    xi_components = xi_tau
  )
}

.combine_tate_inner_fold_info <- function(inner_info1, inner_info0, K) {
  component_names <- names(inner_info1$components)
  if (!identical(component_names, names(inner_info0$components))) {
    stop(".combine_tate_inner_fold_info: treated/control inner-fold names differ.")
  }
  components <- lapply(component_names, function(name) {
    .combine_tate_aggregation_component(
      inner_info1$components[[name]],
      inner_info0$components[[name]],
      K
    )
  })
  names(components) <- component_names

  averaged <- .average_aggregation_components(components)
  averaged$components <- components
  averaged
}

.average_aggregation_components <- function(components) {
  if (length(components) == 0L) {
    stop(".average_aggregation_components: components must be non-empty.")
  }

  raw_fields <- c(
    "mu_pred_ts", "delta_ts",
    "varphi_ot", "zeta_components", "xi_components"
  )
  raw_component_presence <- vapply(components, function(component) {
    all(vapply(raw_fields, function(field) {
      !is.null(component[[field]])
    }, logical(1L)))
  }, logical(1L))
  any_raw_field <- vapply(components, function(component) {
    any(raw_fields %in% names(component))
  }, logical(1L))
  if (any(any_raw_field) && !all(raw_component_presence)) {
    stop(paste0(
      ".average_aggregation_components: raw influence fields must be ",
      "present for either every component or none."
    ))
  }
  has_raw_components <- all(raw_component_presence)

  if (has_raw_components) {
    K <- length(components[[1L]]$zeta_components)
    if (any(vapply(components, function(component) {
      length(component$n_s) != K ||
        length(component$avg_source_est) != K ||
        length(component$mu_pred_ts) != K ||
        length(component$delta_ts) != K
    }, logical(1L)))) {
      stop(".average_aggregation_components: source vectors do not match K.")
    }
    n_t <- vapply(components, `[[`, numeric(1L), "n_t")
    n_s <- do.call(rbind, lapply(components, `[[`, "n_s"))
    if (ncol(n_s) != K || any(n_t <= 0) || any(n_s <= 0)) {
      stop(".average_aggregation_components: invalid raw-component sample sizes.")
    }
    for (component in components) {
      if (length(component$varphi_ot) != component$n_t ||
          length(component$zeta_components) != K ||
          length(component$xi_components) != K ||
          any(vapply(component$zeta_components, length, integer(1L)) !=
                component$n_t) ||
          any(vapply(component$xi_components, length, integer(1L)) !=
                component$n_s)) {
        stop(paste0(
          ".average_aggregation_components: raw influence components do not ",
          "match their target/source sample sizes."
        ))
      }
    }

    varphi_ot <- unlist(
      lapply(components, `[[`, "varphi_ot"), use.names = FALSE
    )
    zeta_components <- lapply(seq_len(K), function(j) {
      unlist(lapply(components, function(component) {
        component$zeta_components[[j]]
      }), use.names = FALSE)
    })
    xi_components <- lapply(seq_len(K), function(j) {
      unlist(lapply(components, function(component) {
        component$xi_components[[j]]
      }), use.names = FALSE)
    })

    target_estimates <- vapply(
      components, `[[`, numeric(1L), "avg_target_est"
    )
    mu_pred_ts <- do.call(
      rbind, lapply(components, `[[`, "mu_pred_ts")
    )
    delta_ts <- do.call(
      rbind, lapply(components, `[[`, "delta_ts")
    )
    return(list(
      V_ot = .centered_second_moment(varphi_ot),
      V_t = vapply(zeta_components, .centered_second_moment, numeric(1L)),
      V_s = vapply(xi_components, .centered_second_moment, numeric(1L)),
      C_ot = vapply(zeta_components, function(zeta) {
        calculate_covariance_term_cpp(varphi_ot, zeta)
      }, numeric(1L)),
      C_cross = .compute_cross_covariance_matrix(zeta_components, K),
      avg_target_est = stats::weighted.mean(target_estimates, n_t),
      avg_source_est = .pool_fold_pairwise_estimates(
        mu_pred_ts, delta_ts, n_t, n_s
      ),
      mu_pred_ts = vapply(seq_len(K), function(j) {
        stats::weighted.mean(mu_pred_ts[, j], n_t)
      }, numeric(1L)),
      delta_ts = vapply(seq_len(K), function(j) {
        stats::weighted.mean(delta_ts[, j], n_s[, j])
      }, numeric(1L)),
      n_t = sum(n_t),
      n_s = colSums(n_s),
      varphi_ot = varphi_ot,
      zeta_components = zeta_components,
      xi_components = xi_components
    ))
  }

  # Backward-compatible summary-only path used by lightweight validation
  # helpers. Production inner-fold components always take the exact raw path
  # above, which correctly handles fold-size rounding.
  list(
    V_ot = mean(vapply(components, `[[`, numeric(1L), "V_ot")),
    V_t = colMeans(do.call(rbind, lapply(components, `[[`, "V_t"))),
    V_s = colMeans(do.call(rbind, lapply(components, `[[`, "V_s"))),
    C_ot = colMeans(do.call(rbind, lapply(components, `[[`, "C_ot"))),
    C_cross = Reduce(`+`, lapply(components, `[[`, "C_cross")) / length(components),
    avg_target_est = mean(vapply(components, `[[`, numeric(1L), "avg_target_est")),
    avg_source_est = colMeans(do.call(rbind, lapply(components, `[[`, "avg_source_est"))),
    n_t = sum(vapply(components, `[[`, numeric(1L), "n_t")),
    n_s = colSums(do.call(rbind, lapply(components, `[[`, "n_s")))
  )
}

.validate_aggregation_grid_component <- function(component, caller) {
  required <- c("V_ot", "C_ot", "avg_target_est", "avg_source_est", "n_t")
  missing_required <- required[!vapply(required, function(field) {
    !is.null(component[[field]])
  }, logical(1L))]
  if (length(missing_required) > 0L) {
    stop(sprintf("%s: component is missing required field(s): %s.",
                 caller, paste(missing_required, collapse = ", ")),
         call. = FALSE)
  }

  estimates <- as.numeric(component$avg_source_est)
  C_ot <- as.numeric(component$C_ot)
  if (length(estimates) == 0L || length(C_ot) != length(estimates)) {
    stop(sprintf("%s: avg_source_est and C_ot must be non-empty vectors with the same length.",
                 caller), call. = FALSE)
  }
  if (!all(is.finite(estimates)) || !all(is.finite(C_ot))) {
    stop(sprintf("%s: avg_source_est and C_ot must be finite.", caller),
         call. = FALSE)
  }

  V_ot <- as.numeric(component$V_ot)
  mu_ot <- as.numeric(component$avg_target_est)
  n_t <- as.numeric(component$n_t)
  if (length(V_ot) != 1L || !is.finite(V_ot) ||
      length(mu_ot) != 1L || !is.finite(mu_ot) ||
      length(n_t) != 1L || !is.finite(n_t) || n_t <= 0) {
    stop(sprintf("%s: V_ot, avg_target_est, and n_t must be finite scalar values with n_t > 0.",
                 caller), call. = FALSE)
  }

  list(
    estimates = pmax(pmin(estimates, ESTIMATE_MAX), -ESTIMATE_MAX),
    C_ot = pmax(pmin(C_ot, ESTIMATE_MAX), -ESTIMATE_MAX),
    V_ot = max(V_ot, VARIANCE_MIN),
    mu_ot = max(min(mu_ot, ESTIMATE_MAX), -ESTIMATE_MAX),
    n_t = max(n_t, 1)
  )
}

.aggregation_lambda_max <- function(component) {
  fields <- .validate_aggregation_grid_component(component, ".aggregation_lambda_max")
  penalty_weights <- (fields$mu_ot - fields$estimates)^2
  grad0 <- 2 * (fields$C_ot - fields$V_ot) / fields$n_t

  penalized <- penalty_weights > 0
  unpenalized_active <- !penalized & abs(grad0) > sqrt(.Machine$double.eps)
  if (any(unpenalized_active)) {
    warning(sprintf(
      ".aggregation_lambda_max: %d source coordinate(s) have zero aggregation penalty weight and non-zero gradient at eta=0; no finite lambda can force those weights to zero. lambda_max is computed over penalized coordinates only.",
      sum(unpenalized_active)
    ), call. = FALSE)
  }
  if (!any(penalized)) {
    warning(".aggregation_lambda_max: all aggregation penalty weights are zero; lambda has no effect on the weighted-L1 penalty. Using LAMBDA_MIN as a degenerate path.",
            call. = FALSE)
    return(LAMBDA_MIN)
  }

  lambda_max <- max(abs(grad0[penalized]) / penalty_weights[penalized])
  if (!is.finite(lambda_max) || lambda_max <= 0) {
    warning(".aggregation_lambda_max: KKT lambda_max is non-positive or non-finite; using LAMBDA_MIN as a degenerate path.",
            call. = FALSE)
    return(LAMBDA_MIN)
  }
  lambda_max <- max(lambda_max, LAMBDA_MIN)
  if (lambda_max > LAMBDA_MAX) {
    warning(sprintf(".aggregation_lambda_max: exact KKT lambda_max=%g exceeds LAMBDA_MAX=%g; clipping the path endpoint.",
                    lambda_max, LAMBDA_MAX), call. = FALSE)
    lambda_max <- LAMBDA_MAX
  }
  lambda_max
}

.aggregation_lambda_grid <- function(component, lambda_grid = NULL) {
  if (!is.null(lambda_grid)) {
    if (!is.numeric(lambda_grid) || length(lambda_grid) == 0L ||
        any(!is.finite(lambda_grid)) || any(lambda_grid <= 0)) {
      stop(".aggregation_lambda_grid: lambda_grid must contain positive finite numeric values.",
           call. = FALSE)
    }
    clipped <- pmax(pmin(lambda_grid, LAMBDA_MAX), LAMBDA_MIN)
    if (!isTRUE(all.equal(clipped, lambda_grid, check.attributes = FALSE))) {
      warning(sprintf(".aggregation_lambda_grid: lambda_grid values were clipped to [%g, %g].",
                      LAMBDA_MIN, LAMBDA_MAX), call. = FALSE)
    }
    if (anyDuplicated(clipped)) {
      stop(
        paste0(
          ".aggregation_lambda_grid: lambda_grid values must remain unique ",
          "after clipping."
        ),
        call. = FALSE
      )
    }
    return(clipped)
  }

  # FACE truncated-Wald penalty factor (Han et al. 2023, Remark 9): lambda is the
  # inverse Wald penalty-activation threshold, not a variance-CV-tuned strength, so the
  # default "grid" is the single pilot-selected value AGG_WALD_LAMBDA (= 1,
  # i.e. activate its penalty when the discrepancy t-statistic exceeds 1).
  # Crossing this threshold need not produce an exact zero. A
  # variance grid + CV-min is inappropriate here (it would pick the penalty-
  # inactive region and retain variance-reducing biased sources); see
  # AGG_WALD_LAMBDA. .aggregation_lambda_max is retained for diagnostics/tests but
  # unused on this default path. A user-supplied lambda_grid (handled above) still
  # overrides for sensitivity analysis.
  grid <- AGG_WALD_LAMBDA
  attr(grid, "lambda_max") <- AGG_WALD_LAMBDA
  attr(grid, "lambda_min") <- AGG_WALD_LAMBDA
  grid
}

.validate_aggregation_lambda_grid <- function(
    lambda_selection, aggregation_lambda_grid, caller) {
  if (is.null(aggregation_lambda_grid)) {
    return(NULL)
  }
  if (!(is.character(lambda_selection) &&
        length(lambda_selection) == 1L &&
        identical(lambda_selection, "cv"))) {
    stop(
      sprintf(
        "%s: aggregation_lambda_grid is used only when lambda_selection = 'cv'.",
        caller
      ),
      call. = FALSE
    )
  }
  .aggregation_lambda_grid(
    component = NULL,
    lambda_grid = aggregation_lambda_grid
  )
}

.validation_aggregation_objective <- function(eta, component) {
  variances <- list(V_ot = component$V_ot, V_t = component$V_t, V_s = component$V_s)
  n_samples <- list(n_t = max(component$n_t, 1), n_s = pmax(component$n_s, 1))
  var_term <- calculate_aggregated_variance(
    eta, variances, component$C_ot, n_samples,
    C_cross = component$C_cross, lambda = 0, mu_ot = component$avg_target_est
  )
  n_val <- component$n_t + sum(component$n_s)
  # FACE Section 3.4: the validation criterion Q(eta) is the PURE validation
  # variance (N_V * Var) with NO penalty term. The aggregation penalty enters
  # only the training-fold optimization (eq:agg_penalized_objective); lambda selection picks
  # the validation-variance-minimizing lambda. Bias detection lives entirely in
  # the scale-free truncated-Wald penalty FACTOR (Han et al. 2023, Remark 9) used
  # inside optimize_weights, not in Q -- so Q stays pure variance per the paper.
  n_val * var_term
}

select_aggregation_lambda_inner_cv <- function(components, verbose = FALSE,
                                               lambda_grid = NULL,
                                               lambda_rule = c("min", "1se")) {
  lambda_rule <- match.arg(lambda_rule)
  if (length(components) < 2L) {
    stop("select_aggregation_lambda_inner_cv: at least two inner folds are required.")
  }

  all_component <- .average_aggregation_components(components)
  lambda_grid <- .aggregation_lambda_grid(all_component, lambda_grid)
  cv_scores <- rep(Inf, length(lambda_grid))
  cv_se <- rep(Inf, length(lambda_grid))
  n_valid_folds <- integer(length(lambda_grid))

  for (li in seq_along(lambda_grid)) {
    lambda <- lambda_grid[li]
    fold_scores <- numeric(length(components))

    for (m in seq_along(components)) {
      train_component <- .average_aggregation_components(components[-m])
      val_component <- components[[m]]

      variances_train <- list(
        V_ot = train_component$V_ot,
        V_t = train_component$V_t,
        V_s = train_component$V_s
      )
      n_train <- list(n_t = max(train_component$n_t, 1),
                      n_s = pmax(train_component$n_s, 1))

      eta_train <- tryCatch(
        optimize_weights(
          train_component$avg_source_est, variances_train, train_component$C_ot,
          n_train, lambda, train_component$avg_target_est,
          train_component$C_cross, clip_weights = FALSE
        ),
        error = function(e) {
          stop(sprintf(
            "select_aggregation_lambda_inner_cv: weight optimization failed for lambda index %d (lambda=%g), validation fold %d/%d: %s",
            li, lambda, m, length(components), conditionMessage(e)
          ), call. = FALSE)
        }
      )
      if (!all(is.finite(eta_train))) {
        stop(sprintf(
          "select_aggregation_lambda_inner_cv: weight optimization returned non-finite weights for lambda index %d (lambda=%g), validation fold %d/%d.",
          li, lambda, m, length(components)
        ), call. = FALSE)
      }

      fold_scores[m] <- .validation_aggregation_objective(eta_train, val_component)
      if (!is.finite(fold_scores[m])) {
        stop(sprintf(
          "select_aggregation_lambda_inner_cv: validation objective was non-finite for lambda index %d (lambda=%g), validation fold %d/%d.",
          li, lambda, m, length(components)
        ), call. = FALSE)
      }
    }

    n_valid_folds[li] <- length(fold_scores)
    cv_scores[li] <- mean(fold_scores)
    cv_se[li] <- stats::sd(fold_scores) / sqrt(length(fold_scores))
  }

  if (any(!is.finite(cv_scores))) {
    stop("select_aggregation_lambda_inner_cv: internal error; non-finite validation score escaped fail-fast validation.",
         call. = FALSE)
  }
  valid <- seq_along(lambda_grid)
  min_score <- min(cv_scores[valid])
  min_idx <- valid[which.min(cv_scores[valid])]
  if (identical(lambda_rule, "min")) {
    best_idx <- min_idx
  } else {
    one_se_cutoff <- min_score + cv_se[min_idx]
    within_1se <- valid[cv_scores[valid] <= one_se_cutoff]
    best_idx <- within_1se[which.max(lambda_grid[within_1se])]
  }
  idx_1se <- valid[cv_scores[valid] <= min_score + cv_se[min_idx]]
  idx_1se <- idx_1se[which.max(lambda_grid[idx_1se])]

  if (verbose) {
    cat(sprintf("  Inner validation min: %.6f at lambda=%.4f\n",
                min_score, lambda_grid[min_idx]))
    cat(sprintf("  Inner validation 1se: lambda=%.4f (cutoff=%.6f, se=%.6f)\n",
                lambda_grid[idx_1se], min_score + cv_se[min_idx], cv_se[min_idx]))
    cat(sprintf("  Selected: %.4f (score=%.6f, rule=%s)\n",
                lambda_grid[best_idx], cv_scores[best_idx], lambda_rule))
  }

  selected <- lambda_grid[best_idx]
  attr(selected, "lambda_min") <- lambda_grid[min_idx]
  attr(selected, "lambda_1se") <- lambda_grid[idx_1se]
  attr(selected, "idx_min") <- min_idx
  attr(selected, "idx_1se") <- idx_1se
  attr(selected, "cv_scores") <- cv_scores
  attr(selected, "cv_se") <- cv_se
  attr(selected, "n_valid_folds") <- n_valid_folds
  attr(selected, "lambda_rule") <- lambda_rule
  selected
}

.compute_phase2_weights <- function(n_folds, inner_fold_info, fold_info,
                                    lambda_selection, K,
                                    verbose, lambda_rule = c("min", "1se"),
                                    aggregation_lambda_grid = NULL,
                                    screening_rule = c(
                                      "soft_penalty", "hard_threshold", "quadratic_bias"
                                    )) {
  lambda_rule <- match.arg(lambda_rule)
  screening_rule <- match.arg(screening_rule)
  aggregation_lambda_grid <- .validate_aggregation_lambda_grid(
    lambda_selection,
    aggregation_lambda_grid,
    ".compute_phase2_weights"
  )
  fold_aggregated_estimates <- numeric(n_folds)
  fold_weights <- array(0, dim = c(n_folds, K))
  fold_lambdas <- numeric(n_folds)
  fold_wald_statistics <- array(NA_real_, dim = c(n_folds, K))
  fold_penalty_coefficients <- array(NA_real_, dim = c(n_folds, K))
  fold_source_included <- array(TRUE, dim = c(n_folds, K))
  fold_training_sample_sizes <- vector("list", n_folds)
  fold_weight_optimizer_iterations <- integer(n_folds)
  fold_weight_psd_ridge <- numeric(n_folds)
  fold_lambda_info <- vector("list", n_folds)

  for (k1 in 1:n_folds) {
    log_info(verbose, "      Phase 2: Inner-fold weights for fold k1=%d/%d", k1, n_folds)

    inner_v <- inner_fold_info[[k1]]
    variances_k1 <- list(V_ot = inner_v$V_ot, V_t = inner_v$V_t, V_s = inner_v$V_s)
    n_samples_k1 <- list(
      n_t = as.numeric(inner_v$n_t),
      n_s = as.numeric(inner_v$n_s)
    )
    if (length(n_samples_k1$n_t) != 1L ||
        !is.finite(n_samples_k1$n_t) || n_samples_k1$n_t <= 0 ||
        length(n_samples_k1$n_s) != K ||
        any(!is.finite(n_samples_k1$n_s)) || any(n_samples_k1$n_s <= 0)) {
      stop(sprintf(
        ".compute_phase2_weights: invalid inner-training sample sizes for outer fold %d.",
        k1
      ), call. = FALSE)
    }
    fold_training_sample_sizes[[k1]] <- n_samples_k1

    lambda_reg_k1 <- if (is.character(lambda_selection) &&
                         length(lambda_selection) == 1L &&
                         identical(lambda_selection, "cv")) {
      select_aggregation_lambda_inner_cv(
        inner_v$components, verbose = FALSE,
        lambda_grid = aggregation_lambda_grid,
        lambda_rule = lambda_rule
      )
    } else if (is.numeric(lambda_selection) &&
               length(lambda_selection) == 1L &&
               is.finite(lambda_selection)) {
      lambda_value <- .validate_lambda_scalar(lambda_selection, "calculate_crossfit_aggregation")
      if (lambda_value > LAMBDA_MAX) {
        stop(sprintf("calculate_crossfit_aggregation: lambda_selection must be <= %g.",
                     LAMBDA_MAX), call. = FALSE)
      }
      lambda_value
    } else {
      stop("calculate_crossfit_aggregation: lambda_selection must be 'cv' or a finite numeric scalar.",
           call. = FALSE)
    }
    fold_lambdas[k1] <- lambda_reg_k1
    discrepancy_variance <- inner_v$V_ot / n_samples_k1$n_t +
      inner_v$V_t / n_samples_k1$n_t +
      inner_v$V_s / n_samples_k1$n_s -
      2 * inner_v$C_ot / n_samples_k1$n_t
    discrepancy_se <- sqrt(pmax(discrepancy_variance, VARIANCE_MIN))
    fold_wald_statistics[k1, ] <-
      abs(inner_v$avg_target_est - inner_v$avg_source_est) /
      discrepancy_se
    fold_penalty_coefficients[k1, ] <- pmax(
      lambda_reg_k1 * fold_wald_statistics[k1, ] - 1,
      0
    )
    fold_lambda_info[[k1]] <- list(
      lambda_used = as.numeric(lambda_reg_k1),
      lambda_min = as.numeric(attr(lambda_reg_k1, "lambda_min") %||% lambda_reg_k1),
      lambda_1se = as.numeric(attr(lambda_reg_k1, "lambda_1se") %||% lambda_reg_k1),
      idx_min = attr(lambda_reg_k1, "idx_min") %||% NA_integer_,
      idx_1se = attr(lambda_reg_k1, "idx_1se") %||% NA_integer_,
      cv_scores = attr(lambda_reg_k1, "cv_scores") %||% numeric(0),
      cv_se = attr(lambda_reg_k1, "cv_se") %||% numeric(0),
      n_valid_folds = attr(lambda_reg_k1, "n_valid_folds") %||% integer(0),
      lambda_rule = attr(lambda_reg_k1, "lambda_rule") %||% "fixed"
    )

    if (identical(screening_rule, "hard_threshold")) {
      if (!is.finite(lambda_reg_k1) || lambda_reg_k1 <= 0) {
        stop(
          paste0(
            ".compute_phase2_weights: hard_threshold requires a strictly ",
            "positive aggregation lambda."
          ),
          call. = FALSE
        )
      }
      included <- fold_wald_statistics[k1, ] <= 1 / lambda_reg_k1
      fold_source_included[k1, ] <- included
      eta_k1 <- numeric(K)
      optimizer_iterations <- 1L
      optimizer_psd_ridge <- 0
      if (any(included)) {
        included_indices <- which(included)
        included_cross <- inner_v$C_cross[
          included_indices, included_indices, drop = FALSE
        ]
        included_weights <- optimize_weights(
          inner_v$avg_source_est[included_indices],
          list(
            V_ot = variances_k1$V_ot,
            V_t = variances_k1$V_t[included_indices],
            V_s = variances_k1$V_s[included_indices]
          ),
          inner_v$C_ot[included_indices],
          list(
            n_t = n_samples_k1$n_t,
            n_s = n_samples_k1$n_s[included_indices]
          ),
          lambda = 0,
          mu_ot = inner_v$avg_target_est,
          C_cross = included_cross,
          clip_weights = FALSE
        )
        eta_k1[included_indices] <- included_weights
        optimizer_iterations <-
          as.integer(attr(included_weights, "optimizer_iterations"))
        optimizer_psd_ridge <- as.numeric(attr(included_weights, "psd_ridge"))
      }
      attr(eta_k1, "optimizer_iterations") <- optimizer_iterations
      attr(eta_k1, "psd_ridge") <- optimizer_psd_ridge
    } else if (identical(screening_rule, "quadratic_bias")) {
      eta_k1 <- .quadratic_bias_weights(c(
        variances_k1,
        list(
          C_ot = inner_v$C_ot,
          C_cross = prepare_cross_matrix(inner_v$C_cross, K),
          n_t = n_samples_k1$n_t,
          n_s = n_samples_k1$n_s,
          avg_target_est = inner_v$avg_target_est,
          avg_source_est = inner_v$avg_source_est
        )
      ))
      # fold_wald_statistics and fold_penalty_coefficients stay the
      # truncated-Wald diagnostics; the quadratic rule does not use them.
      attr(eta_k1, "optimizer_iterations") <- 1L
      attr(eta_k1, "psd_ridge") <- 0
    } else {
      eta_k1 <- optimize_weights(
        inner_v$avg_source_est, variances_k1, inner_v$C_ot, n_samples_k1,
        lambda_reg_k1, inner_v$avg_target_est, inner_v$C_cross,
        clip_weights = FALSE
      )
    }
    fold_weight_optimizer_iterations[k1] <-
      as.integer(attr(eta_k1, "optimizer_iterations"))
    fold_weight_psd_ridge[k1] <- as.numeric(attr(eta_k1, "psd_ridge"))
    if (!is.finite(fold_weight_optimizer_iterations[k1]) ||
        fold_weight_optimizer_iterations[k1] < 1L ||
        !is.finite(fold_weight_psd_ridge[k1]) ||
        fold_weight_psd_ridge[k1] < 0) {
      stop(sprintf(
        ".compute_phase2_weights: invalid optimizer diagnostics for outer fold %d.",
        k1
      ), call. = FALSE)
    }
    fold_weights[k1, ] <- eta_k1

    fold_aggregated_estimates[k1] <- calculate_aggregated_estimate_cpp(
      fold_info[[k1]]$fold_target_estimate,
      fold_info[[k1]]$fold_source_estimates,
      eta_k1
    )

    if (verbose) {
      cat(sprintf("         Fold k1=%d: Weights=[%s]\n",
                  k1, paste(round(eta_k1, 3), collapse = ", ")))
    }
  }

  list(
    fold_weights = fold_weights,
    fold_aggregated_estimates = fold_aggregated_estimates,
    fold_lambdas = fold_lambdas,
    fold_wald_statistics = fold_wald_statistics,
    fold_penalty_coefficients = fold_penalty_coefficients,
    fold_source_included = fold_source_included,
    fold_training_sample_sizes = fold_training_sample_sizes,
    fold_weight_optimizer_iterations = fold_weight_optimizer_iterations,
    fold_weight_psd_ridge = fold_weight_psd_ridge,
    fold_lambda_info = fold_lambda_info
  )
}

.empty_clip_diagnostics <- function() {
  list(
    n_obs = 0L,
    logit_truncated = 0L,
    logit_truncation_fraction = 0,
    weight_min_clipped = 0L,
    weight_max_clipped = 0L,
    ratio_min_clipped = 0L,
    ratio_max_clipped = 0L,
    max_abs_logit = 0,
    max_raw_weight = 0,
    any_truncated = FALSE,
    any_safety_clipped = FALSE,
    any_clipped = FALSE
  )
}

.add_clip_diagnostics <- function(acc, diag) {
  if (is.null(diag)) return(acc)

  for (nm in c("n_obs", "logit_truncated",
               "weight_min_clipped", "weight_max_clipped",
               "ratio_min_clipped", "ratio_max_clipped")) {
    acc[[nm]] <- as.integer(acc[[nm]]) + as.integer(diag[[nm]] %||% 0L)
  }
  # NOTE: we deliberately do NOT pass na.rm = TRUE below. If the C++ clipping
  # diagnostics ever emit NA / NaN for max_abs_logit or max_raw_weight, that is
  # a real signal that something went wrong upstream (e.g., log(0) or Inf
  # propagating into the logit). Letting NA propagate makes the failure visible
  # in the final summary rather than silently absorbing it.
  acc$max_abs_logit  <- max(acc$max_abs_logit,
                            as.numeric(diag$max_abs_logit  %||% 0))
  acc$max_raw_weight <- max(acc$max_raw_weight,
                            as.numeric(diag$max_raw_weight %||% 0))
  acc$logit_truncation_fraction <- if (acc$n_obs > 0L) {
    acc$logit_truncated / acc$n_obs
  } else {
    0
  }
  acc$any_truncated <- acc$logit_truncated > 0L
  acc$any_safety_clipped <-
    (acc$weight_min_clipped + acc$weight_max_clipped +
       acc$ratio_min_clipped + acc$ratio_max_clipped) > 0L
  acc$any_clipped <- acc$any_truncated || acc$any_safety_clipped
  acc
}

.summarize_correction_clipping <- function(fold_results, source_sites) {
  by_source <- setNames(vector("list", length(source_sites)), source_sites)
  total <- .empty_clip_diagnostics()

  for (s in source_sites) {
    source_total <- .empty_clip_diagnostics()
    for (fold_res in fold_results) {
      diag <- fold_res$source_results[[s]]$correction_clip_diagnostics
      source_total <- .add_clip_diagnostics(source_total, diag)
    }
    by_source[[s]] <- source_total
    total <- .add_clip_diagnostics(total, source_total)
  }

  list(total = total, by_source = by_source)
}

.combine_tate_clip_diagnostics <- function(mu1_result, mu0_result) {
  arm_diagnostics <- list(
    treated = mu1_result$clip_diagnostics,
    control = mu0_result$clip_diagnostics
  )
  source_sites <- union(
    names(arm_diagnostics$treated$by_source),
    names(arm_diagnostics$control$by_source)
  )
  by_source <- setNames(vector("list", length(source_sites)), source_sites)
  total <- .empty_clip_diagnostics()
  for (source in source_sites) {
    source_total <- .empty_clip_diagnostics()
    for (arm in arm_diagnostics) {
      source_total <- .add_clip_diagnostics(
        source_total, arm$by_source[[source]]
      )
    }
    by_source[[source]] <- source_total
    total <- .add_clip_diagnostics(total, source_total)
  }
  list(total = total, by_source = by_source, by_arm = arm_diagnostics)
}

.compute_phase3_all_phi <- function(n_folds, fold_weights, fold_info,
                                    n_t, n_source_full, N_all,
                                    source_sites, K, verbose) {
  log_info(verbose, "   Phase 3: Computing per-sample UNCENTERED IF...")

  phi_target <- numeric(n_t)
  phi_source <- lapply(seq_along(source_sites), function(j) numeric(n_source_full[j]))
  target_assignment_count <- integer(n_t)
  source_assignment_count <- lapply(
    n_source_full, function(n_source) integer(n_source)
  )

  for (k1 in 1:n_folds) {
    eta_k1 <- fold_weights[k1, ]
    eta_sum_k1 <- sum(eta_k1)
    info <- fold_info[[k1]]

    phi_ot_raw <- info$varphi_ot + info$fold_target_estimate
    psi_target <- (1 - eta_sum_k1) * phi_ot_raw
    for (j in seq_len(K)) {
      psi_target <- psi_target + eta_k1[j] * (info$zeta_components[[j]] + info$mu_pred_ts[j])
    }
    target_assignment_count <- target_assignment_count +
      .pseudovalue_assignment_counts(
        info$target_idx, psi_target, n_t,
        ".compute_phase3_all_phi target"
      )
    phi_target[info$target_idx] <- (N_all / n_t) * psi_target

    for (j in seq_len(K)) {
      psi_source_j <- eta_k1[j] * (info$xi_components[[j]] + info$delta_ts[j])
      source_assignment_count[[j]] <- source_assignment_count[[j]] +
        .pseudovalue_assignment_counts(
          info$source_idx[[j]], psi_source_j, n_source_full[j],
          sprintf(".compute_phase3_all_phi source %s", source_sites[[j]])
        )
      phi_source[[j]][info$source_idx[[j]]] <- (N_all / n_source_full[j]) * psi_source_j
    }
  }

  if (any(target_assignment_count != 1L) ||
      any(vapply(source_assignment_count, function(counts) {
        any(counts != 1L)
      }, logical(1L)))) {
    stop(
      ".compute_phase3_all_phi: folds must cover every site observation exactly once.",
      call. = FALSE
    )
  }

  all_phi_agg <- c(phi_target, unlist(phi_source))
  stopifnot(length(all_phi_agg) == N_all)
  all_phi_agg
}

#' Compute Cross-fit Aggregation with Sample-Level DML Variance (Version A)
#'
#' Implements the Version A aggregation scheme:
#'
#' \strong{Phase 1:} For each outer fold k, compute per-sample IF components
#'   (varphi_ot, zeta, xi) and variance components (V_ot, V_t, V_s, C_ot,
#'   C_cross) on evaluation fold E_k using out-of-fold nuisances theta^{(-k)}.
#'
#' \strong{Phase 1b:} Inner-fold variance components (Version A Step 2). For each
#'   outer fold k, estimate the variance--covariance components for weight
#'   optimisation via inner M-fold cross-fitting within the training set T_k.
#'   The inner folds coincide with the secondary folds (M = K_f - 1).
#'
#' \strong{Phase 2:} Inner-fold weight estimation (Version A Step 3). For each fold k,
#'   eta^{(-k)} is estimated using variance components from Phase 1b. This
#'   ensures weights are independent of E_k, enabling valid sample-level DML
#'   variance estimation.
#'
#' \strong{Phase 3:} Per-sample UNCENTERED IF. For each observation i in E_k,
#'   compute the uncentered scaled IF Phi_i using fold-specific weights
#'   eta^{(-k)} and out-of-fold nuisances. The point estimate is then
#'   mu_agg = Phi_bar = (1/N) sum Phi_i.
#'
#' \strong{Phase 4:} Sample-level (site-stratified) variance:
#'   Var = (1/N^2) sum_g sum_{i in g} (Phi_i - Phi_bar_g)^2.
#'
#' @param data_split List of site data (from split_data_by_site)
#' @param target_data Target site data list
#' @param source_sites Character vector of source site names
#' @param K Number of source sites
#' @param target_folds List of target fold datasets
#' @param source_folds List of source fold datasets (named by site)
#' @param fold_results List of per-fold estimation results
#' @param n_folds Number of cross-fitting folds
#' @param M_tau Truncation parameter for calibrated losses (training)
#' @param M_tau_inference Truncation radius for inference (default
#'   \code{M_TAU_INFERENCE_DEFAULT} = 5; a single SMMAL-style radius for fitting
#'   and inference, applied to the density-ratio linear predictor in the
#'   correction term and source variance)
#' @param lambda_selection Lambda selection method ("cv" or numeric)
#' @param lambda_rule Nuisance-model CV selection rule: \code{"min"} or
#'   \code{"1se"}.
#' @param aggregation_lambda_grid Optional positive numeric vector of
#'   truncated-Wald multipliers used only when
#'   \code{lambda_selection = "cv"}. Candidate weights are trained and scored
#'   strictly within the outer-training sample.
#' @param verbose Logical; print progress
#' @param final_target_estimate Scalar fold-average target-only estimate,
#'   retained as a diagnostic. The returned target-only point estimate and
#'   variance are computed from observation-level cross-fitted pseudo-values.
#' @param target_estimates Vector of per-fold target estimates
#' @param source_estimates Named numeric vector of source estimates
#' @param source_estimates_matrix Fold-by-source matrix of source estimates
#' @param crossfit_type Character; "two_round" or "one_round"
#' @param algorithm_label Character; algorithm name for the result list
#' @param family_int Integer code for GLM family (0=gaussian, 1=binomial)
#' @param link_int Integer code for link function (0=identity, 1=logit)
#' @param A_val Treatment arm, either 0 or 1.
#' @return List matching the return value of the crossfit algorithms
#' @export
calculate_crossfit_aggregation <- function(data_split, target_data, source_sites, K,
                                         target_folds, source_folds, fold_results,
                                         n_folds, M_tau,
                                         M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                                         lambda_selection, verbose,
                                         lambda_rule = c("min", "1se"),
                                         final_target_estimate, target_estimates,
                                         source_estimates,
                                         source_estimates_matrix,
                                         crossfit_type, algorithm_label,
                                         family_int = 1L, link_int = 1L,
                                         A_val = 1L,
                                         aggregation_lambda_grid = NULL) {
  lambda_rule <- match.arg(lambda_rule)
  n_t <- target_data$n
  source_estimates_matrix <- as.matrix(source_estimates_matrix)
  colnames(source_estimates_matrix) <- source_sites

  # Total sample size across all sites
  n_source_full <- vapply(source_sites, function(s) data_split[[s]]$n, numeric(1L))
  N_all <- n_t + sum(n_source_full)

  # ===========================================================================
  # PHASE 1: Per-fold IF vectors on evaluation folds
  # ===========================================================================
  # For each outer fold k1, compute IF components (varphi_ot, zeta, xi) on
  # E_k1 using averaged nuisance estimates theta^(-k1) trained on T_k1.
  # These IF components are used in Phase 3 for per-sample Phi_i.
  #
  # Note: V/C variance components computed here are diagnostic only.
  # Weight estimation (Phase 2) uses inner-fold variance from Phase 1b.
  # ===========================================================================
  fold_info <- vector("list", n_folds)
  for (k1 in 1:n_folds) {
    log_info(verbose, "      Phase 1: IF components for fold k1=%d/%d", k1, n_folds)
    fold_info[[k1]] <- .compute_phase1_fold_info(
      k1 = k1,
      source_sites = source_sites,
      K = K,
      fold_results = fold_results,
      target_folds = target_folds,
      source_folds = source_folds,
      M_tau_inference = M_tau_inference,
      family_int = family_int,
      link_int = link_int,
      A_val = A_val
    )
  }

  # ===========================================================================
  # PHASE 1b: Inner-fold variance components (Version A Step 2)
  # ===========================================================================
  # For each outer fold k1, estimate the variance--covariance components
  # needed for weight optimisation eta^{(-k1)} via an inner M-fold
  # cross-fitting loop within the training set T_{k1} = D \ D_{k1}.
  #
  # The inner folds coincide with the secondary folds used for nuisance
  # calibration: M = K_f - 1, with inner folds {k2 : k2 != k1}.
  # For each inner fold k2:
  #   - Use per-(k1,k2) nuisances theta^{(-k1,-k2)} (trained on
  #     data excluding both k1 and k2) to compute IF components on fold k2.
  #   - Use inner-fold target-only estimates theta_ot^{(-k1,-k2)} for
  #     varphi_ot on fold k2.
  #   - Compute all variance components (V_ot, V_t, V_s, C_ot, C_cross)
  #     on the inner evaluation fold D_{k2}.
  # The components are averaged across the K_f - 1 inner folds.
  #
  # Because the entire computation uses only data from T_{k1} (with
  # out-of-inner-fold nuisances), eta^{(-k1)} is independent of E_{k1}.
  # ===========================================================================
  inner_fold_info <- vector("list", n_folds)
  for (k1 in 1:n_folds) {
    log_info(verbose, "      Phase 1b: Inner-fold variance for fold k1=%d/%d", k1, n_folds)
    inner_fold_info[[k1]] <- .compute_phase1b_inner_fold_info(
      k1 = k1,
      n_folds = n_folds,
      source_sites = source_sites,
      K = K,
      fold_results = fold_results,
      target_folds = target_folds,
      source_folds = source_folds,
      M_tau_inference = M_tau_inference,
      family_int = family_int,
      link_int = link_int,
      A_val = A_val
    )
  }

  # ===========================================================================
  # PHASE 2: Inner-fold weight estimation (Version A Step 3)
  # ===========================================================================
  # For each fold k1, eta^{(-k1)} is estimated using variance components
  # from the inner M-fold cross-fitting (Phase 1b). This makes eta^{(-k1)}
  # independent of E_k1 data, enabling valid per-sample IF computation
  # on E_k1 in Phase 3.
  #
  # For the weight optimization objective, we use:
  #   - Variance components (V_ot, V_t, V_s, C_ot, C_cross): inner-fold avg
  #   - Point estimates (mu_ot, mu_ts) for penalty: inner-fold avg
  #   - Sample sizes: sizes in T_{k1}, because both the source--target
  #     discrepancy and its standard error are learned inside T_{k1}. With
  #     balanced folds, multiplying every site size by the same training
  #     fraction leaves N * Var(eta) unchanged, while correctly studentizing
  #     the training-fold discrepancy.
  # ===========================================================================
  phase2_res <- .compute_phase2_weights(
    n_folds = n_folds,
    inner_fold_info = inner_fold_info,
    fold_info = fold_info,
    lambda_selection = lambda_selection,
    K = K,
    verbose = verbose,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid
  )
  fold_weights <- phase2_res$fold_weights
  fold_aggregated_estimates <- phase2_res$fold_aggregated_estimates
  fold_lambdas <- phase2_res$fold_lambdas
  fold_wald_statistics <- phase2_res$fold_wald_statistics
  fold_penalty_coefficients <- phase2_res$fold_penalty_coefficients
  fold_training_sample_sizes <- phase2_res$fold_training_sample_sizes
  fold_weight_optimizer_iterations <-
    phase2_res$fold_weight_optimizer_iterations
  fold_weight_psd_ridge <- phase2_res$fold_weight_psd_ridge
  fold_lambda_info <- phase2_res$fold_lambda_info

  colnames(fold_weights) <- source_sites
  colnames(fold_wald_statistics) <- source_sites
  colnames(fold_penalty_coefficients) <- source_sites

  average_weights <- colMeans(fold_weights)

  if (verbose) {
    cat("   FOLD-SPECIFIC WEIGHTS (inner-fold estimation):\n")
    for (k1 in 1:n_folds) {
      cat(sprintf("      Fold k1=%d: [%s]\n", k1,
                  paste(round(fold_weights[k1, ], 3), collapse = ", ")))
    }
    cat(sprintf("   Average weights: [%s]\n",
                paste(round(average_weights, 4), collapse = ", ")))
  }

  # ===========================================================================
  # PHASE 3: Per-sample UNCENTERED influence function
  # ===========================================================================
  # For each observation i in evaluation fold E_{k(i)}, compute the UNCENTERED
  # scaled IF so that the point estimate and variance share the same definition:
  #
  #   mu_agg = (1/N) * sum_i Phi_i  =  Phi_bar
  #   Var    = (1/N^2) * sum_g sum_{i in g} (Phi_i - Phi_bar_g)^2
  #
  #   Target obs i:  Phi_i = (N / n_t) * psi_i^{target}
  #     where psi_i^{target} = (1 - sum_j eta_j) * phi_ot_i
  #                          + sum_j eta_j * psi_j(X_i)
  #     phi_ot_i = AIPW pseudo-outcome (uncentered)
  #     psi_j(X_i) = outcome prediction from source model j (uncentered)
  #
  #   Source obs i (site j):  Phi_i = (N / n_{s_j}) * psi_i^{source}
  #     where psi_i^{source} = eta_j * w_i(Y_i - psi_j(X_i))  (uncentered)
  #
  # This ensures mu_agg = Phi_bar = (1/N) sum Phi_i exactly.
  # Under the multi-sample (site-stratified) regime, the corresponding
  # variance uses within-site centering.
  # ===========================================================================
  all_phi_agg <- .compute_phase3_all_phi(
    n_folds = n_folds,
    fold_weights = fold_weights,
    fold_info = fold_info,
    n_t = n_t,
    n_source_full = n_source_full,
    N_all = N_all,
    source_sites = source_sites,
    K = K,
    verbose = verbose
  )

  # ===========================================================================
  # PHASE 4: Final aggregation with sample-level DML variance
  # ===========================================================================
  log_info(verbose, "   Phase 4: Sample-level DML variance...")

  final_agg <- aggregate_fold_estimates(
    fold_aggregated_estimates, all_phi_agg,
    N_all, n_t, n_source_full,
    fold_weights, n_folds, verbose,
    fold_info = fold_info,
    inner_fold_info = inner_fold_info,
    fold_lambdas = fold_lambdas,
    fold_weight_psd_ridge = fold_weight_psd_ridge,
    screening_rule = "soft_penalty"
  )

  final_estimate <- final_agg$final_estimate
  final_variance <- final_agg$final_variance
  se <- final_agg$se
  ci_lower <- final_agg$ci_lower
  ci_upper <- final_agg$ci_upper
  weights <- final_agg$average_weights
  names(weights) <- source_sites
  clip_diagnostics <- .summarize_correction_clipping(fold_results, source_sites)
  if (isTRUE(clip_diagnostics$total$any_safety_clipped)) {
    warning(sprintf(
      "Inference clipping was triggered: weight_min=%d, weight_max=%d, ratio_min=%d, ratio_max=%d. Inspect result$clip_diagnostics.",
      clip_diagnostics$total$weight_min_clipped,
      clip_diagnostics$total$weight_max_clipped,
      clip_diagnostics$total$ratio_min_clipped,
      clip_diagnostics$total$ratio_max_clipped
    ))
  }
  phase1_summary <- list(
    V_ot = vapply(fold_info, function(x) x$V_ot, numeric(1L)),
    V_t = do.call(rbind, lapply(fold_info, function(x) x$V_t)),
    V_s = do.call(rbind, lapply(fold_info, function(x) x$V_s)),
    C_ot = do.call(rbind, lapply(fold_info, function(x) x$C_ot)),
    fold_target_estimate = vapply(fold_info, function(x) x$fold_target_estimate, numeric(1L)),
    fold_source_estimates = do.call(rbind, lapply(fold_info, function(x) x$fold_source_estimates)),
    mu_pred_ts = do.call(rbind, lapply(fold_info, function(x) x$mu_pred_ts)),
    delta_ts = do.call(rbind, lapply(fold_info, function(x) x$delta_ts))
  )
  phase1b_summary <- list(
    V_ot = vapply(inner_fold_info, function(x) x$V_ot, numeric(1L)),
    V_t = do.call(rbind, lapply(inner_fold_info, function(x) x$V_t)),
    V_s = do.call(rbind, lapply(inner_fold_info, function(x) x$V_s)),
    C_ot = do.call(rbind, lapply(inner_fold_info, function(x) x$C_ot)),
    avg_target_est = vapply(inner_fold_info, function(x) x$avg_target_est, numeric(1L)),
    avg_source_est = do.call(rbind, lapply(inner_fold_info, function(x) x$avg_source_est)),
    mu_pred_ts = do.call(rbind, lapply(inner_fold_info, function(x) x$mu_pred_ts)),
    delta_ts = do.call(rbind, lapply(inner_fold_info, function(x) x$delta_ts)),
    n_t = vapply(inner_fold_info, function(x) x$n_t, numeric(1L)),
    n_s = do.call(rbind, lapply(inner_fold_info, function(x) x$n_s))
  )
  colnames(phase1b_summary$avg_source_est) <- source_sites
  target_only_phi <- .assemble_target_pseudovalues(
    fold_info, n_t, "calculate_crossfit_aggregation"
  )
  target_only_estimate <- mean(target_only_phi)
  target_only_variance <- .multisite_pseudovalue_variance(
    target_only_phi, as.integer(n_t)
  )

  results <- list(
    estimate = final_estimate,
    se = se,
    variance = final_variance,
    se_fixed_weights = final_agg$se_fixed_weights,
    variance_fixed_weights = final_agg$variance_fixed_weights,
    weight_layer = final_agg$weight_layer,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    target_only = list(
      estimate = target_only_estimate,
      se = sqrt(target_only_variance),
      variance = target_only_variance,
      fold_average_estimate = final_target_estimate
    ),
    source_estimates = source_estimates,
    weights = weights,
    fold_weights = fold_weights,
    fold_lambdas = fold_lambdas,
    fold_wald_statistics = fold_wald_statistics,
    fold_penalty_coefficients = fold_penalty_coefficients,
    fold_training_sample_sizes = fold_training_sample_sizes,
    fold_weight_optimizer_iterations = fold_weight_optimizer_iterations,
    fold_weight_psd_ridge = fold_weight_psd_ridge,
    fold_lambda_info = fold_lambda_info,
    fold_aggregated_estimates = fold_aggregated_estimates,
    aggregation_lambda_selection = lambda_selection,
    aggregation_lambda_grid = aggregation_lambda_grid,
    aggregation_lambda_rule = lambda_rule,
    clip_diagnostics = clip_diagnostics,
    n_sites = K + 1,
    n_folds = n_folds,
    N_all = N_all,
    all_phi_agg = all_phi_agg,
    intermediates = list(
      target_estimates = target_estimates,
      target_only_phi = target_only_phi,
      source_estimates_matrix = source_estimates_matrix,
      fold_info = fold_info,
      inner_fold_info = inner_fold_info,
      phase1 = phase1_summary,
      phase1b = phase1b_summary,
      fold_lambdas = fold_lambdas,
      fold_wald_statistics = fold_wald_statistics,
      fold_penalty_coefficients = fold_penalty_coefficients,
      fold_training_sample_sizes = fold_training_sample_sizes,
      fold_weight_optimizer_iterations = fold_weight_optimizer_iterations,
      fold_weight_psd_ridge = fold_weight_psd_ridge,
      fold_lambda_info = fold_lambda_info,
      fold_weights = fold_weights,
      fold_aggregated_estimates = fold_aggregated_estimates,
      aggregation_lambda_rule = lambda_rule,
      clip_diagnostics = clip_diagnostics,
      sample_sizes = list(n_t = n_t, n_source = n_source_full, N_all = N_all)
    ),
    method = algorithm_label,    # consistent with comparison methods
    fold_results = fold_results
  )

  if (verbose) {
    cat(paste0(paste(rep("=", 70), collapse = ""), "\n"))
    cat(sprintf("==> %s COMPLETED!\n", toupper(gsub("_", " ", algorithm_label))))
    cat(paste0(paste(rep("=", 70), collapse = ""), "\n\n"))
  }

  return(results)
}

#' TATE aggregation from arm-specific cross-fitted nuisance fits
#'
#' Combines treated- and control-arm cross-fitting results before aggregation.
#' A single source-weight vector is learned by minimizing the estimated
#' variance of the TATE contrast, so all treated/control covariance terms enter
#' both the weight objective and the final Wald variance.
#'
#' @param data_split Site-stratified data list.
#' @param mu1_result Result returned by \code{run_crossfit(..., A_val = 1)}.
#' @param mu0_result Result returned by \code{run_crossfit(..., A_val = 0)}.
#' @param lambda_selection Aggregation Wald-penalty factor. The default
#'   \code{AGG_WALD_LAMBDA = 1} is the pilot-selected cutoff \eqn{c=1}
#'   specified in \code{main.tex}; \code{"cv"} remains available for
#'   sensitivity analyses.
#' @param lambda_rule Aggregation inner-validation rule.
#' @param aggregation_lambda_grid Optional positive numeric vector of
#'   truncated-Wald multipliers used only when
#'   \code{lambda_selection = "cv"}.
#' @param screening_rule Source-screening rule. \code{"soft_penalty"} uses
#'   the manuscript's truncated-Wald L1 penalty. \code{"hard_threshold"} is
#'   an explicit diagnostic that excludes sources above the foldwise Wald
#'   cutoff before variance minimization. \code{"quadratic_bias"} is the
#'   pre-specified smooth sensitivity rule with penalty
#'   \eqn{n_t^{3/4}\sum_j \delta_j^2\eta_j^2} (main.tex
#'   rem:quadratic_bias_rule); it never sets a weight exactly to zero.
#' @param verbose Print progress.
#' @return A RoCE result list for the TATE estimator.
#' @export
calculate_tate_crossfit_aggregation <- function(
    data_split, mu1_result, mu0_result,
    lambda_selection = AGG_WALD_LAMBDA,
    lambda_rule = c("min", "1se"),
    verbose = TRUE,
    aggregation_lambda_grid = NULL,
    screening_rule = c("soft_penalty", "hard_threshold", "quadratic_bias")) {
  lambda_rule <- match.arg(lambda_rule)
  screening_rule <- match.arg(screening_rule)
  required_intermediates <- c("fold_info", "inner_fold_info")
  for (arm_result in list(mu1_result, mu0_result)) {
    missing <- required_intermediates[!vapply(required_intermediates, function(name) {
      !is.null(arm_result$intermediates[[name]])
    }, logical(1L))]
    if (length(missing) > 0L) {
      stop(sprintf(
        "calculate_tate_crossfit_aggregation: arm result is missing aggregation intermediate(s): %s.",
        paste(missing, collapse = ", ")
      ))
    }
  }

  n_folds <- mu1_result$n_folds
  if (!identical(as.integer(n_folds), as.integer(mu0_result$n_folds)) ||
      length(mu1_result$intermediates$fold_info) !=
        length(mu0_result$intermediates$fold_info)) {
    stop("calculate_tate_crossfit_aggregation: treated/control results use different fold counts.")
  }
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  K <- length(source_sites)
  n_t <- target_data$n
  n_source_full <- vapply(source_sites, function(s) data_split[[s]]$n, numeric(1L))
  N_all <- n_t + sum(n_source_full)

  fold_info <- vector("list", n_folds)
  inner_fold_info <- vector("list", n_folds)
  for (k1 in seq_len(n_folds)) {
    log_info(verbose, "      TATE Phase 1: fold k1=%d/%d", k1, n_folds)
    fold_info[[k1]] <- .combine_tate_fold_info(
      mu1_result$intermediates$fold_info[[k1]],
      mu0_result$intermediates$fold_info[[k1]],
      K
    )

    log_info(verbose, "      TATE Phase 1b: fold k1=%d/%d", k1, n_folds)
    inner_fold_info[[k1]] <- .combine_tate_inner_fold_info(
      mu1_result$intermediates$inner_fold_info[[k1]],
      mu0_result$intermediates$inner_fold_info[[k1]],
      K
    )
  }

  phase2 <- .compute_phase2_weights(
    n_folds = n_folds,
    inner_fold_info = inner_fold_info,
    fold_info = fold_info,
    lambda_selection = lambda_selection,
    K = K,
    verbose = verbose,
    lambda_rule = lambda_rule,
    aggregation_lambda_grid = aggregation_lambda_grid,
    screening_rule = screening_rule
  )
  colnames(phase2$fold_weights) <- source_sites
  colnames(phase2$fold_wald_statistics) <- source_sites
  colnames(phase2$fold_penalty_coefficients) <- source_sites
  colnames(phase2$fold_source_included) <- source_sites

  all_phi_tau <- .compute_phase3_all_phi(
    n_folds = n_folds,
    fold_weights = phase2$fold_weights,
    fold_info = fold_info,
    n_t = n_t,
    n_source_full = n_source_full,
    N_all = N_all,
    source_sites = source_sites,
    K = K,
    verbose = verbose
  )
  final <- aggregate_fold_estimates(
    phase2$fold_aggregated_estimates, all_phi_tau,
    N_all, n_t, n_source_full,
    phase2$fold_weights, n_folds, verbose,
    fold_info = fold_info,
    inner_fold_info = inner_fold_info,
    fold_lambdas = phase2$fold_lambdas,
    fold_weight_psd_ridge = phase2$fold_weight_psd_ridge,
    screening_rule = screening_rule
  )

  target_tau_raw <- .assemble_target_pseudovalues(
    fold_info, n_t, "calculate_tate_crossfit_aggregation"
  )
  target_tau_est <- mean(target_tau_raw)
  target_tau_var <- .multisite_pseudovalue_variance(
    target_tau_raw, as.integer(n_t)
  )
  target_estimates <- vapply(
    fold_info, function(x) x$fold_target_estimate, numeric(1L)
  )
  source_estimates_matrix <- do.call(
    rbind, lapply(fold_info, function(x) x$fold_source_estimates)
  )
  colnames(source_estimates_matrix) <- source_sites
  target_fold_sizes <- vapply(
    fold_info, function(info) length(info$target_idx), numeric(1L)
  )
  source_fold_sizes <- do.call(rbind, lapply(fold_info, function(info) {
    vapply(info$source_idx, length, numeric(1L))
  }))
  source_estimates <- stats::setNames(
    .pool_fold_pairwise_estimates(
      do.call(rbind, lapply(fold_info, `[[`, "mu_pred_ts")),
      do.call(rbind, lapply(fold_info, `[[`, "delta_ts")),
      target_fold_sizes,
      source_fold_sizes
    ),
    source_sites
  )
  clip_diagnostics <- .combine_tate_clip_diagnostics(
    mu1_result, mu0_result
  )

  phase1_summary <- list(
    V_ot = vapply(fold_info, `[[`, numeric(1L), "V_ot"),
    V_t = do.call(rbind, lapply(fold_info, `[[`, "V_t")),
    V_s = do.call(rbind, lapply(fold_info, `[[`, "V_s")),
    C_ot = do.call(rbind, lapply(fold_info, `[[`, "C_ot")),
    C_cross = lapply(fold_info, `[[`, "C_cross"),
    fold_target_estimate = target_estimates,
    fold_source_estimates = source_estimates_matrix,
    mu_pred_ts = do.call(rbind, lapply(fold_info, `[[`, "mu_pred_ts")),
    delta_ts = do.call(rbind, lapply(fold_info, `[[`, "delta_ts"))
  )
  phase1b_summary <- list(
    V_ot = vapply(inner_fold_info, `[[`, numeric(1L), "V_ot"),
    V_t = do.call(rbind, lapply(inner_fold_info, `[[`, "V_t")),
    V_s = do.call(rbind, lapply(inner_fold_info, `[[`, "V_s")),
    C_ot = do.call(rbind, lapply(inner_fold_info, `[[`, "C_ot")),
    C_cross = lapply(inner_fold_info, `[[`, "C_cross"),
    avg_target_est = vapply(inner_fold_info, `[[`, numeric(1L), "avg_target_est"),
    avg_source_est = do.call(rbind, lapply(inner_fold_info, `[[`, "avg_source_est")),
    mu_pred_ts = do.call(rbind, lapply(inner_fold_info, `[[`, "mu_pred_ts")),
    delta_ts = do.call(rbind, lapply(inner_fold_info, `[[`, "delta_ts")),
    n_t = vapply(inner_fold_info, `[[`, numeric(1L), "n_t"),
    n_s = do.call(rbind, lapply(inner_fold_info, `[[`, "n_s"))
  )
  colnames(phase1b_summary$avg_source_est) <- source_sites

  list(
    estimate = final$final_estimate,
    se = final$se,
    variance = final$final_variance,
    se_fixed_weights = final$se_fixed_weights,
    variance_fixed_weights = final$variance_fixed_weights,
    weight_layer = final$weight_layer,
    ci_lower = final$ci_lower,
    ci_upper = final$ci_upper,
    target_only = list(
      estimate = target_tau_est,
      se = sqrt(target_tau_var),
      variance = target_tau_var
    ),
    source_estimates = source_estimates,
    weights = stats::setNames(final$average_weights, source_sites),
    fold_weights = phase2$fold_weights,
    fold_lambdas = phase2$fold_lambdas,
    fold_wald_statistics = phase2$fold_wald_statistics,
    fold_penalty_coefficients = phase2$fold_penalty_coefficients,
    fold_source_included = phase2$fold_source_included,
    fold_training_sample_sizes = phase2$fold_training_sample_sizes,
    fold_weight_optimizer_iterations = phase2$fold_weight_optimizer_iterations,
    fold_weight_psd_ridge = phase2$fold_weight_psd_ridge,
    fold_lambda_info = phase2$fold_lambda_info,
    fold_aggregated_estimates = phase2$fold_aggregated_estimates,
    aggregation_lambda_selection = lambda_selection,
    aggregation_lambda_grid = aggregation_lambda_grid,
    aggregation_lambda_rule = lambda_rule,
    aggregation_screening_rule = screening_rule,
    clip_diagnostics = clip_diagnostics,
    n_sites = K + 1L,
    n_folds = n_folds,
    N_all = N_all,
    all_phi_agg = all_phi_tau,
    all_phi_tau = all_phi_tau,
    estimand = "TATE",
    method = "direct_tate_two_layer_crossfit",
    arm_results = list(mu1 = mu1_result, mu0 = mu0_result),
    intermediates = list(
      target_estimates = target_estimates,
      target_only_phi = target_tau_raw,
      source_estimates_matrix = source_estimates_matrix,
      fold_info = fold_info,
      inner_fold_info = inner_fold_info,
      phase1 = phase1_summary,
      phase1b = phase1b_summary,
      fold_lambdas = phase2$fold_lambdas,
      fold_wald_statistics = phase2$fold_wald_statistics,
      fold_penalty_coefficients = phase2$fold_penalty_coefficients,
      fold_source_included = phase2$fold_source_included,
      fold_training_sample_sizes = phase2$fold_training_sample_sizes,
      fold_weight_optimizer_iterations = phase2$fold_weight_optimizer_iterations,
      fold_weight_psd_ridge = phase2$fold_weight_psd_ridge,
      fold_lambda_info = phase2$fold_lambda_info,
      fold_weights = phase2$fold_weights,
      fold_aggregated_estimates = phase2$fold_aggregated_estimates,
      clip_diagnostics = clip_diagnostics,
      sample_sizes = list(n_t = n_t, n_source = n_source_full, N_all = N_all)
    )
  )
}
