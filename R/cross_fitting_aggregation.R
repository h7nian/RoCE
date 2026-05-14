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
#
#' Aggregate fold-level estimates into a single cross-fitted estimate with
#' sample-level (DML-style) variance
#'
#' This is the final step of Version A aggregation.  It receives UNCENTERED
#' per-sample IF values from all outer folds (each observation appears exactly
#' once) and computes both the point estimate and variance from the same
#' set of values:
#'   mu_agg = Phi_bar = (1/N) sum Phi_i
#'   Var    = (1/N^2) sum_g sum_{i in group g} (Phi_i - Phi_bar_g)^2
#' where groups are the target and each source site.
#'
#' Under the multi-sample (site-stratified) asymptotic regime (fixed sites,
#' independent samples per site with fixed sampling fractions), the variance of
#' the overall mean depends on within-site variation only; global centering would
#' add an ANOVA between-site term that does not correspond to sampling
#' variability.
#'
#' This keeps a unified estimating-equation workflow: both point estimate and
#' variance are computed from the same set {Phi_i}.
#'
#' @param fold_aggregated_estimates Vector of fold-specific aggregated estimates
#' @param all_phi_agg Numeric vector of per-sample aggregated IF values
#'   (length = N_all, each observation appearing exactly once across all folds)
#' @param N_all Total sample size across all sites
#' @param n_t Target site sample size
#' @param n_source_full Numeric vector of source-site sample sizes (length K)
#' @param fold_weights Matrix of fold-specific weights (n_folds x K)
#' @param n_folds Number of cross-fitting folds
#' @param verbose Print progress
#' @return List with final_estimate, final_variance, se, ci_lower, ci_upper, average_weights
aggregate_fold_estimates <- function(fold_aggregated_estimates, all_phi_agg,
                                         N_all, n_t, n_source_full,
                                         fold_weights, n_folds, verbose) {
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

  wss <- 0
  start <- 1
  for (g_size in group_sizes) {
    if (g_size > 0) {
      end <- start + g_size - 1
      phi_g <- all_phi_agg[start:end]
      phi_bar_g <- mean(phi_g)
      wss <- wss + stable_sum_kahan((phi_g - phi_bar_g)^2)
      start <- end + 1
    }
  }
  final_variance <- wss / (N_all^2)

  # Ensure non-negative variance
  final_variance <- max(final_variance, 0)

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
    cat(sprintf("      - Variance (within-site): %.6f (SE=%.4f)\n", final_variance, se))
    cat(paste0("   95% CI: [", round(ci_lower, 4), ", ", round(ci_upper, 4), "]\n"))
  }

  return(list(
    final_estimate = final_estimate,
    final_variance = final_variance,
    se = se,
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

.compute_phase1b_inner_fold_info <- function(k1, n_folds, source_sites, K,
                                             fold_results, target_folds,
                                             source_folds, M_tau_inference,
                                             family_int, link_int, A_val) {
  secondary_folds <- setdiff(1:n_folds, k1)
  n_inner <- length(secondary_folds)

  V_ot_sum <- 0
  V_t_sum <- numeric(K)
  V_s_sum <- numeric(K)
  C_ot_sum <- numeric(K)
  C_cross_sum <- matrix(0, nrow = K, ncol = K)
  target_est_sum <- 0
  source_est_sum <- numeric(K)

  for (k2 in secondary_folds) {
    k2_key <- paste0("k2_", k2)

    to_inner <- fold_results[[k1]]$target_only_inner[[k2_key]]
    varphi_ot_inner <- to_inner$varphi_ot
    V_ot_sum <- V_ot_sum + to_inner$V_ot
    target_est_sum <- target_est_sum + to_inner$estimate

    target_fold_k2 <- materialize_fold(target_folds, k2)
    zeta_list_inner <- vector("list", K)

    for (i in seq_along(source_sites)) {
      s <- source_sites[i]
      src <- fold_results[[k1]]$source_results[[s]]
      gamma_k2 <- src$per_k2_gamma[[k2_key]]
      alpha_k2 <- src$per_k2_alpha[[k2_key]]

      if (is.null(gamma_k2) || is.null(alpha_k2)) next

      source_fold_k2 <- materialize_fold(source_folds[[s]], k2)

      correction_inner <- calculate_correction_term_cpp(
        source_fold_k2$Z_site, source_fold_k2$A, source_fold_k2$Y,
        gamma_k2, alpha_k2, source_fold_k2$W_outcome, M_tau_inference,
        family_int, link_int, A_val
      )
      delta_inner <- correction_inner$delta_ts

      mu_pred_inner <- mean(predict_glm_cpp(
        target_fold_k2$W_outcome, alpha_k2, family_int, link_int
      ))

      source_est_sum[i] <- source_est_sum[i] + mu_pred_inner + delta_inner

      sv_inner <- calculate_source_variance_cpp(
        source_fold_k2$Z_site, source_fold_k2$A, source_fold_k2$Y,
        gamma_k2, alpha_k2, delta_inner,
        source_fold_k2$W_outcome, M_tau_inference,
        family_int, link_int, A_val
      )
      V_s_sum[i] <- V_s_sum[i] + sv_inner$V_s

      tv_inner <- calculate_target_variance_cpp(
        target_fold_k2$W_outcome, alpha_k2, mu_pred_inner,
        family_int, link_int
      )
      V_t_sum[i] <- V_t_sum[i] + tv_inner$V_t
      zeta_list_inner[[i]] <- tv_inner$zeta_components

      C_ot_sum[i] <- C_ot_sum[i] + calculate_covariance_term_cpp(
        varphi_ot_inner, tv_inner$zeta_components
      )
    }

    C_cross_sum <- C_cross_sum + .compute_cross_covariance_matrix(zeta_list_inner, K)
  }

  list(
    V_ot = V_ot_sum / n_inner,
    V_t = V_t_sum / n_inner,
    V_s = V_s_sum / n_inner,
    C_ot = C_ot_sum / n_inner,
    C_cross = C_cross_sum / n_inner,
    avg_target_est = target_est_sum / n_inner,
    avg_source_est = source_est_sum / n_inner
  )
}

.compute_phase2_weights <- function(n_folds, inner_fold_info, fold_info,
                                    n_t, n_source_full,
                                    lambda_selection, crossfit_type, K,
                                    verbose, lambda_rule = c("min", "1se")) {
  lambda_rule <- match.arg(lambda_rule)
  fold_aggregated_estimates <- numeric(n_folds)
  fold_weights <- array(0, dim = c(n_folds, K))
  fold_lambdas <- numeric(n_folds)

  for (k1 in 1:n_folds) {
    log_info(verbose, "      Phase 2: Inner-fold weights for fold k1=%d/%d", k1, n_folds)

    inner_v <- inner_fold_info[[k1]]
    variances_k1 <- list(V_ot = inner_v$V_ot, V_t = inner_v$V_t, V_s = inner_v$V_s)
    n_samples_k1 <- list(n_t = as.numeric(n_t), n_s = as.numeric(n_source_full))

    lambda_reg_k1 <- if (lambda_selection == "cv") {
      select_lambda_cv_crossfit(
        n_folds, inner_v$avg_target_est, inner_v$avg_source_est,
        variances_k1, inner_v$C_ot, n_samples_k1, inner_v$C_cross, verbose = FALSE,
        crossfit_type = crossfit_type,
        lambda_rule = lambda_rule
      )
    } else if (is.numeric(lambda_selection)) {
      lambda_selection
    } else {
      LAMBDA_DEFAULT
    }
    fold_lambdas[k1] <- lambda_reg_k1

    eta_k1 <- optimize_weights(
      inner_v$avg_source_est, variances_k1, inner_v$C_ot, n_samples_k1,
      lambda_reg_k1, inner_v$avg_target_est, inner_v$C_cross, clip_weights = FALSE
    )
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
    fold_lambdas = fold_lambdas
  )
}

.empty_clip_diagnostics <- function() {
  list(
    n_obs = 0L,
    weight_min_clipped = 0L,
    weight_max_clipped = 0L,
    ratio_min_clipped = 0L,
    ratio_max_clipped = 0L,
    max_abs_logit = 0,
    max_raw_weight = 0,
    any_clipped = FALSE
  )
}

.add_clip_diagnostics <- function(acc, diag) {
  if (is.null(diag)) return(acc)

  for (nm in c("n_obs", "weight_min_clipped", "weight_max_clipped",
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
  acc$any_clipped <- (acc$weight_min_clipped + acc$weight_max_clipped +
                        acc$ratio_min_clipped + acc$ratio_max_clipped) > 0L
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

.compute_phase3_all_phi <- function(n_folds, fold_weights, fold_info,
                                    n_t, n_source_full, N_all,
                                    source_sites, K, verbose) {
  log_info(verbose, "   Phase 3: Computing per-sample UNCENTERED IF...")

  phi_target <- numeric(n_t)
  phi_source <- lapply(seq_along(source_sites), function(j) numeric(n_source_full[j]))

  for (k1 in 1:n_folds) {
    eta_k1 <- fold_weights[k1, ]
    eta_sum_k1 <- sum(eta_k1)
    info <- fold_info[[k1]]

    phi_ot_raw <- info$varphi_ot + info$fold_target_estimate
    psi_target <- (1 - eta_sum_k1) * phi_ot_raw
    for (j in seq_len(K)) {
      psi_target <- psi_target + eta_k1[j] * (info$zeta_components[[j]] + info$mu_pred_ts[j])
    }
    phi_target[info$target_idx] <- (N_all / n_t) * psi_target

    for (j in seq_len(K)) {
      psi_source_j <- eta_k1[j] * (info$xi_components[[j]] + info$delta_ts[j])
      phi_source[[j]][info$source_idx[[j]]] <- (N_all / n_source_full[j]) * psi_source_j
    }
  }

  all_phi_agg <- c(phi_target, unlist(phi_source))
  stopifnot(length(all_phi_agg) == N_all)
  all_phi_agg
}

#' Compute Cross-fit Aggregation with Sample-Level DML Variance (Version A)
#'
#' Implements the Version A aggregation scheme:
#'
#' **Phase 1:** For each outer fold k, compute per-sample IF components
#'   (varphi_ot, zeta, xi) and variance components (V_ot, V_t, V_s, C_ot,
#'   C_cross) on evaluation fold E_k using out-of-fold nuisances theta^{(-k)}.
#'
#' **Phase 1b:** Inner-fold variance components (Version A Step 2). For each
#'   outer fold k, estimate the variance--covariance components for weight
#'   optimisation via inner M-fold cross-fitting within the training set T_k.
#'   The inner folds coincide with the secondary folds (M = K_f - 1).
#'
#' **Phase 2:** Inner-fold weight estimation (Version A Step 3). For each fold k,
#'   eta^{(-k)} is estimated using variance components from Phase 1b. This
#'   ensures weights are independent of E_k, enabling valid sample-level DML
#'   variance estimation.
#'
#' **Phase 3:** Per-sample UNCENTERED IF. For each observation i in E_k,
#'   compute the uncentered scaled IF Phi_i using fold-specific weights
#'   eta^{(-k)} and out-of-fold nuisances. The point estimate is then
#'   mu_agg = Phi_bar = (1/N) sum Phi_i.
#'
#' **Phase 4:** Sample-level (site-stratified) variance:
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
#' @param M_tau_inference Truncation parameter for inference (default Inf = no truncation)
#' @param lambda_selection Lambda selection method ("cv" or numeric)
#' @param verbose Logical; print progress
#' @param final_target_estimate Scalar; cross-fitted target-only estimate
#' @param target_estimates Vector of per-fold target estimates
#' @param source_estimates Named numeric vector of source estimates
#' @param source_estimates_matrix Fold-by-source matrix of source estimates
#' @param crossfit_type Character; "two_round" or "one_round"
#' @param algorithm_label Character; algorithm name for the result list
#' @param family_int Integer code for GLM family (0=gaussian, 1=binomial)
#' @param link_int Integer code for link function (0=identity, 1=logit)
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
                                         A_val = 1L) {
  lambda_rule <- match.arg(lambda_rule)
  n_t <- target_data$n

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
  #   - Sample sizes: FULL sizes (n_t, n_sj) for overall variance optimization
  # ===========================================================================
  phase2_res <- .compute_phase2_weights(
    n_folds = n_folds,
    inner_fold_info = inner_fold_info,
    fold_info = fold_info,
    n_t = n_t,
    n_source_full = n_source_full,
    lambda_selection = lambda_selection,
    crossfit_type = crossfit_type,
    K = K,
    verbose = verbose,
    lambda_rule = lambda_rule
  )
  fold_weights <- phase2_res$fold_weights
  fold_aggregated_estimates <- phase2_res$fold_aggregated_estimates
  fold_lambdas <- phase2_res$fold_lambdas

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
    fold_weights, n_folds, verbose
  )

  final_estimate <- final_agg$final_estimate
  final_variance <- final_agg$final_variance
  se <- final_agg$se
  ci_lower <- final_agg$ci_lower
  ci_upper <- final_agg$ci_upper
  weights <- final_agg$average_weights
  clip_diagnostics <- .summarize_correction_clipping(fold_results, source_sites)
  if (isTRUE(clip_diagnostics$total$any_clipped)) {
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
    avg_source_est = do.call(rbind, lapply(inner_fold_info, function(x) x$avg_source_est))
  )

  results <- list(
    estimate = final_estimate,
    se = se,
    variance = final_variance,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    target_only = list(
      estimate = final_target_estimate,
      se = sd(target_estimates) / sqrt(n_folds)
    ),
    source_estimates = source_estimates,
    weights = weights,
    fold_weights = fold_weights,
    fold_lambdas = fold_lambdas,
    fold_aggregated_estimates = fold_aggregated_estimates,
    aggregation_lambda_rule = lambda_rule,
    clip_diagnostics = clip_diagnostics,
    n_sites = K + 1,
    n_folds = n_folds,
    N_all = N_all,
    all_phi_agg = all_phi_agg,
    intermediates = list(
      target_estimates = target_estimates,
      source_estimates_matrix = source_estimates_matrix,
      phase1 = phase1_summary,
      phase1b = phase1b_summary,
      fold_lambdas = fold_lambdas,
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
