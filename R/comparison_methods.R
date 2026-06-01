# comparison_methods.R - Multi-site comparison estimators
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# All methods estimate the POTENTIAL OUTCOME MEAN:
#   μ¹_t = E_t[Y(1)]  (expected outcome under treatment at target site)
# =============================================================================
#
# This file contains multi-site estimators that aggregate across sites:
#   1. Sample-size weighted (SS)
#   2. Inverse-variance weighted (IVW) with DerSimonian-Laird heterogeneity
#   3. Tilted AIPW (exponential tilting with unpenalized MLE nuisances)
#   4. Federated DR-AIPW (density-ratio corrected, federated)
#   5. Pooled DR-AIPW (density-ratio corrected, centralized)
#   6. run_all_comparisons (orchestrator)
#
# Helper functions are in estimators_helpers.R.
# Target-only estimators are in estimators_target.R.
# Oracle estimator is in estimators_oracle.R.


#' Compute site-level AIPW fits once for naive comparison methods
#'
#' Internal helper used by sample-size and inverse-variance baselines to avoid
#' duplicated nuisance-model fitting work.
#'
#' @param data_split split data by site
#' @param family GLM family
#' @param use_rcal Logical
#' @param use_crossfit Logical
#' @param n_folds Optional number of folds
#' @param A_val Treatment value
#' @return Named list of per-site fit_site_aipw outputs
.fit_site_aipw_all_sites <- function(data_split, family = "binomial",
                                     use_rcal = FALSE, use_crossfit = TRUE,
                                     n_folds = NULL, A_val = 1L) {
  sites <- names(data_split)
  site_fits <- setNames(vector("list", length(sites)), sites)
  for (site in sites) {
    site_data <- data_split[[site]]
    site_fits[[site]] <- fit_site_aipw(
      site_data,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }
  site_fits
}


#' Sample-size adjusted estimator (SS) with correct variance estimation
#' 
#' This estimator computes a sample-size weighted average of site-specific
#' AIPW estimates. By default (\code{variance_method = "bootstrap"}) the standard
#' error is the multiplier (wild) bootstrap, i.e. the fixed-effects sampling variance
#' of the weighted average; the DerSimonian-Laird random-effects (\eqn{\tau^2}) variance
#' is available via \code{variance_method = "analytic"} and is retained in
#' \code{components} for reference.
#' 
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   for variance-valid inference (use_rcal is ignored in this mode).
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param site_fits Optional precomputed per-site outputs from
#'   \\code{fit_site_aipw}. When provided, avoids refitting nuisances.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @return estimate with variance
#' @export
estimate_sample_size_weighted <- function(data_split, family = "binomial",
                                          use_rcal = FALSE, use_crossfit = TRUE,
                                          n_folds = NULL, A_val = 1L,
                                          site_fits = NULL,
                                          variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  sites <- names(data_split)
  site_estimates <- list()
  total_n <- 0

  if (is.null(site_fits)) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }
  
  for (site in sites) {
    res <- site_fits[[site]]
    site_estimates[[site]] <- list(
      estimate = res$estimate,
      variance = res$variance,
      n = res$n,
      V_ot = res$V_ot
    )
    total_n <- total_n + res$n
  }
  
  # ==========================================================================
  # VARIANCE ESTIMATION WITH HETEROGENEITY CORRECTION
  # ==========================================================================
  # Similar to IVW, we use a random-effects model to account for between-site
  # heterogeneity. When sites have different covariate distributions (covariate
  # shift), they may estimate different population parameters, introducing
  # extra variance that the naive pooled-IF formula doesn't capture.
  #
  # Random-effects model:
  #   μ̂_j ~ N(μ, Var_j + τ²)
  # where τ² is estimated via DerSimonian-Laird using Cochran's Q.
  #
  # The random-effects variance for sample-size weighted estimator is:
  #   Var_RE = Σ_j w_j² × (Var_j + τ²)
  # ==========================================================================
  
  K <- length(site_estimates)
  
  # Compute sample-size weighted estimate
  weighted_sum <- 0
  for (site in names(site_estimates)) {
    w <- site_estimates[[site]]$n / total_n
    weighted_sum <- weighted_sum + w * site_estimates[[site]]$estimate
  }
  
  # ==========================================================================
  # Estimate between-site heterogeneity using DerSimonian-Laird
  # ==========================================================================
  site_est_vec <- sapply(site_estimates, function(x) x$estimate)
  site_variances <- sapply(site_estimates, function(x) x$V_ot / x$n)
  het <- calculate_dl_heterogeneity(site_est_vec, site_variances)
  tau_sq <- het$tau_sq
  
  # ==========================================================================
  # Compute variance with heterogeneity correction
  # For sample-size weights: Var_RE = Σ_j (n_j/N)² × (Var_j + τ²)
  # ==========================================================================
  
  var_fe <- 0   # Fixed-effects variance (naive)
  var_re <- 0   # Random-effects variance (with τ²)
  
  for (site in names(site_estimates)) {
    n_j <- site_estimates[[site]]$n
    w_j <- n_j / total_n
    var_j <- site_variances[site]
    
    var_fe <- var_fe + w_j^2 * var_j
    var_re <- var_re + w_j^2 * (var_j + tau_sq)
  }
  
  # Variance: bootstrap (default) over per-site influence blocks with sample-size
  # weights, or the analytic random-effects (tau^2) variance. The bootstrap reproduces
  # the fixed-effects sampling variance and avoids the tau^2 inflation.
  boot_blocks <- lapply(names(site_estimates), function(site) {
    list(influence = site_fits[[site]]$varphi_ot,
         weight = site_estimates[[site]]$n / total_n)
  })
  var_res <- .resolve_comparison_variance(
    analytic_variance = var_re, blocks = boot_blocks, variance_method = variance_method
  )

  return(list(
    estimate = weighted_sum,
    variance = var_res$variance,
    se = var_res$se,
    method = "sample_size",
    n = total_n,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      var_fixed_effects = as.numeric(var_fe),
      var_random_effects = as.numeric(var_re),
      tau_squared = as.numeric(tau_sq),
      Q_statistic = het$Q,
      df = K - 1,
      I_squared = het$I_squared
    )
  ))
}

#' Inverse-variance weighted estimator (IVW) with heterogeneity correction
#' 
#' This estimator computes an inverse-variance weighted average of site-specific
#' AIPW estimates. By default (\code{variance_method = "bootstrap"}) the standard
#' error is the multiplier (wild) bootstrap, i.e. the fixed-effects sampling variance
#' \eqn{1/\sum_j \mathrm{precision}_j}; the random-effects (\eqn{\tau^2}) variance is
#' available via \code{variance_method = "analytic"} and retained in \code{components}.
#'
#' The analytic random-effects path is motivated as follows. When sites have different
#' covariate distributions, they may estimate slightly different population parameters,
#' and the fixed-effects IVW variance (\eqn{1/\sum_j \mathrm{precision}_j}) would then
#' understate uncertainty ABOUT THAT COMMON PARAMETER; Cochran's Q estimates the
#' between-site component (\eqn{\tau^2}). For the fixed target estimand, however, that
#' between-site spread is transport bias rather than sampling variance, so the bootstrap
#' (fixed-effects) SE is the honest default and the \eqn{\tau^2} inflation is opt-in.
#' 
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   for variance-valid inference (use_rcal is ignored in this mode).
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param site_fits Optional precomputed per-site outputs from
#'   \\code{fit_site_aipw}. When provided, avoids refitting nuisances.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @return estimate with variance
#' @export
estimate_inverse_variance_weighted <- function(data_split, family = "binomial",
                                               use_rcal = FALSE, use_crossfit = TRUE,
                                               n_folds = NULL, A_val = 1L,
                                               site_fits = NULL,
                                               variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  sites <- names(data_split)
  site_results <- list()
  total_n <- 0

  if (is.null(site_fits)) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }
  
  for (site in sites) {
    res <- site_fits[[site]]
    safe_variance <- max(res$variance, VARIANCE_MIN)
    safe_theta <- max(res$V_ot, VARIANCE_MIN)
    site_results[[site]] <- list(
      estimate = res$estimate,
      variance = safe_variance,
      precision = 1 / safe_variance,
      theta = safe_theta,
      n = res$n
    )
    total_n <- total_n + res$n
  }
  
  # ==========================================================================
  # VARIANCE ESTIMATION WITH HETEROGENEITY CORRECTION
  # ==========================================================================
  # Standard IVW (fixed-effects) assumes all sites estimate the same parameter.
  # When there's covariate shift, sites may estimate different parameters,
  # introducing between-site heterogeneity that inflates variance.
  #
  # We use a random-effects approach:
  #   μ̂_j ~ N(μ, Var_j + τ²)
  # where τ² is the between-site variance estimated via DerSimonian-Laird.
  #
  # Random-effects variance:
  #   Var_RE = 1 / Σ_j (1 / (Var_j + τ²))
  # ==========================================================================
  
  K <- length(site_results)
  
  # Calculate fixed-effects IVW weights
  total_precision <- sum(sapply(site_results, function(x) x$precision))
  
  # Fixed-effects weighted estimate
  weighted_sum <- 0
  for (site in names(site_results)) {
    w <- site_results[[site]]$precision / total_precision
    weighted_sum <- weighted_sum + w * site_results[[site]]$estimate
  }
  
  # Fixed-effects variance (base variance)
  var_fe <- 1 / total_precision
  
  # ==========================================================================
  # Estimate between-site heterogeneity using DerSimonian-Laird
  # ==========================================================================
  site_est_vec <- sapply(site_results, function(x) x$estimate)
  site_var_vec <- sapply(site_results, function(x) x$variance)
  het <- calculate_dl_heterogeneity(site_est_vec, site_var_vec)
  tau_sq <- het$tau_sq
  
  # ==========================================================================
  # Random-effects variance
  # Var_RE = 1 / Σ_j (1 / (Var_j + τ²))
  # ==========================================================================
  if (tau_sq > 0) {
    # Random-effects weights
    re_precision <- 0
    for (site in names(site_results)) {
      var_j <- site_results[[site]]$variance
      re_precision <- re_precision + 1 / (var_j + tau_sq)
    }
    var_re <- 1 / re_precision
    
    # Re-compute weighted estimate with random-effects weights
    weighted_sum_re <- 0
    for (site in names(site_results)) {
      var_j <- site_results[[site]]$variance
      w_re <- (1 / (var_j + tau_sq)) / re_precision
      weighted_sum_re <- weighted_sum_re + w_re * site_results[[site]]$estimate
    }
    
    # Use random-effects estimate and variance
    final_estimate <- weighted_sum_re
    final_variance <- var_re
  } else {
    # No heterogeneity detected, use fixed-effects
    final_estimate <- weighted_sum
    final_variance <- var_fe
  }
  
  # ==========================================================================
  # Additional pooled IF variance as a check (like sample_size method)
  # This provides an alternative variance estimate based on pooled IFs
  # ==========================================================================
  
  # Compute weights for pooled variance (use RE weights if tau_sq > 0)
  weights <- numeric(K)
  for (idx in seq_along(site_results)) {
    site <- names(site_results)[idx]
    if (tau_sq > 0) {
      var_j <- site_results[[site]]$variance
      re_precision_total <- sum(sapply(site_results, function(x) 1/(x$variance + tau_sq)))
      weights[idx] <- (1 / (var_j + tau_sq)) / re_precision_total
    } else {
      weights[idx] <- site_results[[site]]$precision / total_precision
    }
  }
  
  # Pooled IF variance with between-site component
  within_site_var <- 0
  between_site_var <- 0
  
  for (idx in seq_along(site_results)) {
    site <- names(site_results)[idx]
    w_j <- weights[idx]
    n_j <- site_results[[site]]$n
    theta_j <- site_results[[site]]$theta
    mu_j <- site_results[[site]]$estimate
    
    # Within-site: w_j² × Var_j = w_j² × θ_j / n_j
    within_site_var <- within_site_var + w_j^2 * theta_j / n_j
    
    # Between-site: w_j² × (μ̂_j - μ̂)²
    between_site_var <- between_site_var + w_j^2 * (mu_j - final_estimate)^2
  }
  
  pooled_variance <- within_site_var + between_site_var
  
  # Use the larger of RE variance and pooled variance for robustness (analytic path)
  final_variance <- max(final_variance, pooled_variance)

  # Variance: bootstrap (default) over per-site influence blocks with fixed-effects
  # precision weights, or the analytic random-effects variance above.
  boot_blocks <- lapply(names(site_results), function(site) {
    list(influence = site_fits[[site]]$varphi_ot,
         weight = site_results[[site]]$precision / total_precision)
  })
  var_res <- .resolve_comparison_variance(
    analytic_variance = final_variance, blocks = boot_blocks, variance_method = variance_method
  )

  return(list(
    estimate = final_estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "inverse_variance",
    n = total_n,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      var_fixed_effects = var_fe,
      var_random_effects = if (tau_sq > 0) var_re else var_fe,
      var_pooled = pooled_variance,
      tau_squared = tau_sq,
      Q_statistic = het$Q,
      df = K - 1,
      I_squared = het$I_squared,
      within_site_var = within_site_var,
      between_site_var = between_site_var
    )
  ))
}

#' Tilted AIPW estimator with unpenalized MLE nuisances
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#'        Supports both "binomial" (logit) and "gaussian" (identity).
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @return estimate with variance
#' @export
estimate_tilted_aipw <- function(data_split, family = "binomial", A_val = 1L,
                                 variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)
  glm_spec <- resolve_glm_family(family)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")

  # Target-only estimate using unpenalized MLE for strict variance derivation.
  # validate_algorithm_inputs() already guarantees at least one A == A_val unit
  # in every site, so the inline emptiness checks live in the validator now.
  target_x <- as.matrix(target_data$W_outcome)
  target_y <- as.numeric(target_data$Y)
  target_a <- as.numeric(target_data$A)
  target_n <- target_data$n

  # Fit nuisance models and compute AIPW estimate
  target_est <- tryCatch({
    ps_fit_t <- fit_logit_mle(target_x, target_a)
    or_fit_t <- fit_glm_mle(target_x[target_a == A_val, , drop = FALSE],
                            target_y[target_a == A_val],
                            family = glm_spec$family)
    eta_t <- as.numeric(cbind(1, target_x) %*% or_fit_t$coefficients)
    m_hat_t <- switch(glm_spec$link,
                      "logit" = 1 / (1 + exp(-eta_t)),
                      "identity" = eta_t,
                      eta_t)
    m_hat_t <- clip_outcome_pred(m_hat_t, glm_spec$family)
    target_aipw <- calculate_aipw_influence(target_y, target_a, target_x, m_hat_t,
                                            ps_fit_t$fitted, A_val = A_val,
                                            family = glm_spec$family)
    list(
      estimate = target_aipw$estimate,
      variance = target_aipw$variance,
      varphi_ot = target_aipw$influence
    )
  }, error = function(e) {
    stop(sprintf("estimate_tilted_aipw: target-site AIPW nuisance fitting failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })

  # Calculate source site estimates with unpenalized MLE nuisances
  source_estimates <- list()
  Z_target <- as.matrix(target_data$Z_site)
  mean_phi_target <- c(1, colMeans(Z_target))
  # Pre-compute target-side matrices (constant across source sites)
  Z_target_int <- cbind(1, Z_target)
  Z_target_centered <- sweep(Z_target_int, 2, mean_phi_target, "-")
  
  for (site in source_sites) {
    source_data <- data_split[[site]]
    
    y_source <- source_data$Y
    tr_source <- source_data$A
    x_source <- as.matrix(source_data$W_outcome)
    n_source <- source_data$n

    # validate_algorithm_inputs() ensures every source site has >= 1 unit with
    # A == A_val, so treated_idx is guaranteed non-empty here.
    treated_idx <- which(tr_source == A_val)

    # Fit nuisance models and compute weighted AIPW estimate
    source_estimates[[site]] <- tryCatch({
      X_matrix <- as.matrix(x_source)
      Z_source <- as.matrix(source_data$Z_site)

      # Step 1: Fit outcome model (unpenalized MLE)
      or_fit <- fit_glm_mle(X_matrix[treated_idx, , drop = FALSE],
                as.numeric(y_source[treated_idx]),
                family = glm_spec$family)
      eta_source <- as.numeric(cbind(1, X_matrix) %*% or_fit$coefficients)
      m1_pred <- switch(glm_spec$link,
            "logit" = 1 / (1 + exp(-eta_source)),
            "identity" = eta_source,
            eta_source)
      m1_pred <- clip_outcome_pred(m1_pred, glm_spec$family)

      # Step 2: Fit propensity model (unpenalized MLE)
      ps_fit <- fit_logit_mle(X_matrix, as.numeric(tr_source))
      prop_scores <- clip_propensity(ps_fit$fitted)

      # Step 3: Fit density ratio model via exponential tilting (unpenalized)
      # γ_{s,A_val} is the arm-specific density ratio
      alpha <- fit_initial_density_ratio(Z_source, tr_source, mean_phi_target, lambda = 0.0,
                                         A_val = A_val)

      # Density ratio weights
      Z_int <- cbind(1, Z_source)
      eta_i <- as.numeric(Z_int %*% alpha)
      w_i <- exp(-eta_i)
      w_i <- w_i / mean(w_i)

      # Step 4: AIPW estimate with nuisance-adjusted influence
      aipw_res <- calculate_aipw_influence(as.numeric(y_source), as.numeric(tr_source),
                                         X_matrix, m1_pred, prop_scores, w_i,
                                         A_val = A_val, family = glm_spec$family)

      weighted_estimate <- aipw_res$estimate
      phi_i <- aipw_res$phi

      # Density ratio adjustment for influence function
      d_alpha <- aipw_res$d
      phi_centered <- phi_i - weighted_estimate
      
      # Center source Z by mean_phi_target for consistent score function
      # This ensures source and target influence components are on same scale
      Z_int_centered <- sweep(Z_int, 2, mean_phi_target, "-")
      
      A_s <- -1 / d_alpha * colMeans(as.numeric(w_i) * Z_int_centered * as.numeric(phi_centered))
      # Use A-weighted score to match the density ratio estimating equation
      aw_i <- as.numeric(tr_source == A_val) * as.numeric(w_i)
      M_alpha <- t(Z_int_centered) %*% (Z_int_centered * aw_i) / n_source
      adj_alpha <- solve_with_ridge(M_alpha) %*% A_s
      infl_alpha_source <- as.numeric(aw_i * (Z_int_centered %*% adj_alpha))

      influence <- aipw_res$influence + infl_alpha_source
      weighted_variance <- mean(influence^2) / n_source

      # Target-side influence component for cross-site covariance
      target_if_component <- -as.numeric(Z_target_centered %*% adj_alpha)

      # Return the result (will be assigned to source_estimates[[site]])
      list(
        estimate = weighted_estimate,
        variance = weighted_variance,
        varphi_ot = influence,
        weights = w_i,
        n = n_source,
        A_s = A_s,
        alpha = alpha,
        d_alpha = d_alpha,
        target_if_component = target_if_component
      )
      
    }, error = function(e) {
      stop(sprintf("estimate_tilted_aipw: source site '%s' nuisance fitting failed: %s",
                   site, conditionMessage(e)), call. = FALSE)
    })
  }
  
  # Aggregate using sample size weights
  total_n <- target_data$n + sum(sapply(source_estimates, function(x) x$n))
  
  # Calculate sample size weights
  w_target <- target_data$n / total_n
  w_sources <- sapply(source_estimates, function(x) x$n / total_n)
  
  # Aggregated estimate
  weighted_sum <- w_target * target_est$estimate
  
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    weighted_sum <- weighted_sum + w * source_estimates[[site]]$estimate
  }
  
  # ============================================================
  # Calculate variance for tilted AIPW using influence functions
  # ============================================================
  # The aggregated estimator is a weighted sum of target-only and source
  # estimators. We combine target-side IF components (shared target data)
  # and source-side IF components (site-specific data) to capture
  # cross-site covariance induced by shared target information.
  # ============================================================
  
  target_if <- target_est$varphi_ot
  if (length(target_if) != target_n) {
    target_if <- rep(0, target_n)
  }
  
  aggregated_target_if <- w_target * target_if
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    target_if_component <- source_estimates[[site]]$target_if_component
    if (is.null(target_if_component) || length(target_if_component) != target_n) {
      target_if_component <- rep(0, target_n)
    }
    aggregated_target_if <- aggregated_target_if + w * target_if_component
  }
  
  aggregated_target_if <- aggregated_target_if - mean(aggregated_target_if)
  var_target_component <- mean(aggregated_target_if^2) / max(1, target_n)
  
  var_source_component <- 0
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    n_source <- source_estimates[[site]]$n
    source_if <- source_estimates[[site]]$varphi_ot
    if (!is.null(source_if) && length(source_if) > 0 && n_source > 0) {
      source_if <- source_if - mean(source_if)
      var_source_component <- var_source_component +
        (w^2) * mean(source_if^2) / n_source
    }
  }
  
  var_total <- max(var_target_component + var_source_component, VARIANCE_MIN)

  # Variance: bootstrap (default) over the shared-target influence block plus per-source
  # influence blocks, or the analytic influence-function variance above.
  boot_blocks <- c(
    list(list(influence = aggregated_target_if, weight = 1)),
    lapply(names(source_estimates), function(site) {
      list(influence = source_estimates[[site]]$varphi_ot,
           weight = source_estimates[[site]]$n / total_n)
    })
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = var_total, blocks = boot_blocks, variance_method = variance_method
  )

  return(list(
    estimate = weighted_sum,
    variance = var_res$variance,
    se = var_res$se,
    method = "tilted_aipw",
    n = total_n,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      target_weight = w_target,
      source_weights = w_sources,
      target_variance_component = var_target_component,
      source_variance_component = var_source_component
    )
  ))
}


#' Federated DR-AIPW estimator
#' 
#' Each source site uses density ratio weighting to estimate E_t[Y(1)].
#' All sites estimate the SAME estimand, then aggregate with IVW.
#' This is the theoretically correct federated baseline.
#' 
#' @param data_split Split data by site
#' @param dr_lambda Regularization for density ratio. If NULL, selected via CV.
#' @param dr_lambda_rule CV selection rule when \code{dr_lambda = NULL}:
#'   \code{"min"} (default) selects \code{lambda.min}; \code{"1se"} selects
#'   \code{lambda.1se}.
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @return List with estimate, variance, se
#' @export
estimate_federated_dr <- function(data_split, dr_lambda = NULL,
                                  dr_lambda_rule = c("min", "1se"),
                                  A_val = 1L, family = "binomial",
                                  variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)
  dr_lambda_rule <- .match_nuisance_lambda_rule(
    dr_lambda_rule, "estimate_federated_dr", arg = "dr_lambda_rule"
  )
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  Z_target <- as.matrix(target_data$Z_site)
  n_target <- target_data$n
  mean_phi_target <- c(1, colMeans(Z_target))
  Z_target_int <- cbind(1, Z_target)
  Z_target_centered <- sweep(Z_target_int, 2, mean_phi_target, "-")
  
  site_results <- list()
  
  # Target site: standard AIPW
  target_res <- calculate_weighted_site_aipw(
    y = target_data$Y, a = target_data$A,
    X = as.matrix(target_data$W_outcome), weights = NULL, family = family, A_val = A_val
  )
  site_results[["t"]] <- list(
    estimate = target_res$estimate,
    variance = target_res$variance,
    n = target_data$n,
    varphi_ot = target_res$varphi_ot
  )
  
  # Source sites: DR-weighted AIPW
  for (site in source_sites) {
    source_data <- data_split[[site]]
    Z_source <- as.matrix(source_data$Z_site)
    
    dr_weights <- calculate_dr_weights(
      Z_source, Z_target, lambda = dr_lambda, lambda_rule = dr_lambda_rule
    )
    
    source_res <- calculate_weighted_site_aipw(
      y = source_data$Y, a = source_data$A,
      X = as.matrix(source_data$W_outcome), weights = dr_weights, family = family, A_val = A_val
    )

    # First-order correction for density-ratio estimation effect.
    # Mirrors the score-adjustment structure used in tilted_aipw.
    n_source <- source_data$n
    Z_int <- cbind(1, Z_source)
    Z_int_centered <- sweep(Z_int, 2, mean_phi_target, "-")
    d_alpha <- max(mean(dr_weights), DIVISION_FLOOR)
    phi_centered <- source_res$phi - source_res$estimate
    A_s <- -1 / d_alpha * colMeans(as.numeric(dr_weights) * Z_int_centered * as.numeric(phi_centered))
    aw_i <- as.numeric(source_data$A == A_val) * as.numeric(dr_weights)
    M_alpha <- t(Z_int_centered) %*% (Z_int_centered * aw_i) / max(1, n_source)
    adj_alpha <- solve_with_ridge(M_alpha) %*% A_s
    infl_alpha_source <- as.numeric(aw_i * (Z_int_centered %*% adj_alpha))
    target_if_component <- -as.numeric(Z_target_centered %*% adj_alpha)
    varphi_source <- as.numeric(source_res$varphi_ot) + infl_alpha_source
    varphi_source <- varphi_source - mean(varphi_source)
    source_variance <- mean(varphi_source^2) / max(1, n_source)
    
    site_results[[site]] <- list(
      estimate = source_res$estimate,
      variance = source_variance,
      n = source_data$n,
      varphi_ot = varphi_source,
      target_if_component = target_if_component
    )
  }
  
  # Filter valid sites and aggregate with IVW
  valid_sites <- names(site_results)[sapply(site_results, function(x) 
    !is.na(x$estimate) && is.finite(x$variance) && x$variance > 0)]
  
  if (length(valid_sites) == 0) {
    return(list(estimate = NA, variance = Inf, se = Inf, method = "federated_dr"))
  }
  
  precisions <- sapply(site_results[valid_sites], function(x) 1/x$variance)
  total_precision <- sum(precisions)
  ivw_weights <- precisions / total_precision
  
  estimate <- sum(ivw_weights * sapply(site_results[valid_sites], function(x) x$estimate))
  
  # Shared-target + source-side variance decomposition (same style as main methods)
  weight_by_site <- setNames(as.numeric(ivw_weights), valid_sites)
  aggregated_target_if <- rep(0, n_target)
  
  if ("t" %in% valid_sites) {
    aggregated_target_if <- aggregated_target_if + weight_by_site[["t"]] * site_results[["t"]]$varphi_ot
  }
  for (site in setdiff(valid_sites, "t")) {
    target_if_component <- site_results[[site]]$target_if_component
    if (!is.null(target_if_component) && length(target_if_component) == n_target) {
      aggregated_target_if <- aggregated_target_if + weight_by_site[[site]] * target_if_component
    }
  }
  aggregated_target_if <- aggregated_target_if - mean(aggregated_target_if)
  var_target_component <- mean(aggregated_target_if^2) / max(1, n_target)

  var_source_component <- 0
  for (site in setdiff(valid_sites, "t")) {
    w_site <- weight_by_site[[site]]
    n_site <- site_results[[site]]$n
    varphi_site <- site_results[[site]]$varphi_ot
    if (!is.null(varphi_site) && length(varphi_site) > 0) {
      varphi_site <- varphi_site - mean(varphi_site)
      var_source_component <- var_source_component + w_site^2 * mean(varphi_site^2) / max(1, n_site)
    }
  }

  variance <- max(var_target_component + var_source_component, VARIANCE_MIN)

  # Variance: bootstrap (default) over the shared-target influence block plus per-source
  # influence blocks (IVW weights), or the analytic influence-function variance above.
  boot_blocks <- c(
    list(list(influence = aggregated_target_if, weight = 1)),
    lapply(setdiff(valid_sites, "t"), function(site) {
      list(influence = site_results[[site]]$varphi_ot,
           weight = weight_by_site[[site]])
    })
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = variance, blocks = boot_blocks, variance_method = variance_method
  )

  return(list(
    estimate = estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "federated_dr",
    n = sum(sapply(site_results[valid_sites], function(x) x$n)),
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      ivw_weights = weight_by_site,
      var_target_component = var_target_component,
      var_source_component = var_source_component
    )
  ))
}

#' Pooled DR-AIPW estimator (centralized version)
#' 
#' Pools all data and uses density ratio weighting for each observation.
#' This is a centralized (non-federated) method for comparison.
#' 
#' @param data_split Split data by site
#' @param dr_lambda Regularization for density ratio. If NULL, selected via CV.
#' @param dr_lambda_rule CV selection rule when \code{dr_lambda = NULL}:
#'   \code{"min"} (default) selects \code{lambda.min}; \code{"1se"} selects
#'   \code{lambda.1se}.
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @return List with estimate, variance, se
#' @export
estimate_pooled_dr <- function(data_split, dr_lambda = NULL,
                               dr_lambda_rule = c("min", "1se"),
                               A_val = 1L, family = "binomial",
                               variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)
  dr_lambda_rule <- .match_nuisance_lambda_rule(
    dr_lambda_rule, "estimate_pooled_dr", arg = "dr_lambda_rule"
  )
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  Z_target <- as.matrix(target_data$Z_site)
  n_target <- target_data$n
  mean_phi_target <- c(1, colMeans(Z_target))
  Z_target_int <- cbind(1, Z_target)
  Z_target_centered <- sweep(Z_target_int, 2, mean_phi_target, "-")
  
  # Pool all data
  all_y <- target_data$Y
  all_a <- target_data$A
  all_X <- as.matrix(target_data$W_outcome)
  all_weights <- rep(1, target_data$n)
  all_sites <- rep("t", target_data$n)
  source_meta <- list()
  row_cursor <- target_data$n
  
  for (site in source_sites) {
    source_data <- data_split[[site]]
    Z_source <- as.matrix(source_data$Z_site)
    
    dr_weights <- calculate_dr_weights(
      Z_source, Z_target, lambda = dr_lambda, lambda_rule = dr_lambda_rule
    )

    source_res <- calculate_weighted_site_aipw(
      y = source_data$Y, a = source_data$A,
      X = as.matrix(source_data$W_outcome), weights = dr_weights, family = family, A_val = A_val
    )

    # First-order correction for density-ratio estimation uncertainty
    # (source-side and shared-target components).
    n_source <- source_data$n
    Z_int <- cbind(1, Z_source)
    Z_int_centered <- sweep(Z_int, 2, mean_phi_target, "-")
    d_alpha <- max(mean(dr_weights), DIVISION_FLOOR)
    phi_centered <- source_res$phi - source_res$estimate
    A_s <- -1 / d_alpha * colMeans(as.numeric(dr_weights) * Z_int_centered * as.numeric(phi_centered))
    aw_i <- as.numeric(source_data$A == A_val) * as.numeric(dr_weights)
    M_alpha <- t(Z_int_centered) %*% (Z_int_centered * aw_i) / max(1, n_source)
    adj_alpha <- solve_with_ridge(M_alpha) %*% A_s
    infl_alpha_source <- as.numeric(aw_i * (Z_int_centered %*% adj_alpha))
    target_if_component <- -as.numeric(Z_target_centered %*% adj_alpha)
    
    idx_start <- row_cursor + 1L
    idx_end <- row_cursor + n_source
    row_cursor <- idx_end
    source_meta[[site]] <- list(
      idx = idx_start:idx_end,
      infl_alpha_source = infl_alpha_source,
      target_if_component = target_if_component,
      lambda = as.numeric(attr(dr_weights, "lambda_used") %||% dr_lambda %||% NA_real_)
    )
    
    all_y <- c(all_y, source_data$Y)
    all_a <- c(all_a, source_data$A)
    all_X <- rbind(all_X, as.matrix(source_data$W_outcome))
    all_weights <- c(all_weights, dr_weights)
    all_sites <- c(all_sites, rep(site, n_source))
  }
  
  result <- calculate_weighted_site_aipw(y = all_y, a = all_a, X = all_X, weights = all_weights, family = family, A_val = A_val)

  # Stacked IF: base weighted AIPW IF + DR-weight estimation correction
  N_all <- length(all_y)
  total_correction <- rep(0, N_all)
  for (site in names(source_meta)) {
    meta <- source_meta[[site]]
    total_correction[meta$idx] <- total_correction[meta$idx] + meta$infl_alpha_source
    total_correction[seq_len(n_target)] <- total_correction[seq_len(n_target)] + meta$target_if_component
  }
  varphi_total <- as.numeric(result$varphi_ot) + total_correction

  # Site-stratified variance (same style as main cross-fitting estimators):
  #   Var = (1/N_all^2) * sum_g sum_{i in g} (varphi_i - mean_g(varphi))^2
  wss <- 0
  for (site in unique(all_sites)) {
    idx <- which(all_sites == site)
    if (length(idx) == 0) next
    varphi_g <- varphi_total[idx]
    varphi_g <- varphi_g - mean(varphi_g)
    wss <- wss + sum(varphi_g^2)
  }
  final_variance <- max(wss / (N_all^2), VARIANCE_MIN)

  # Variance: bootstrap (default) over per-site blocks of the stacked influence function
  # (one block per site, weight n_site/N_all so the wild bootstrap reproduces the
  # site-stratified within-group variance above), or that analytic variance.
  boot_blocks <- lapply(unique(all_sites), function(site) {
    idx <- which(all_sites == site)
    list(influence = varphi_total[idx], weight = length(idx) / N_all)
  })
  var_res <- .resolve_comparison_variance(
    analytic_variance = final_variance, blocks = boot_blocks, variance_method = variance_method
  )

  source_lambdas <- setNames(
    sapply(source_sites, function(site) {
      meta <- source_meta[[site]]
      if (is.null(meta)) NA_real_ else as.numeric(meta$lambda)
    }),
    source_sites
  )
  
  return(list(
    estimate = result$estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "pooled_dr",
    n = N_all,
    varphi_ot = varphi_total,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_target = n_target,
      n_source_total = N_all - n_target,
      mean_weight = mean(all_weights),
      source_lambdas = source_lambdas,
      variance_wss = wss
    )
  ))
}

#' Run all comparison methods
#' 
#' Returns both:
#' - Naive methods (sample_size, inverse_variance): may have bias due to covariate shift
#' - DR-corrected methods (federated_dr, pooled_dr): theoretically correct
#' 
#' @param data_split split data by site
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances.
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param variance_method Standard-error method for the comparison baselines
#'   (\code{sample_size}, \code{inverse_variance}, \code{federated_dr},
#'   \code{pooled_dr}, \code{tilted_aipw}): \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap; \code{"analytic"} uses each method's analytic
#'   variance. \code{target_only} always uses its analytic influence-function variance.
#' @return Named list of per-method results. The five comparison baselines additionally
#'   carry \code{components$variance_method}, \code{components$variance_analytic},
#'   \code{components$se_analytic}, and \code{components$se_bootstrap}; the
#'   \code{target_only} benchmark does not. Consumers should read \code{$estimate} and
#'   \code{$se} uniformly and treat those four \code{components} fields as baseline-only.
#' @export
run_all_comparisons <- function(data_split, use_rcal = FALSE,
                                use_crossfit = TRUE, n_folds = NULL,
                                family = "binomial", A_val = 1L,
                                variance_method = c("bootstrap", "analytic")) {
  variance_method <- match.arg(variance_method)

  precomputed_site_fits <- .fit_site_aipw_all_sites(
    data_split = data_split,
    family = family,
    use_rcal = use_rcal,
    use_crossfit = use_crossfit,
    n_folds = n_folds,
    A_val = A_val
  )

  results <- list(
    # Target-only (benchmark; analytic influence-function variance)
    target_only = estimate_target_only(
      data_split, family, use_rcal = use_rcal,
      use_crossfit = use_crossfit, n_folds = n_folds, A_val = A_val
    ),

    # Naive methods (may have bias due to different estimands)
    sample_size = estimate_sample_size_weighted(
      data_split, family, use_rcal = use_rcal,
      use_crossfit = use_crossfit, n_folds = n_folds,
      A_val = A_val, site_fits = precomputed_site_fits,
      variance_method = variance_method
    ),
    inverse_variance = estimate_inverse_variance_weighted(
      data_split, family, use_rcal = use_rcal,
      use_crossfit = use_crossfit, n_folds = n_folds,
      A_val = A_val, site_fits = precomputed_site_fits,
      variance_method = variance_method
    ),

    # DR-corrected methods (theoretically correct, lambda selected via CV)
    federated_dr = estimate_federated_dr(
      data_split, dr_lambda = NULL, A_val = A_val, family = family,
      variance_method = variance_method
    ),
    pooled_dr = estimate_pooled_dr(
      data_split, dr_lambda = NULL, A_val = A_val, family = family,
      variance_method = variance_method
    ),

    # Tilted AIPW
    tilted_aipw = estimate_tilted_aipw(data_split, family, A_val = A_val,
                                       variance_method = variance_method)
  )

  return(results)
}
