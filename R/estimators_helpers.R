# estimators_helpers.R - Shared helper functions for estimators
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# All methods estimate the POTENTIAL OUTCOME MEAN:
#   μ¹_t = E_t[Y(1)]  (expected outcome under treatment at target site)
# =============================================================================
#
# Contents:
#   1. safe_solve, calculate_aipw_pseudo_outcome, fit_logit_mle
#   2. calculate_aipw_influence
#   3. fit_site_aipw (nuisance model fitting helper)
#   4. calculate_dl_heterogeneity (DerSimonian-Laird)
#   5. calculate_dr_weights, calculate_weighted_site_aipw

#' Safe matrix inversion with ridge stabilization
#' @param mat square matrix
#' @param ridge ridge term to stabilize inversion
#' @return inverse of mat with ridge regularization
safe_solve <- function(mat, ridge = RIDGE_DEFAULT) {
  mat_reg <- mat + diag(ridge, nrow(mat))
  tryCatch({
    solve(mat_reg)
  }, error = function(e) {
    qr.solve(mat_reg)
  })
}

#' Compute AIPW pseudo-outcomes for potential outcome estimation
#'
#' Unified helper for the doubly robust pseudo-outcome:
#'   phi_i = m(X_i) + I(A_i = a) * (Y_i - m(X_i)) / P(A = a | X_i)
#'
#' @param y Outcome vector
#' @param a Treatment vector
#' @param m_hat Outcome model predictions E[Y|A=a, X]
#' @param p_a Propensity P(A=a|X), already computed for the correct arm
#' @param A_val Treatment value (default 1)
#' @return Numeric vector of pseudo-outcomes
calculate_aipw_pseudo_outcome <- function(y, a, m_hat, p_a, A_val = 1L) {
  p_a <- pmax(p_a, PROP_SCORE_LOWER)
  indicator <- as.numeric(a == A_val)
  phi <- m_hat + indicator * (y - m_hat) / p_a
  phi
}

#' Fit unpenalized GLM MLE with intercept
#' @param X covariate matrix (n x p)
#' @param y response vector
#' @param family GLM family ("binomial" or "gaussian")
#' @return list with coefficients and fitted values on response scale
fit_glm_mle <- function(X, y, family = "binomial") {
  glm_spec <- resolve_glm_family(family)
  X_int <- cbind(1, X)
  fit <- glm.fit(
    x = X_int,
    y = y,
    family = switch(glm_spec$family,
      "binomial" = binomial(),
      "gaussian" = gaussian()
    )
  )
  if (any(is.na(fit$coefficients))) {
    stop(sprintf(
      "fit_glm_mle: GLM produced %d NA coefficient(s) out of %d. "
      , sum(is.na(fit$coefficients)), length(fit$coefficients)),
      "This may indicate perfect separation or collinearity (n=",
      nrow(X), ", p=", ncol(X), ").")
  }
  eta <- as.numeric(X_int %*% fit$coefficients)
  fitted <- switch(glm_spec$link,
    "logit" = 1 / (1 + exp(-eta)),
    "identity" = eta,
    fit$fitted.values
  )
  list(coefficients = fit$coefficients, fitted = fitted, converged = fit$converged)
}

#' Fit unpenalized logit MLE with intercept
#' @param X covariate matrix (n x p)
#' @param y binary response
#' @return list with coefficients and fitted probabilities
fit_logit_mle <- function(X, y) {
  fit_glm_mle(X = X, y = y, family = "binomial")
}

#' Calculate AIPW estimate and influence function with nuisance adjustment
#' @param y outcome vector
#' @param a treatment vector
#' @param x covariate matrix for outcome/propensity
#' @param m_hat outcome model predictions
#' @param pi_hat propensity scores P(A=1|X)
#' @param w optional density ratio weights (defaults to 1)
#' @param A_val treatment value (1 or 0, default 1)
#' @param family GLM family string for outcome clipping/variance ("binomial", "gaussian", etc.)
#' @return list with estimate, influence function, variance, and intermediate terms
calculate_aipw_influence <- function(y, a, x, m_hat, pi_hat, w = NULL, A_val = 1L,
                                     family = "binomial") {
  n <- length(y)
  X_int <- cbind(1, x)
  if (is.null(w)) {
    w <- rep(1, n)
  }
  pi_hat <- clip_propensity(pi_hat)
  m_hat <- clip_outcome_pred(m_hat, family)

  # General AIPW for E[Y(a)]:
  #   A_val=1: phi = m + A*(Y-m)/pi
  #   A_val=0: phi = m + (1-A)*(Y-m)/(1-pi)
  indicator <- as.numeric(a == A_val)
  p_a <- if (A_val == 1) pi_hat else (1 - pi_hat)
  p_a <- pmax(p_a, PROP_SCORE_LOWER)  # safety floor

  phi <- m_hat + indicator * (y - m_hat) / p_a
  d <- mean(w)
  mu_hat <- mean(w * phi) / d

  # Response derivative h'(η) = dμ/dη — depends on GLM family/link.
  # For canonical links, h'(η) = V(μ) (GLM variance function).
  # Used for outcome model adjustment (influence function correction).
  m_prime <- switch(family,
    "binomial" = m_hat * (1 - m_hat),      # V(μ) = μ(1-μ) for binomial
    "gaussian" = rep(1, length(m_hat)),     # V(μ) = 1 for gaussian
    m_hat * (1 - m_hat)                     # default: binomial
  )

  # Outcome model adjustment
  s_beta <- X_int * (indicator * (y - m_hat))
  A_beta <- colMeans(w * (1 - indicator / p_a) * m_prime * X_int) / d
  M_beta <- t(X_int) %*% (X_int * (indicator * m_prime)) / n
  adj_beta <- safe_solve(M_beta) %*% A_beta
  infl_beta <- as.numeric(s_beta %*% adj_beta)

  # Propensity model adjustment
  # NOTE: s_gamma_ps / adj_gamma_ps refer to the PROPENSITY SCORE model
  # coefficients, NOT the density ratio γ_{s_j,a} used in FACE-HD.
  # This function is for standard AIPW (target-only / comparison methods).
  # Score function for propensity: s_gamma_ps = X * (a - pi)
  s_gamma_ps <- X_int * (as.numeric(a) - pi_hat)
  # Sensitivity of phi to pi:
  #   A_val=1: d(phi)/d(pi) = -A*(Y-m)/pi^2, times d(pi)/d(gamma) = pi*(1-pi)*X
  #     => A * (Y-m) * (-(1-pi)/pi) * X
  #   A_val=0: d(phi)/d(pi) = (1-A)*(Y-m)/(1-pi)^2, times pi*(1-pi)*X
  #     => (1-A) * (Y-m) * (pi/(1-pi)) * X
  if (A_val == 1) {
    dphi_term <- -indicator * (y - m_hat) * (1 - pi_hat) / pi_hat
  } else {
    dphi_term <- indicator * (y - m_hat) * pi_hat / (1 - pi_hat)
  }
  A_gamma_ps <- colMeans(w * dphi_term * X_int) / d
  M_gamma_ps <- t(X_int) %*% (X_int * (pi_hat * (1 - pi_hat))) / n
  adj_gamma_ps <- safe_solve(M_gamma_ps) %*% A_gamma_ps
  infl_gamma_ps <- as.numeric(s_gamma_ps %*% adj_gamma_ps)

  influence <- w * (phi - mu_hat) / d - infl_beta - infl_gamma_ps
  # NOTE: Use mean(IF^2)/n instead of var(IF)/n to avoid n vs n-1 bias
  variance <- mean(influence^2) / n

  list(
    estimate = mu_hat,
    influence = influence,
    variance = variance,
    phi = phi,
    X_int = X_int,
    m_hat = m_hat,
    pi_hat = pi_hat,
    d = d
  )
}

# =============================================================================
# SHARED GLMNET CV HELPER
# =============================================================================

#' Fit L1-penalized GLM via glmnet cross-validation and predict on new data
#'
#' Shared helper that encapsulates the repeated pattern:
#'   \code{cv.glmnet(x_train, y_train) -> predict(x_predict) -> clip}
#' with automatic fold count selection and error handling.
#'
#' @param x_train Training covariate matrix (n_train x p).
#' @param y_train Training response vector (length n_train).
#' @param x_predict Prediction covariate matrix (n_predict x p). Columns must
#'   match \code{x_train}; number of rows may differ.
#' @param family glmnet family string (\code{"binomial"}, \code{"gaussian"}, etc.).
#' @param fallback_value Numeric scalar: if fitting fails, return
#'   \code{rep(fallback_value, nrow(x_predict))}. \code{NULL} returns \code{NA}s.
#' @param clip_fn Optional unary function applied element-wise to the
#'   predictions (e.g., \code{clip_propensity}, \code{clip_outcome_pred}).
#' @param nlambda Number of lambda grid points for \code{cv.glmnet}.
#' @param min_per_fold Minimum observations per fold passed to
#'   \code{get_cv_fold_count} (default 10).
#' @param caller_name Caller function name for warning messages.
#' @param model_name Short model label for warning messages
#'   (e.g., \code{"PS"}, \code{"OR"}).
#' @param fold_id Optional integer fold index included in the warning.
#' @return Numeric vector of predictions (length \code{nrow(x_predict)}).
#' @keywords internal
fit_glmnet_cv <- function(x_train, y_train, x_predict,
                          family = "binomial",
                          fallback_value = NULL,
                          clip_fn = NULL,
                          nlambda = LAMBDA_GRID_SIZE_FAST,
                          min_per_fold = 10L,
                          caller_name = "", model_name = "",
                          fold_id = NULL) {
  n_predict <- nrow(x_predict)

  pred <- tryCatch({
    n_cv_folds <- get_cv_fold_count(nrow(x_train), min_per_fold = min_per_fold)
    if (n_cv_folds < 4L) {
      stop(sprintf("insufficient sample size for cv.glmnet: n_train=%d yields nfolds=%d (<4 minimum valid folds)",
                   nrow(x_train), n_cv_folds))
    }
    cv_fit <- glmnet::cv.glmnet(
      x = x_train, y = y_train, family = family,
      alpha = 1, nfolds = n_cv_folds, nlambda = nlambda
    )
    as.numeric(predict(cv_fit, newx = x_predict,
                       s = "lambda.min", type = "response"))
  }, error = function(e) {
    fold_msg <- if (!is.null(fold_id)) sprintf(" on fold %d", fold_id) else ""
    warning(sprintf("%s: %s model failed%s (%s); using small-sample fallback prediction.",
                    caller_name, model_name, fold_msg, conditionMessage(e)))
    rep(if (!is.null(fallback_value)) fallback_value else NA_real_, n_predict)
  })

  if (!is.null(clip_fn)) pred <- clip_fn(pred)
  pred
}

# =============================================================================
# SHARED NUISANCE MODEL FITTING HELPER
# =============================================================================

#' Fit nuisance models (propensity + outcome) and compute AIPW estimate for a single site
#'
#' This helper extracts the common pattern of:
#'   1. Fit propensity score model P(A=1|X)
#'   2. Fit outcome regression model E[Y|A=1,X]
#'   3. Compute doubly robust AIPW pseudo-outcomes and estimate
#'
#' @param site_data List with components Y, A, W_outcome, n
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Whether to use RCAL (default FALSE, uses glmnet)
#' @param use_crossfit Whether to use cross-fitting (default TRUE)
#' @param n_folds Number of cross-fitting folds (NULL = data-driven)
#' @param A_val Treatment value to estimate potential outcome for (default 1)
#' @return List with estimate, variance, se, phi_i, V_ot, varphi_ot, n,
#'         prop_scores_range, n_treated. Returns NULL on total failure.
fit_site_aipw <- function(site_data, family = "binomial", use_rcal = FALSE,
                          use_crossfit = TRUE, n_folds = NULL, A_val = 1L) {
  glm_spec <- resolve_glm_family(family)
  y <- site_data$Y
  tr <- site_data$A
  x <- as.matrix(site_data$W_outcome)
  n <- site_data$n
  
  treated_idx <- which(tr == A_val)
  if (length(treated_idx) == 0) {
    return(list(
      estimate = 0, variance = 1, se = 1, n = n,
      phi_i = rep(0, n), V_ot = 1,
      varphi_ot = rep(0, n),
      prop_scores_range = c(NA, NA), n_treated = 0
    ))
  }
  
  # Try cross-fitting first (preferred for valid inference)
  if (isTRUE(use_crossfit)) {
    n_cv_folds <- if (is.null(n_folds)) get_cv_fold_count(n) else n_folds
    cf_res <- tryCatch({
      estimate_target_only_crossfit(site_data, n_folds = n_cv_folds, family = family, A_val = A_val)
    }, error = function(e) {
      warning(sprintf("fit_site_aipw: cross-fitted estimation failed (%s), falling back to non-cross-fitted estimator.",
                      conditionMessage(e)))
      NULL
    })
    if (!is.null(cf_res)) {
      phi_i <- if (!is.null(cf_res$varphi_ot)) {
        cf_res$varphi_ot + cf_res$estimate
      } else {
        rep(cf_res$estimate, n)
      }
      V_ot <- mean((phi_i - cf_res$estimate)^2)
      return(list(
        estimate = as.numeric(cf_res$estimate),
        variance = as.numeric(cf_res$variance),
        se = as.numeric(sqrt(cf_res$variance)),
        n = n,
        phi_i = phi_i,
        V_ot = V_ot,
        varphi_ot = as.numeric(cf_res$varphi_ot),
        prop_scores_range = if (!is.null(cf_res$prop_scores)) range(cf_res$prop_scores) else c(NA, NA),
        n_treated = length(treated_idx)
      ))
    }
  }
  
  # Non-cross-fitted estimation (fallback or when use_crossfit=FALSE)
  tryCatch({
    X_matrix <- as.matrix(x)
    if (any(is.na(X_matrix)) || any(is.na(tr)) || any(is.na(y))) {
      stop(sprintf("fit_site_aipw: Missing values detected — X has %d NA(s), A has %d NA(s), Y has %d NA(s). Remove or impute before calling.",
                   sum(is.na(X_matrix)), sum(is.na(tr)), sum(is.na(y))))
    }
    n_cv_folds <- get_cv_fold_count(n)
    
    # Step 1: Fit propensity score model P(A=1|X)
    prop_scores <- NULL
    if (use_rcal && requireNamespace("RCAL", quietly = TRUE)) {
      prop_scores <- tryCatch({
        ps_result <- RCAL::glm.regu.cv(
          fold = n_cv_folds, y = as.numeric(tr), x = X_matrix,
          loss = "cal", nrho = LAMBDA_GRID_SIZE_FAST
        )
        if (!is.null(ps_result$sel.fit) && !any(is.na(ps_result$sel.fit[, 1]))) {
          ps_result$sel.fit[, 1]
        } else { NULL }
      }, error = function(e) {
        warning(sprintf("fit_site_aipw: RCAL propensity score model failed (%s), falling back to glmnet.", conditionMessage(e)))
        NULL
      })
    }
    if (is.null(prop_scores)) {
      cv_fit <- glmnet::cv.glmnet(
        x = X_matrix, y = as.numeric(tr), family = "binomial",
        alpha = 1, nfolds = n_cv_folds, nlambda = LAMBDA_GRID_SIZE_FAST
      )
      prop_scores <- as.numeric(predict(cv_fit, newx = X_matrix,
                                        s = "lambda.min", type = "response"))
    }
    prop_scores <- clip_propensity(prop_scores)
    
    # Step 2: Fit outcome regression model E[Y|A=a, X]
    X_treated <- X_matrix[treated_idx, , drop = FALSE]
    y_treated <- as.numeric(y[treated_idx])
    n_cv_folds <- get_cv_fold_count(length(y_treated))
    
    m1_pred <- NULL
    if (use_rcal && requireNamespace("RCAL", quietly = TRUE)) {
      or_result <- tryCatch({
        loss_type <- switch(family,
          "binomial" = "ml",
          "gaussian" = "gaus",
          "gaus"
        )
        RCAL::glm.regu.cv(
          fold = n_cv_folds, y = y_treated, x = X_treated,
          loss = loss_type, nrho = LAMBDA_GRID_SIZE_FAST
        )
      }, error = function(e) {
        warning(sprintf("fit_site_aipw: RCAL outcome regression failed (%s), falling back to glmnet.", conditionMessage(e)))
        NULL
      })
      if (!is.null(or_result) && !is.null(or_result$sel.bet)) {
        beta_coef <- or_result$sel.bet[, 1]
        X_design <- cbind(1, X_matrix)
        if (ncol(X_design) == length(beta_coef)) {
          m1_pred <- as.numeric(X_design %*% beta_coef)
        }
      }
    }
    if (is.null(m1_pred)) {
      family_type <- glm_spec$glmnet_family
      cv_fit_or <- glmnet::cv.glmnet(
        x = X_treated, y = y_treated, family = family_type,
        alpha = 1, nfolds = n_cv_folds, nlambda = LAMBDA_GRID_SIZE_FAST
      )
      beta_coef <- as.vector(coef(cv_fit_or, s = "lambda.min"))
      X_design <- cbind(1, X_matrix)
      m1_pred <- as.numeric(X_design %*% beta_coef)
    }
    # Apply response function: linear predictor → response scale
    # RCAL and manual coef() extraction return linear predictors (η = X^T β),
    # but AIPW needs predictions on the response scale (μ = h(η)).
    # Note: glmnet predict(..., type="response") already applies this transform,
    # but the manual extraction path above does not.
    if (family == "binomial") {
      m1_pred <- 1 / (1 + exp(-m1_pred))         # logistic: μ = 1/(1+e^{-η})
    }
    # gaussian/identity: no transform needed (μ = η)
    m1_pred <- clip_outcome_pred(m1_pred, family)
    
    # Step 3: Compute doubly robust AIPW pseudo-outcomes
    p_a <- if (A_val == 1) prop_scores else (1 - prop_scores)
    phi_i <- calculate_aipw_pseudo_outcome(as.numeric(y), tr, m1_pred, p_a, A_val = A_val)
    
    if (any(is.na(phi_i)) || any(is.infinite(phi_i))) {
      phi_i[is.na(phi_i) | is.infinite(phi_i)] <- mean(y[treated_idx])
    }
    
    estimate <- mean(phi_i)
    influence_function <- phi_i - estimate
    V_ot <- mean(influence_function^2)
    variance <- V_ot / n
    
    return(list(
      estimate = as.numeric(estimate),
      variance = as.numeric(variance),
      se = as.numeric(sqrt(variance)),
      n = n,
      phi_i = phi_i,
      V_ot = V_ot,
      varphi_ot = as.numeric(influence_function),
      prop_scores_range = range(prop_scores),
      n_treated = length(treated_idx)
    ))
    
  }, error = function(e) {
    # Fallback: simple mean of treated outcomes
    est <- mean(y[treated_idx])
    centered <- y[treated_idx] - est
    theta_hat <- mean(centered^2)
    var_est <- theta_hat / length(treated_idx)
    phi_i <- rep(est, n)
    phi_i[treated_idx] <- y[treated_idx]
    return(list(
      estimate = est, variance = var_est, se = sqrt(var_est), n = n,
      phi_i = phi_i, V_ot = theta_hat,
      varphi_ot = phi_i - est,
      prop_scores_range = c(NA, NA), n_treated = length(treated_idx),
      error = paste("Estimation failed:", e$message)
    ))
  })
}
# =============================================================================
# 4. HETEROGENEITY ESTIMATION HELPER
# =============================================================================

#' DerSimonian-Laird heterogeneity estimator
#'
#' Computes Cochran's Q statistic, between-site variance τ² via DerSimonian-Laird,
#' and I² heterogeneity proportion. Used by both sample-size and inverse-variance
#' weighted estimators.
#'
#' @param site_estimates Named numeric vector of site-specific point estimates μ̂_j.
#' @param site_variances Named numeric vector of within-site variances Var_j.
#' @return List with Q (Cochran's Q), tau_sq (between-site variance),
#'   I_squared (% heterogeneity), precisions (1/Var_j), fe_estimate (fixed-effects estimate).
#' @export
calculate_dl_heterogeneity <- function(site_estimates, site_variances) {
  K <- length(site_estimates)
  precisions <- 1 / pmax(site_variances, VARIANCE_MIN)
  total_precision <- sum(precisions)
  
  # Fixed-effects weighted estimate (using inverse-variance weights)
  fe_estimate <- sum(precisions * site_estimates) / total_precision
  
  # Cochran's Q statistic: Q = Σ_j w_j (μ̂_j - μ̂_FE)²
  Q <- sum(precisions * (site_estimates - fe_estimate)^2)
  
  # DerSimonian-Laird estimator: τ² = max(0, (Q - (K-1)) / C)
  C <- total_precision - sum(precisions^2) / total_precision
  tau_sq <- 0
  if (K > 1 && C > 0) {
    tau_sq <- max(0, (Q - (K - 1)) / C)
  }
  
  # I² statistic (proportion of variance due to heterogeneity)
  I_squared <- if (Q > K - 1) (Q - (K - 1)) / Q * 100 else 0
  
  return(list(
    Q = Q,
    tau_sq = tau_sq,
    I_squared = I_squared,
    precisions = precisions,
    total_precision = total_precision,
    fe_estimate = fe_estimate
  ))
}


# =============================================================================
# THEORETICALLY CORRECT BASELINE METHODS
# =============================================================================
# The methods below use density ratio weighting to ensure ALL methods
# estimate the SAME target estimand E_t[Y(1)].
# =============================================================================

#' Select optimal density ratio lambda via cross-validation
#'
#' Uses K-fold CV on the source site to select the regularization parameter
#' for density ratio estimation. The CV criterion is the density ratio loss:
#'   L(gamma) = mean_phi_target' gamma + E_s[exp(-Z' gamma)]
#' evaluated on held-out source data.
#'
#' @param Z_source Site assignment covariates for source site (n_source x p)
#' @param Z_target Site assignment covariates for target site (n_target x p)
#' @param lambda_grid Numeric vector of candidate lambda values.
#'   Default: geometric grid from 1e-4 to 1, length 20.
#' @param n_cv_folds Number of CV folds (default 5)
#' @return Scalar: selected lambda value
select_dr_lambda_cv <- function(Z_source, Z_target,
                                lambda_grid = NULL,
                                n_cv_folds = N_CV_FOLDS_LAMBDA) {
  Z_source <- as.matrix(Z_source)
  Z_target <- as.matrix(Z_target)

  mean_phi_target <- c(1, colMeans(Z_target))
  A_dummy <- rep(1, nrow(Z_source))

  if (is.null(lambda_grid)) {
    lmax <- compute_lambda_max_initial_dr(Z_source, A_dummy, mean_phi_target, A_val = 1L)
    lambda_min_ratio <- if (nrow(Z_source) > ncol(Z_source)) 1e-4 else 0.01
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio)
  }

  # Use existing C++ CV function for initial density ratio selection
  cv_result <- tryCatch({
    select_lambda_cv_initial_density_ratio_cpp(
      Z_source, A_dummy, mean_phi_target,
      lambda_grid, n_cv_folds,
      MAX_ITER_DEFAULT, TOL_DEFAULT, 1L
    )
  }, error = function(e) {
    warning(sprintf("DR lambda CV failed: %s. Using default lambda=%.4f.",
                    conditionMessage(e), COMPARISON_DR_LAMBDA_DEFAULT))
    list(best_lambda = COMPARISON_DR_LAMBDA_DEFAULT)
  })

  return(cv_result$best_lambda)
}

#' Compute density ratio weights for source site relative to target
#' Uses exponential tilting: w(X) = exp(-Z'α)
#' @param Z_source Site assignment covariates for source site
#' @param Z_target Site assignment covariates for target site
#' @param lambda Regularization parameter
#' @return Vector of normalized density ratio weights
calculate_dr_weights <- function(Z_source, Z_target, lambda = COMPARISON_DR_LAMBDA_DEFAULT) {
  Z_source <- as.matrix(Z_source)
  Z_target <- as.matrix(Z_target)

  if (ncol(Z_source) != ncol(Z_target)) {
    stop(sprintf("calculate_dr_weights: Z_source has %d columns but Z_target has %d. They must match.",
                 ncol(Z_source), ncol(Z_target)))
  }

  n_source <- nrow(Z_source)
  
  mean_phi_target <- c(1, colMeans(Z_target))
  A_dummy <- rep(1, n_source)
  
  alpha <- tryCatch({
    fit_initial_density_ratio(Z_source, A_dummy, mean_phi_target, lambda = lambda)
  }, error = function(e) {
    warning(sprintf("calculate_dr_weights: density ratio fitting failed (%s). Using uniform weights.",
                    conditionMessage(e)))
    rep(0, length(mean_phi_target))
  })
  
  Z_int <- cbind(1, Z_source)
  eta <- as.numeric(Z_int %*% alpha)
  w <- exp(-eta)
  w <- w / mean(w)
  w <- pmax(pmin(w, DR_WEIGHT_UPPER), DR_WEIGHT_LOWER)
  
  return(w)
}

#' Compute AIPW estimate for a single site with optional density ratio weights
#'
#' Unlike \code{fit_site_aipw} (which uses cross-fitting and returns rich
#' diagnostics), this is a lightweight helper for DR-corrected methods that
#' need weighted AIPW with external density ratio weights.
#'
#' @param y Outcome vector
#' @param a Treatment vector
#' @param X Covariate matrix
#' @param weights Optional density ratio weights
#' @param family GLM family for outcome model ("binomial", "gaussian"). Default "binomial".
#' @param A_val Treatment value to estimate (default 1)
#' @return List with estimate, variance, influence function
calculate_weighted_site_aipw <- function(y, a, X, weights = NULL, family = "binomial", A_val = 1L) {
  n <- length(y)
  X <- as.matrix(X)
  if (is.null(weights)) weights <- rep(1, n)
  
  treated_idx <- which(a == A_val)
  if (length(treated_idx) < MIN_TREATED_FOR_MODEL) {
    return(list(estimate = NA, variance = Inf, psi = rep(0, n)))
  }
  
  glm_spec <- resolve_glm_family(family)
  
  # Fit propensity score via shared helper
  pi_hat <- fit_glmnet_cv(
    x_train = X, y_train = a, x_predict = X,
    family = "binomial", fallback_value = mean(a),
    clip_fn = clip_propensity,
    caller_name = "calculate_weighted_site_aipw", model_name = "PS"
  )
  
  # P(A = A_val | X)
  p_a <- if (A_val == 1L) pi_hat else (1 - pi_hat)
  p_a <- pmax(p_a, PROP_SCORE_LOWER)
  
  # Fit outcome model on A_val arm via shared helper
  X_treated <- X[treated_idx, , drop = FALSE]
  y_treated <- y[treated_idx]
  
  m_hat <- fit_glmnet_cv(
    x_train = X_treated, y_train = y_treated, x_predict = X,
    family = glm_spec$glmnet_family, fallback_value = mean(y_treated),
    clip_fn = function(pred) clip_outcome_pred(pred, family),
    caller_name = "calculate_weighted_site_aipw", model_name = "OR"
  )
  
  # Nuisance-adjusted IF-based weighted AIPW (theory-aligned with helpers used
  # by other methods). This treats density-ratio weights as fixed inputs.
  aipw_if <- calculate_aipw_influence(
    y = as.numeric(y),
    a = as.numeric(a),
    x = X,
    m_hat = m_hat,
    pi_hat = pi_hat,
    w = as.numeric(weights),
    A_val = A_val,
    family = family
  )

  return(list(
    estimate = as.numeric(aipw_if$estimate),
    variance = as.numeric(aipw_if$variance),
    psi = as.numeric(aipw_if$influence),
    varphi_ot = as.numeric(aipw_if$influence),
    phi = as.numeric(aipw_if$phi),
    weights = as.numeric(weights)
  ))
}
