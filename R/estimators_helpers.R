# estimators_helpers.R - Shared helper functions for estimators
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# The helpers in this file are arm-specific: for A_val = a they estimate
#   mu^a_t = E_t[Y(a)].
# TATE entry points pair the A_val = 1 and A_val = 0 influence blocks before
# calculating the contrast variance, so within-site cross-arm covariance is
# retained.
# =============================================================================
#
# Contents:
#   1. solve_with_ridge, calculate_aipw_pseudo_outcome, fit_logit_mle
#   2. calculate_aipw_influence
#   3. fit_site_aipw (nuisance model fitting helper)
#   4. calculate_dl_heterogeneity (DerSimonian-Laird)
#   5. calculate_dr_weights, calculate_weighted_site_aipw

#' Matrix inversion with ridge stabilization
#'
#' Inverts `mat + ridge * I` via `solve()`. If `solve()` errors, fail with
#' context instead of substituting a different linear solver, since a singular
#' adjustment matrix changes the variance/influence calculation.
#'
#' @param mat square numeric matrix
#' @param ridge non-negative ridge term to stabilize inversion
#' @return inverse of `mat + ridge * I`
solve_with_ridge <- function(mat, ridge = RIDGE_DEFAULT) {
  mat_reg <- mat + diag(ridge, nrow(mat))
  tryCatch({
    solve(mat_reg)
  }, error = function(e) {
    stop(sprintf(
      "solve_with_ridge: solve() failed on %dx%d ridge-regularized matrix (ridge=%g); refusing to switch to another solver. Increase ridge or inspect the adjustment design for collinearity. Original error: %s",
      nrow(mat_reg), ncol(mat_reg), ridge, conditionMessage(e)
    ), call. = FALSE)
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
      "gaussian" = stats::gaussian()
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
    "logit" = logistic(eta),
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

  # Response derivative h'(eta) = d mu / d eta depends on GLM family/link.
  # For canonical links, h'(eta) = V(mu) (GLM variance function).
  # Used for outcome model adjustment (influence function correction).
  m_prime <- switch(family,
    "binomial" = m_hat * (1 - m_hat),      # V(mu) = mu(1-mu) for binomial
    "gaussian" = rep(1, length(m_hat)),     # V(mu) = 1 for gaussian
    m_hat * (1 - m_hat)                     # default: binomial
  )

  # Outcome model adjustment
  s_beta <- X_int * (indicator * (y - m_hat))
  A_beta <- colMeans(w * (1 - indicator / p_a) * m_prime * X_int) / d
  M_beta <- t(X_int) %*% (X_int * (indicator * m_prime)) / n
  adj_beta <- solve_with_ridge(M_beta) %*% A_beta
  infl_beta <- as.numeric(s_beta %*% adj_beta)

  # Propensity model adjustment
  # NOTE: s_gamma_ps / adj_gamma_ps refer to the PROPENSITY SCORE model
  # coefficients, NOT the density ratio gamma_{s_j,a} used in RoCE.
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
  adj_gamma_ps <- solve_with_ridge(M_gamma_ps) %*% A_gamma_ps
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
# GRACEFUL-DEGRADATION DIAGNOSTICS (process-local)
# =============================================================================
# When a nuisance model cannot be fit on a fold (e.g. a single-class binomial
# outcome at high dimension, where every treated outcome is 1), the outcome
# model falls back to a well-defined degenerate value (the constant empirical
# mean) instead of aborting the whole Monte Carlo run. Following the convention
# of mature resampling frameworks -- tidymodels/caret record per-fold model
# failures and report them; base glm warns on separation rather than failing --
# the event must be VISIBLE, not silently swallowed. run_single_simulation()
# resets this counter and reads it back to surface a per-setting "degenerate
# fold" rate via summarize_results().
.roce_fit_diag <- new.env(parent = emptyenv())
.roce_fit_diag$or_degenerate_folds <- 0L

#' Reset the per-simulation fit-diagnostics counters.
#' @keywords internal
.reset_fit_diagnostics <- function() {
  .roce_fit_diag$or_degenerate_folds <- 0L
  invisible(NULL)
}

#' Record one outcome-model fold that fell back to the constant degenerate
#' nuisance because its binomial response was single-class.
#' @keywords internal
.record_or_degenerate_fold <- function() {
  .roce_fit_diag$or_degenerate_folds <- .roce_fit_diag$or_degenerate_folds + 1L
  invisible(NULL)
}

#' Read the current fit-diagnostics counters.
#' @return Named list with \code{or_degenerate_folds}.
#' @keywords internal
.get_fit_diagnostics <- function() {
  list(or_degenerate_folds = .roce_fit_diag$or_degenerate_folds)
}

# =============================================================================
# SHARED GLMNET CV HELPER
# =============================================================================

#' Fit L1-penalized GLM via glmnet cross-validation and predict on new data
#'
#' Shared helper that encapsulates the repeated pattern:
#'   \code{cv.glmnet(x_train, y_train) -> predict(x_predict) -> clip}
#' with automatic fold count selection and contextual errors.
#'
#' @param x_train Training covariate matrix (n_train x p).
#' @param y_train Training response vector (length n_train).
#' @param x_predict Prediction covariate matrix (n_predict x p). Columns must
#'   match \code{x_train}; number of rows may differ.
#' @param family glmnet family string (\code{"binomial"}, \code{"gaussian"}, etc.).
#' @param clip_fn Optional unary function applied element-wise to the
#'   predictions (e.g., \code{clip_propensity}, \code{clip_outcome_pred}).
#' @param nlambda Number of lambda grid points for \code{cv.glmnet}.
#' @param min_per_fold Minimum observations per fold passed to
#'   \code{get_cv_fold_count} (default 10).
#' @param caller_name Caller function name for warning messages.
#' @param model_name Short model label for warning messages
#'   (e.g., \code{"PS"}, \code{"OR"}).
#' @param fold_id Optional integer fold index included in the warning.
#' @param on_degenerate_response How to handle a single-class binomial response
#'   (a class with fewer than 2 observations), which \code{glmnet} cannot fit.
#'   \code{"error"} (default) fails loudly -- appropriate for the propensity
#'   model, where it signals a positivity failure. \code{"constant"} returns the
#'   constant empirical mean -- the correct degenerate nuisance for the outcome
#'   model when, e.g., a saturated high-dimensional fold has all-1 outcomes.
#' @param lambda_rule Cross-validation rule: \code{"min"} uses
#'   \code{lambda.min}; \code{"1se"} uses \code{lambda.1se}.
#' @param cv_group_id Optional positive integer origin/group identifier aligned
#'   with the training rows. Repeated IDs are kept intact within nuisance-CV
#'   folds. \code{NULL}, or an all-unique vector, preserves the legacy path.
#' @return Numeric vector of predictions (length \code{nrow(x_predict)}). The
#'   integer attribute \code{outcome_degenerate} is set to 1 when the
#'   constant-outcome fallback was used and is otherwise absent.
#' @keywords internal
fit_glmnet_cv <- function(x_train, y_train, x_predict,
                          family = "binomial",
                          clip_fn = NULL,
                          nlambda = LAMBDA_GRID_SIZE_STANDARD,
                          min_per_fold = 10L,
                          caller_name = "", model_name = "",
                          fold_id = NULL,
                          on_degenerate_response = c("error", "constant"),
                          lambda_rule = c("min", "1se"),
                          cv_group_id = NULL) {
  on_degenerate_response <- match.arg(on_degenerate_response)
  lambda_rule <- .match_nuisance_lambda_rule(
    lambda_rule, "fit_glmnet_cv"
  )
  n_predict <- nrow(x_predict)
  cv_group_caller <- paste0(
    if (nzchar(caller_name)) caller_name else "fit_glmnet_cv",
    if (nzchar(model_name)) paste0(" ", model_name) else ""
  )
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, nrow(x_train), cv_group_caller
  )

  # Saturated/near-separable binomial outcome folds: below glmnet's own
  # <8-per-class threshold the penalized logistic CV is unstable (single-class
  # sub-folds, or a non-conformable failure from inconsistent per-fold lambda
  # paths on near-separable data). For the outcome model the constant empirical
  # mean is the correct degenerate nuisance. The propensity model keeps its
  # strict behaviour -- treatment is ~balanced, so this never triggers there.
  if (family == "binomial" && on_degenerate_response == "constant" &&
      min(table(factor(y_train, levels = c(0, 1)))) < 8L) {
    .record_or_degenerate_fold()
    pred <- rep(mean(y_train), n_predict)
    if (!is.null(clip_fn)) pred <- clip_fn(pred)
    attr(pred, "outcome_degenerate") <- 1L
    return(pred)
  }

  pred <- tryCatch({
    n_cv_folds <- get_cv_fold_count(nrow(x_train), min_per_fold = min_per_fold)
    if (n_cv_folds < 4L) {
      stop(sprintf("insufficient sample size for cv.glmnet: n_train=%d yields nfolds=%d (<4 minimum valid folds)",
                   nrow(x_train), n_cv_folds))
    }
    nuisance_fold_id <- .make_nuisance_cv_fold_id(
      cv_group_id, n_cv_folds, cv_group_caller
    )
    cv_fit <- if (is.null(nuisance_fold_id)) {
      glmnet::cv.glmnet(
        x = x_train, y = y_train, family = family,
        alpha = 1, nfolds = n_cv_folds, nlambda = nlambda,
        maxit = GLMNET_MAX_ITER
      )
    } else {
      glmnet::cv.glmnet(
        x = x_train, y = y_train, family = family,
        alpha = 1, foldid = nuisance_fold_id, nlambda = nlambda,
        maxit = GLMNET_MAX_ITER
      )
    }
    selected_lambda <- if (identical(lambda_rule, "1se")) {
      "lambda.1se"
    } else {
      "lambda.min"
    }
    as.numeric(predict(
      cv_fit, newx = x_predict, s = selected_lambda, type = "response"
    ))
  }, error = function(e) {
    msg <- conditionMessage(e)
    # A (sub-)fold can carry a single-class binomial response -- e.g. a saturated
    # high-dimensional outcome where all treated outcomes are 1 -- which glmnet
    # refuses to fit ("one ... class has 1 or 0 observations"). For the outcome
    # model the correct degenerate nuisance is the constant empirical mean; the
    # propensity model keeps the default "error" (a positivity failure).
    if (on_degenerate_response == "constant" &&
        grepl("1 or 0 observations|0 or 1 observations|non-conformable", msg)) {
      .record_or_degenerate_fold()
      fallback <- rep(mean(y_train), n_predict)
      attr(fallback, "outcome_degenerate") <- 1L
      return(fallback)
    }
    fold_msg <- if (!is.null(fold_id)) sprintf(" on fold %d", fold_id) else ""
    stop(sprintf("%s: %s model failed%s: %s",
                 caller_name, model_name, fold_msg, msg),
         call. = FALSE)
  })

  used_degenerate_fallback <- identical(
    attr(pred, "outcome_degenerate"), 1L
  )
  if (!is.null(clip_fn)) pred <- clip_fn(pred)
  if (used_degenerate_fallback) {
    attr(pred, "outcome_degenerate") <- 1L
  }
  if (length(pred) != n_predict || any(!is.finite(pred))) {
    stop(sprintf("%s: %s model produced invalid predictions (expected length %d, got %d; non-finite=%d).",
                 caller_name, model_name, n_predict, length(pred),
                 sum(!is.finite(pred))), call. = FALSE)
  }
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
#' @param site_data List with components Y, A, W_outcome, optional Z_site, n
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Whether to use RCAL (default FALSE, uses glmnet)
#' @param use_crossfit Whether to use cross-fitting (default TRUE)
#' @param n_folds Number of cross-fitting folds (NULL = data-driven)
#' @param A_val Treatment value to estimate potential outcome for (default 1)
#' @return List with estimate, variance, se, phi_i, V_ot, varphi_ot, n,
#'         prop_scores_range, n_treated. Throws an error on estimation failure.
fit_site_aipw <- function(site_data, family = "binomial", use_rcal = FALSE,
                          use_crossfit = TRUE, n_folds = NULL, A_val = 1L) {
  glm_spec <- resolve_glm_family(family)
  y <- site_data$Y
  tr <- site_data$A
  x_or <- as.matrix(site_data$W_outcome)
  x_ps <- if (!is.null(site_data$Z_site)) {
    as.matrix(site_data$Z_site)
  } else {
    x_or
  }
  n <- site_data$n
  
  treated_idx <- which(tr == A_val)
  if (length(treated_idx) == 0) {
    stop(sprintf("fit_site_aipw: no observations with A_val=%d; AIPW estimator is not identifiable for this site.",
                 A_val), call. = FALSE)
  }
  
  # Try cross-fitting first (preferred for valid inference)
  if (isTRUE(use_crossfit)) {
    n_cv_folds <- if (is.null(n_folds)) get_cv_fold_count(n) else n_folds
    cf_res <- estimate_target_only_crossfit(site_data, n_folds = n_cv_folds,
                                            family = family, A_val = A_val)
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
  
  # Non-cross-fitted estimation, used only when explicitly requested.
    X_or_matrix <- as.matrix(x_or)
    X_ps_matrix <- as.matrix(x_ps)
    if (any(is.na(X_or_matrix)) || any(is.na(X_ps_matrix)) ||
        any(is.na(tr)) || any(is.na(y))) {
      stop(sprintf("fit_site_aipw: Missing values detected: W_outcome has %d NA(s), Z_site has %d NA(s), A has %d NA(s), Y has %d NA(s). Remove or impute before calling.",
                   sum(is.na(X_or_matrix)), sum(is.na(X_ps_matrix)),
                   sum(is.na(tr)), sum(is.na(y))))
    }
    n_cv_folds <- get_cv_fold_count(n)
    
    # Step 1: Fit propensity score model P(A=1|Z)
    prop_scores <- NULL
    if (use_rcal) {
      if (!requireNamespace("RCAL", quietly = TRUE)) {
        stop("fit_site_aipw: use_rcal=TRUE but package 'RCAL' is not available.",
             call. = FALSE)
      }
      prop_scores <- tryCatch({
        ps_result <- RCAL::glm.regu.cv(
          fold = n_cv_folds, y = as.numeric(tr), x = X_ps_matrix,
          loss = "cal", nrho = LAMBDA_GRID_SIZE_STANDARD
        )
        if (!is.null(ps_result$sel.fit) && !any(is.na(ps_result$sel.fit[, 1]))) {
          ps_result$sel.fit[, 1]
        } else {
          stop("RCAL propensity fit did not return finite selected fitted values.")
        }
      }, error = function(e) {
        stop(sprintf("fit_site_aipw: RCAL propensity score model failed: %s",
                     conditionMessage(e)), call. = FALSE)
      })
    }
    if (is.null(prop_scores)) {
      cv_fit <- glmnet::cv.glmnet(
        x = X_ps_matrix, y = as.numeric(tr), family = "binomial",
        alpha = 1, nfolds = n_cv_folds,
        nlambda = LAMBDA_GRID_SIZE_STANDARD,
        maxit = GLMNET_MAX_ITER
      )
      prop_scores <- as.numeric(predict(cv_fit, newx = X_ps_matrix,
                                        s = "lambda.min", type = "response"))
    }
    prop_scores <- clip_propensity(prop_scores)
    
    # Step 2: Fit outcome regression model E[Y|A=a, W]
    X_treated <- X_or_matrix[treated_idx, , drop = FALSE]
    y_treated <- as.numeric(y[treated_idx])
    n_cv_folds <- get_cv_fold_count(length(y_treated))
    
    m1_pred <- NULL
    if (use_rcal) {
      if (!requireNamespace("RCAL", quietly = TRUE)) {
        stop("fit_site_aipw: use_rcal=TRUE but package 'RCAL' is not available.",
             call. = FALSE)
      }
      or_result <- tryCatch({
        loss_type <- switch(family,
          "binomial" = "ml",
          "gaussian" = "gaus",
          "gaus"
        )
        RCAL::glm.regu.cv(
          fold = n_cv_folds, y = y_treated, x = X_treated,
          loss = loss_type, nrho = LAMBDA_GRID_SIZE_STANDARD
        )
      }, error = function(e) {
        stop(sprintf("fit_site_aipw: RCAL outcome regression failed: %s",
                     conditionMessage(e)), call. = FALSE)
      })
      if (!is.null(or_result) && !is.null(or_result$sel.bet)) {
        beta_coef <- or_result$sel.bet[, 1]
        X_design <- cbind(1, X_or_matrix)
        if (ncol(X_design) == length(beta_coef)) {
          m1_pred <- as.numeric(X_design %*% beta_coef)
        } else {
          stop(sprintf("fit_site_aipw: RCAL outcome coefficient length %d does not match design columns %d.",
                       length(beta_coef), ncol(X_design)), call. = FALSE)
        }
      } else {
        stop("fit_site_aipw: RCAL outcome regression did not return selected coefficients.",
             call. = FALSE)
      }
    }
    if (is.null(m1_pred)) {
      family_type <- glm_spec$glmnet_family
      cv_fit_or <- glmnet::cv.glmnet(
        x = X_treated, y = y_treated, family = family_type,
        alpha = 1, nfolds = n_cv_folds,
        nlambda = LAMBDA_GRID_SIZE_STANDARD,
        maxit = GLMNET_MAX_ITER
      )
      beta_coef <- as.vector(coef(cv_fit_or, s = "lambda.min"))
      X_design <- cbind(1, X_or_matrix)
      m1_pred <- as.numeric(X_design %*% beta_coef)
    }
    # Apply response function: linear predictor to response scale.
    # RCAL and manual coef() extraction return linear predictors
    # (eta = X^T beta), but AIPW needs response-scale predictions.
    # Note: glmnet predict(..., type="response") already applies this transform,
    # but the manual extraction path above does not.
    if (family == "binomial") {
      m1_pred <- logistic(m1_pred)
    }
    # gaussian/identity: no transform needed (mu = eta)
    m1_pred <- clip_outcome_pred(m1_pred, family)
    
    # Step 3: Compute doubly robust AIPW pseudo-outcomes
    p_a <- if (A_val == 1) prop_scores else (1 - prop_scores)
    phi_i <- calculate_aipw_pseudo_outcome(as.numeric(y), tr, m1_pred, p_a, A_val = A_val)
    
    if (any(is.na(phi_i)) || any(is.infinite(phi_i))) {
      stop(sprintf("fit_site_aipw: AIPW pseudo-outcome produced %d non-finite value(s).",
                   sum(is.na(phi_i) | is.infinite(phi_i))), call. = FALSE)
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
}
# =============================================================================
# 4. HETEROGENEITY ESTIMATION HELPER
# =============================================================================

#' DerSimonian-Laird heterogeneity estimator
#'
#' Computes Cochran's Q statistic and the DerSimonian-Laird between-site
#' variance and heterogeneity proportion. Used by both sample-size and inverse-variance
#' weighted estimators.
#'
#' @param site_estimates Named numeric vector of site-specific point estimates.
#' @param site_variances Named numeric vector of within-site variances Var_j.
#' @return List with Q (Cochran's Q), tau_sq (between-site variance),
#'   I_squared (\% heterogeneity), precisions (1/Var_j), fe_estimate (fixed-effects estimate).
#' @export
calculate_dl_heterogeneity <- function(site_estimates, site_variances) {
  K <- length(site_estimates)
  precisions <- 1 / pmax(site_variances, VARIANCE_MIN)
  total_precision <- sum(precisions)
  
  # Fixed-effects weighted estimate (using inverse-variance weights)
  fe_estimate <- sum(precisions * site_estimates) / total_precision
  
  # Cochran's Q statistic.
  Q <- sum(precisions * (site_estimates - fe_estimate)^2)
  
  # DerSimonian-Laird estimator: tau2 = max(0, (Q - (K-1)) / C)
  C <- total_precision - sum(precisions^2) / total_precision
  tau_sq <- 0
  if (K > 1 && C > 0) {
    tau_sq <- max(0, (Q - (K - 1)) / C)
  }
  
  # I-squared statistic (proportion of variance due to heterogeneity)
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
# 5. BOOTSTRAP VARIANCE HELPER (comparison-method standard errors)
# =============================================================================

#' Multiplier (wild) bootstrap standard error for a weighted sum of influence blocks
#'
#' Each comparison estimator can be written as \eqn{\sum_b w_b \, \tildeE_b[\varphi_b]},
#' a weighted sum of block means of (mean-zero) per-observation influence functions
#' \eqn{\varphi_b}, where blocks are mutually independent samples (e.g. one block per
#' site, or a shared-target block plus per-source blocks). This computes the standard
#' error by a wild bootstrap: for each replicate it draws Rademacher multipliers
#' \eqn{\omega_i \in \{-1,+1\}} per observation and forms
#' \eqn{\sum_b w_b \, \overline{\omega \, \varphi_b}}; the SE is the standard deviation
#' across replicates. In expectation this equals the fixed-effects influence-function
#' variance \eqn{\sum_b w_b^2 \, \overline{\varphi_b^2} / n_b}, i.e. an honest sampling
#' variance that, unlike a random-effects (\eqn{\tau^2}) correction, does NOT inflate
#' for systematic between-block (covariate-shift) differences.
#'
#' The global RNG state is saved and restored so the bootstrap draws do not perturb
#' any downstream random number generation.
#'
#' @param blocks List of blocks; each block is a list with `influence` (numeric vector,
#'   centered internally) and `weight` (scalar combination weight). Empty-influence
#'   blocks are dropped.
#' @param n_bootstrap Number of bootstrap replicates.
#' @return Bootstrap standard error (numeric scalar), or \code{NA_real_} if no block has
#'   any observation.
#' @keywords internal
.validate_bootstrap_replicates <- function(n_bootstrap, caller) {
  if (length(n_bootstrap) != 1L || is.na(n_bootstrap) ||
      !is.finite(n_bootstrap) || n_bootstrap < 2 ||
      n_bootstrap %% 1 != 0 || n_bootstrap > .Machine$integer.max) {
    stop(
      sprintf("%s: n_bootstrap must be one integer >= 2.", caller),
      call. = FALSE
    )
  }
  as.integer(n_bootstrap)
}

.multiplier_bootstrap_se <- function(blocks, n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, ".multiplier_bootstrap_se"
  )
  blocks <- Filter(function(b) length(b$influence) > 0L, blocks)
  if (length(blocks) == 0L) {
    return(NA_real_)
  }

  centered_influence <- lapply(blocks, function(b) {
    infl <- as.numeric(b$influence)
    infl - mean(infl)
  })
  weights <- vapply(blocks, function(b) as.numeric(b$weight), numeric(1L))

  # Keep the bootstrap draws from perturbing the global RNG stream, whether or not
  # the RNG had already been initialized before this call.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    saved_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", saved_seed, envir = globalenv()), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      },
      add = TRUE
    )
  }

  replicates <- vapply(seq_len(n_bootstrap), function(.rep) {
    total <- 0
    for (j in seq_along(centered_influence)) {
      infl <- centered_influence[[j]]
      multipliers <- sample(c(-1, 1), length(infl), replace = TRUE)
      total <- total + weights[[j]] * mean(multipliers * infl)
    }
    total
  }, numeric(1L))

  stats::sd(replicates)
}

#' Resolve the reported variance of a comparison estimator (bootstrap or analytic)
#'
#' Centralises the bootstrap-versus-analytic choice so every comparison method reports
#' variance consistently. When \code{variance_method = "bootstrap"} (the default), the
#' reported variance is the wild-bootstrap SE squared (see
#' \code{.multiplier_bootstrap_se});
#' the analytic variance is always retained alongside for reference. Falls back to the
#' analytic variance if the bootstrap cannot be computed (no usable blocks).
#'
#' @param analytic_variance Method-specific analytic variance estimate.
#' @param blocks Influence-function blocks for the bootstrap (see
#'   \code{.multiplier_bootstrap_se}).
#' @param variance_method Either \code{"bootstrap"} (default) or \code{"analytic"}.
#' @param n_bootstrap Number of bootstrap replicates.
#' @return List with `variance`, `se`, `variance_method` (the method actually used),
#'   `variance_analytic`, `se_analytic`, `se_bootstrap`, and `n_bootstrap`.
#' @keywords internal
.resolve_comparison_variance <- function(analytic_variance, blocks,
                                         variance_method = c("bootstrap", "analytic"),
                                         n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  variance_method <- match.arg(variance_method)
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, ".resolve_comparison_variance"
  )
  analytic_variance <- max(as.numeric(analytic_variance), VARIANCE_MIN)

  se_bootstrap <- if (identical(variance_method, "bootstrap")) {
    .multiplier_bootstrap_se(blocks, n_bootstrap = n_bootstrap)
  } else {
    NA_real_
  }

  # A non-finite or non-positive bootstrap SE (e.g. degenerate all-constant influence)
  # is treated as undefined: fall back to the analytic variance rather than silently
  # collapsing to VARIANCE_MIN while still claiming method = "bootstrap".
  use_bootstrap <- identical(variance_method, "bootstrap") &&
    is.finite(se_bootstrap) && se_bootstrap > 0
  variance <- if (use_bootstrap) max(se_bootstrap^2, VARIANCE_MIN) else analytic_variance

  list(
    variance = variance,
    se = sqrt(variance),
    variance_method = if (use_bootstrap) "bootstrap" else "analytic",
    variance_analytic = analytic_variance,
    se_analytic = sqrt(analytic_variance),
    se_bootstrap = se_bootstrap,
    n_bootstrap = n_bootstrap
  )
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
#'   Default: glmnet-style geometric path from lambda_max down to
#'   1e-4 * lambda_max.
#' @param n_cv_folds Number of CV folds. If \code{NULL}, an adaptive fold count
#'   is chosen from the source sample size.
#' @param lambda_rule CV selection rule: \code{"min"} uses \code{lambda.min};
#'   \code{"1se"} uses \code{lambda.1se}.
#' @return Scalar: selected lambda value
select_dr_lambda_cv <- function(Z_source, Z_target,
                                lambda_grid = NULL,
                                n_cv_folds = NULL,
                                lambda_rule = c("min", "1se")) {
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "select_dr_lambda_cv")
  Z_source <- as.matrix(Z_source)
  Z_target <- as.matrix(Z_target)
  if (nrow(Z_source) < 2L) {
    stop(sprintf("select_dr_lambda_cv: at least two source observations are required for CV; found %d.",
                 nrow(Z_source)), call. = FALSE)
  }
  if (is.null(n_cv_folds)) {
    n_cv_folds <- get_cv_fold_count(nrow(Z_source))
  } else if (length(n_cv_folds) != 1L || is.na(n_cv_folds) || n_cv_folds < 2L) {
    stop("select_dr_lambda_cv: n_cv_folds must be NULL or an integer >= 2.",
         call. = FALSE)
  } else {
    n_cv_folds <- as.integer(min(n_cv_folds, nrow(Z_source)))
  }

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
      MAX_ITER_DEFAULT, TOL_DEFAULT, 1L, M_tau = Inf
    )
  }, error = function(e) {
    stop(sprintf("select_dr_lambda_cv: DR lambda CV failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })

  .select_nuisance_cv_lambda(cv_result, lambda_rule, "select_dr_lambda_cv")
}

#' Compute density ratio weights for source site relative to target
#' Uses exponential tilting: w(X) = exp(-Z' alpha)
#' @param Z_source Site assignment covariates for source site
#' @param Z_target Site assignment covariates for target site
#' @param lambda Regularization parameter. \code{NULL} (default) selects
#'   \code{lambda.min} by CV; numeric values use a fixed lambda.
#' @param lambda_rule CV selection rule when \code{lambda = NULL}:
#'   \code{"min"} selects \code{lambda.min}; \code{"1se"} selects
#'   \code{lambda.1se}.
#' @return Vector of normalized and bounded density-ratio weights. Attributes
#'   retain the selected penalty and pre-clipping diagnostics for downstream
#'   overlap audits.
calculate_dr_weights <- function(Z_source, Z_target, lambda = NULL,
                                 lambda_rule = c("min", "1se")) {
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "calculate_dr_weights")
  Z_source <- as.matrix(Z_source)
  Z_target <- as.matrix(Z_target)

  if (ncol(Z_source) != ncol(Z_target)) {
    stop(sprintf("calculate_dr_weights: Z_source has %d columns but Z_target has %d. They must match.",
                 ncol(Z_source), ncol(Z_target)))
  }

  n_source <- nrow(Z_source)
  if (n_source < 2L) {
    stop(sprintf("calculate_dr_weights: at least two source observations are required; found %d.",
                 n_source), call. = FALSE)
  }
  if (nrow(Z_target) == 0L) {
    stop("calculate_dr_weights: target site has 0 observations.", call. = FALSE)
  }
  
  mean_phi_target <- c(1, colMeans(Z_target))
  A_dummy <- rep(1, n_source)
  
  alpha <- tryCatch({
    fit_initial_density_ratio(
      Z_source, A_dummy, mean_phi_target,
      lambda = lambda, M_tau = Inf, lambda_rule = lambda_rule
    )
  }, error = function(e) {
    stop(sprintf("calculate_dr_weights: density ratio fitting failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })
  
  Z_int <- cbind(1, Z_source)
  eta <- as.numeric(Z_int %*% alpha)
  w <- .normalize_log_weights(-eta, "calculate_dr_weights")
  n_below <- sum(w < DR_WEIGHT_LOWER)
  n_above <- sum(w > DR_WEIGHT_UPPER)
  preclip_min <- min(w)
  preclip_max <- max(w)
  w <- pmax(pmin(w, DR_WEIGHT_UPPER), DR_WEIGHT_LOWER)
  if (n_below + n_above > 0L) {
    warning(sprintf(
      "calculate_dr_weights: clipped %d of %d density-ratio weights to [%g, %g] (below: %d, above: %d). A large clipped fraction signals an unreliable density-ratio model (covariate-shift extrapolation).",
      n_below + n_above, length(w), DR_WEIGHT_LOWER, DR_WEIGHT_UPPER,
      n_below, n_above
    ), call. = FALSE)
  }
  attr(w, "lambda_used") <- attr(alpha, "lambda_used") %||% as.numeric(lambda)
  attr(w, "lambda_rule") <- attr(alpha, "lambda_rule") %||% "fixed"
  attr(w, "clipping_diagnostics") <- list(
    n = length(w),
    n_below = n_below,
    n_above = n_above,
    n_clipped = n_below + n_above,
    fraction_clipped = (n_below + n_above) / length(w),
    preclip_min = preclip_min,
    preclip_max = preclip_max
  )

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
    warning(sprintf(
      "calculate_weighted_site_aipw: only %d unit(s) with A == %s out of n=%d (need >= %d for nuisance fitting). Returning NA estimate, Inf variance, and zero influence function.",
      length(treated_idx), format(A_val), n, MIN_TREATED_FOR_MODEL
    ), call. = FALSE)
    return(list(estimate = NA_real_, variance = Inf, psi = rep(0, n)))
  }
  
  glm_spec <- resolve_glm_family(family)
  
  # Fit propensity score via shared helper
  pi_hat <- fit_glmnet_cv(
    x_train = X, y_train = a, x_predict = X,
    family = "binomial",
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
    family = glm_spec$glmnet_family,
    clip_fn = function(pred) clip_outcome_pred(pred, family),
    caller_name = "calculate_weighted_site_aipw", model_name = "OR",
    on_degenerate_response = "constant"
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
