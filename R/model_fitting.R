# model_fitting.R - Unified model fitting functions for FACE algorithm
#
# This file provides both R-based (glmnet) and C++ accelerated implementations
# of model fitting functions for causal inference in federated settings.
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# These functions support estimation of the POTENTIAL OUTCOME MEAN:
#   μ^a_t = E_t[Y(a)]  where a ∈ {0, 1}
#
# By default, we estimate μ¹_t = E_t[Y(1)] (treatment effect for treated).
# Set A_val = 0 to estimate μ⁰_t = E_t[Y(0)] (outcome for control).
# ATE = μ¹_t - μ⁰_t requires running the algorithm twice.
#
# ============================================================================
# ESTIMAND: μ¹ = E_t[Y(1)] - Potential Outcome Mean (NOT ATE)
# ============================================================================
# The FACE algorithm estimates the potential outcome mean for the treated group
# on the target site population: μ¹ = E[Y(1) | R = t]
#
# To estimate ATE = E[Y(1) - Y(0)], run with A_val=1 and A_val=0 separately.
#
# ============================================================================
# NOTATION MAPPING (Code Variable → Paper Symbol in main.tex)
# ============================================================================
# - alpha / alpha_init → α (outcome model parameters)
# - gamma / gamma_s → γ (site/treatment model / density ratio parameters)
# - alpha_ot → α_{ot} (target-only outcome model parameters)
# - X / W_outcome → φ(X) (basis expansion of covariates for outcome model)
# - Z_site → φ(X) (basis expansion for site/treatment model)
# - A_val → a ∈ {0,1} (treatment indicator value)
# - delta_ts → δ_{t,s_j} (correction term from source site)
# - mu_pred_ts → Ẽ_t[ψ(φ(X); α_{t,s_j})] (outcome model prediction on target)
#
# ============================================================================
# FUNCTION CATEGORIES
# ============================================================================
# 1. R-based GLM fitting (using glmnet):
#    - fit_initial_outcome(): Initial outcome model with L1 regularization
#
# 2. C++ accelerated fitting (coordinate descent):
#    - fit_initial_density_ratio(): Initial density ratio (eq:gamma_init in main.tex)
#    - fit_unified_density_ratio(): Refined/Calibrated density ratio
#    - fit_unified_outcome(): Refined/Calibrated outcome model
#
# 3. Variance and aggregation functions:
#    - optimize_weights(): Optimal aggregation weights (eq:final_opt in main.tex)
#    - calculate_aggregated_variance(): Aggregated variance with cross-site cov
#
# ============================================================================

# Dependencies: constants.R, numerical_utils.R
# All loaded automatically by the R package system (see DESCRIPTION Collate field).

# ============================================================================
# R-BASED GLM FITTING (using glmnet)
# ============================================================================

#' Fit initial outcome model using L1-regularized GLM
#'
#' Implements the initial outcome model estimation for the two-level cross-fitting
#' algorithm. Following eq:alpha_init in main.tex:
#'
#' $$
#' \ell(\boldsymbol{\alpha}_{t,s_j}^{1}) = \widetilde{E}_{s_j}\left[I(A=1)
#'   \ell_{\text{GLM}}\left(Y, \psi(\phi(\mathbf{X}); \boldsymbol{\alpha}_{t,s_j}^1)\right)\right]
#'   + \lambda_\alpha \|\boldsymbol{\alpha}_{t,s_j}^{1}\|_1
#' $$
#'
#' This function uses `glmnet` for efficient L1-regularized GLM fitting.
#' Lambda is selected via cross-validation if not provided.
#'
#' @param W_outcome Covariate matrix (n x p), the design matrix phi(X)
#' @param Y Outcome vector (n x 1)
#' @param A Treatment indicator vector (n x 1), binary \{0, 1\}
#' @param A_val Treatment value to fit (default 1 for treated group)
#' @param lambda L1 regularization parameter. If NULL, selected via CV.
#' @param nlambda Number of lambda values in glmnet path (default 100).
#'        Lower values (e.g., 20) speed up CV with minimal precision loss.
#' @param family GLM family: "gaussian" or "binomial".
#'        Default "binomial".
#'
#' @return Vector of outcome model parameters alpha (including intercept as first element)
#'
#' @details
#' The function:
#' 1. Filters data to the specified treatment group (A = A_val)
#' 2. Uses cv.glmnet for lambda selection if lambda is NULL
#' 3. Fits a GLM with L1 penalty using the specified family
#' 4. Returns coefficients including the intercept
#'
#' @seealso \code{\link{fit_unified_outcome}} for C++ accelerated outcome fitting
#' @export
fit_initial_outcome <- function(W_outcome, Y, A, A_val = 1L, lambda = NULL, 
                                nlambda = LAMBDA_GRID_SIZE_STANDARD, family = "binomial") {
  
  # Resolve GLM family for glmnet
  glm_spec <- resolve_glm_family(family)
  glmnet_family <- glm_spec$glmnet_family
  
  # Filter to treatment group
  treated_idx <- which(A == A_val)
  if (length(treated_idx) == 0) {
    warning(sprintf("fit_initial_outcome: No observations with A=%d (n=%d). Returning zero coefficients.",
                    A_val, nrow(W_outcome)))
    return(rep(0, ncol(W_outcome) + 1))
  }
  
  X_treated <- W_outcome[treated_idx, , drop = FALSE]
  Y_treated <- Y[treated_idx]
  
  # Ensure minimum sample size for CV
  n_treated <- length(Y_treated)
  
  # ---------- Standard path: glmnet L1-regularized GLM -----------------------
  if (is.null(lambda)) {
    # Use cross-validation to select lambda
    lambda_use <- tryCatch({
      # Adaptive number of folds based on sample size
      n_cv_folds <- get_cv_fold_count(n_treated)
      
      cv_fit <- glmnet::cv.glmnet(
        x = X_treated,
        y = Y_treated,
        family = glmnet_family,
        alpha = 1,  # Lasso penalty
        standardize = TRUE,
        nfolds = n_cv_folds,
        nlambda = nlambda  # Control lambda path length
      )
      # Use lambda.min (minimises CV error) instead of lambda.1se
      # (1-SE rule).  For causal inference the priority is consistency
      # (bias → 0) rather than prediction parsimony.  lambda.1se
      # over-shrinks, especially on nonlinear X_dagger features,
      # breaking the doubly-robust bias cancellation.
      cv_fit$lambda.min
    }, error = function(e) {
      # Fallback lambda if CV fails
      warning(sprintf("cv.glmnet failed (family=%s, n_treated=%d, nlambda=%d): %s. Using LAMBDA_DEFAULT=%.4f.",
                      glmnet_family, n_treated, nlambda, conditionMessage(e), LAMBDA_DEFAULT))
      LAMBDA_DEFAULT
    })
  } else {
    lambda_use <- lambda
  }
  
  # Fit final model with selected lambda
  fit <- glmnet(
    x = X_treated,
    y = Y_treated,
    family = glmnet_family,
    alpha = 1,
    lambda = lambda_use,
    standardize = TRUE
  )
  
  # Extract coefficients (includes intercept as first element)
  alpha <- as.vector(coef(fit, s = lambda_use))
  
  # Return both alpha and the selected lambda (for caching)
  attr(alpha, "lambda_used") <- lambda_use
  
  return(alpha)
}

#' Calculate true potential outcome mean for binary outcomes
#'
#' Computes the true potential outcome mean E[Y(a)] for a given treatment
#' level using the logistic model structure from the data generation process.
#'
#' @param site_data Data for the site containing W_outcome, A, Y, n
#' @param a Treatment level (0 or 1)
#'
#' @return Scalar: estimated E[Y(a)] for the site
#'
#' @details
#' For binary outcomes with logistic model:
#' E[Y(a)] = E[ψ(φ(X)^T α_a)]
# ============================================================================
# C++ ACCELERATED FUNCTIONS
# ============================================================================

# C++ functions (fit_unified_density_ratio_cpp, etc.) are compiled and loaded
# automatically by the package build system via useDynLib(FACEC) in NAMESPACE.
# The src/ C++ files are compiled during R CMD INSTALL or devtools::load_all().

#' Fit INITIAL density ratio function (γ_init parameters)
#'
#' Implements the INITIAL loss function for density ratio estimation.
#' Following eq:gamma_init in main.tex:
#' $$
#' \ell(\boldsymbol{\gamma}_{s_j,1}) = \left(\widetilde{E}_t[\boldsymbol{\phi}(\mathbf{X})]\right)^T \boldsymbol{\gamma}_{s_j,1}
#'   + \widetilde{E}_{s_j}\left[I(A=1) \exp(-\boldsymbol{\phi}(\mathbf{X})^T \boldsymbol{\gamma}_{s_j,1})\right]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}_{s_j,1}\|_1
#' $$
#'
#' NOTE: This does NOT include ψ'(α_init) term - only uses mean_phi from target
#'
#' @param Z_site Covariate matrix (n × p)
#' @param A Treatment indicator vector (n × 1)
#' @param mean_phi Mean of φ(X) from target site: Ẽ_t[φ(X)]
#' @param lambda L1 regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Vector of initial density ratio parameters
#' @export
fit_initial_density_ratio <- function(Z_site, A, mean_phi, lambda = LAMBDA_DEFAULT, 
                                       max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                       A_val = 1L, warm_start = NULL) {
  # Handle lambda = NULL case: use CV selection
  if (is.null(lambda)) {
    lmax <- compute_lambda_max_initial_dr(Z_site, A, mean_phi, A_val)
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) 1e-4 else 0.01
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio)
    cv_result <- select_lambda_cv_initial_density_ratio_cpp(
      Z_site, A, mean_phi, lambda_grid, N_CV_FOLDS_LAMBDA, max_iter, tol, A_val
    )
    lambda <- cv_result$best_lambda
  }
  
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  cpp_result <- fit_initial_density_ratio_cpp(Z_site, A, mean_phi, lambda, max_iter, tol, A_val, ws)
  result <- cpp_result$gamma
  attr(result, "lambda_used") <- lambda
  return(result)
}

#' Fit UNIFIED density ratio function (γ parameters)
#'
#' Implements the density ratio estimation loss function.
#' 
#' **Refined Loss (untruncated counterpart of eq:gamma_calibrated_loss, calibrated = FALSE):**
#' $$
#' \ell(\boldsymbol{\gamma}) = \widetilde{E}_t[\nabla_\alpha \psi]^T \boldsymbol{\gamma}
#'   + \widetilde{E}_{s_j}[I(A=1) \exp(-\phi_{site}^T \boldsymbol{\gamma}) \psi'(\phi_{outcome}^T \boldsymbol{\alpha})]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}\|_1
#' $$
#'
#' **Calibrated Loss (eq:gamma_calibrated_loss in main.tex, calibrated = TRUE):**
#' $$
#' \ell_{\gamma,\text{cal}}(\boldsymbol{\gamma}) = \widetilde{E}_t[\nabla_\alpha \psi]^T \boldsymbol{\gamma}
#'   + \widetilde{E}_{s_j}[I(A=1) \exp(-\phi_{site}^T \boldsymbol{\gamma}) \psi'(\mathcal{T}(\phi_{outcome}^T \boldsymbol{\alpha}))]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}\|_1
#' $$
#'
#' @param Z_site Site model features for gamma parameterization (n x p_site)
#' @param A Treatment indicator vector (n x 1)
#' @param mean_grad_psi Mean gradient of ψ from target site, i.e., Ẽ_t[∇_α ψ(ϕ(X); α̂_init)]
#' @param alpha_init Initial outcome model parameters
#' @param lambda L1 regularization parameter
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param calibrated Whether to use calibrated loss with truncation \mathcal{T}(·)
#'   (TRUE for calibrated loss eq:gamma_calibrated_loss; FALSE for its untruncated refined counterpart).
#'   Replaces the previous `use_truncation` parameter name for clarity.
#' @param M_tau Truncation threshold (only used if calibrated = TRUE)
#' @param W_outcome Outcome model features for psi' computation (n x p_outcome)
#'        Must be provided explicitly and have the same number of rows as \code{Z_site}.
#' @return Vector of density ratio parameters
#' @export
fit_unified_density_ratio <- function(Z_site, A, mean_grad_psi, alpha_init,
                                      lambda = LAMBDA_DEFAULT, max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                      calibrated = FALSE, M_tau = M_TAU_DEFAULT,
                                      W_outcome = NULL, A_val = 1L,
                                      family_int = 1L, link_int = 1L,
                                      warm_start = NULL) {
  if (is.null(W_outcome)) {
    stop("fit_unified_density_ratio: W_outcome must be provided explicitly.")
  }
  W_outcome_use <- as.matrix(W_outcome)
  if (nrow(W_outcome_use) != nrow(Z_site)) {
    stop(sprintf("fit_unified_density_ratio: row mismatch between Z_site (%d) and W_outcome (%d).",
                 nrow(Z_site), nrow(W_outcome_use)))
  }
  
  # Handle lambda = NULL case
  if (is.null(lambda)) {
    lmax <- compute_lambda_max_refined_dr(Z_site, A, mean_grad_psi, alpha_init,
                                           A_val = A_val,
                                           family_int = family_int, link_int = link_int,
                                           W_outcome = W_outcome,
                                           calibrated = calibrated, M_tau = M_tau)
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) 1e-4 else 0.01
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio)
    if (calibrated) {
      # Forward W_outcome to CV so ψ' uses correct features
      cv_result <- select_lambda_cv_calibrated_density_ratio_cpp(
        Z_site, A, mean_grad_psi, alpha_init, lambda_grid, N_CV_FOLDS_LAMBDA, max_iter, tol, M_tau,
        W_outcome_use, A_val, family_int, link_int
      )
    } else {
      # Refined DR should use ψ'(W_outcome^T alpha_init) with no truncation.
      # Reuse the calibrated CV kernel with M_tau = Inf to keep feature usage
      # (Z_site for gamma, W_outcome for ψ') consistent with final fitting.
      cv_result <- select_lambda_cv_calibrated_density_ratio_cpp(
        Z_site, A, mean_grad_psi, alpha_init, lambda_grid, N_CV_FOLDS_LAMBDA, max_iter, tol, Inf,
        W_outcome_use, A_val, family_int, link_int
      )
    }
    lambda <- cv_result$best_lambda
  }
  
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  # C++ param name: mean_grad_psi (positional), calibrated (positional)
  cpp_result <- fit_unified_density_ratio_cpp(
    Z_site, A, mean_grad_psi, alpha_init, lambda, max_iter, tol, 
    calibrated, M_tau, W_outcome_use, A_val, family_int, link_int, ws
  )
  result <- cpp_result$gamma
  attr(result, "lambda_used") <- lambda
  return(result)
}

#' Fit UNIFIED outcome model function (α parameters)
#'
#' Implements the weighted GLM loss for outcome regression following main.tex.
#' Supports GLM families gaussian and binomial.
#'
#' **Refined Loss (untruncated counterpart of eq:alpha_calibrated_loss, calibrated = FALSE):**
#' $$
#' \ell(\boldsymbol{\alpha}) = \widetilde{E}_{s_j}\left[\frac{I(A=1)}{\exp(\boldsymbol{\phi}^T \hat{\boldsymbol{\gamma}})} 
#'   (-Y \cdot \eta + \Psi(\eta))\right] + \lambda \|\boldsymbol{\alpha}\|_1
#' $$
#'
#' **Calibrated Loss (eq:alpha_calibrated_loss in main.tex, calibrated = TRUE):**
#' $$
#' \ell(\boldsymbol{\alpha}) = \widetilde{E}_{s_j}\left[\frac{I(A=1)}{\exp(\mathcal{T}(\boldsymbol{\phi}^T \hat{\boldsymbol{\gamma}}))} 
#'   (-Y \cdot \eta + \Psi(\eta))\right] + \lambda \|\boldsymbol{\alpha}\|_1
#' $$
#'
#' @param W_outcome Covariate matrix for outcome model (n × p)
#' @param Y Outcome vector (n × 1)
#' @param A Treatment indicator vector (n × 1)
#' @param A_val Treatment value to focus on (default 1)
#' @param gamma_s Density ratio parameters from previous step
#' @param lambda L1 regularization parameter (NULL for CV)
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param family GLM family ("gaussian", "binomial").
#'        Ignored if family_int and link_int are provided.
#' @param link Link function ("identity", "logit").
#'        Ignored if link_int is provided.
#' @param family_int Integer code for GLM family (0=gaussian, 1=binomial).
#'        If provided with link_int, skips string-based resolution for efficiency.
#' @param link_int Integer code for link function (0=identity, 1=logit).
#'        If provided with family_int, skips string-based resolution.
#' @param calibrated Whether to use calibrated loss with truncation \mathcal{T}(·)
#'   (TRUE for calibrated loss eq:alpha_calibrated_loss; FALSE for its untruncated refined counterpart).
#' @param M_tau Truncation threshold (only used if calibrated = TRUE)
#' @param Z_site Feature matrix for density ratio (Z_site). Must be provided explicitly.
#' @return Vector of outcome model parameters
#' @export
fit_unified_outcome <- function(W_outcome, Y, A, A_val = 1, gamma_s, lambda = NULL,
                                max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                family = "binomial", link = NULL,
                                family_int = NULL, link_int = NULL,
                                calibrated = FALSE, M_tau = M_TAU_DEFAULT, Z_site = NULL,
                                warm_start = NULL) {
  
  # Resolve GLM family/link: use integer codes directly if provided, else map from strings
  if (!is.null(family_int) && !is.null(link_int)) {
    family_int <- as.integer(family_int)
    link_int <- as.integer(link_int)
  } else {
    glm_spec <- resolve_glm_family(family)
    if (!is.null(link) && tolower(link) != glm_spec$link) {
      stop(sprintf("Explicit link override is not supported. family='%s' requires link='%s'.",
                   glm_spec$family, glm_spec$link))
    }
    family_int <- glm_spec$family_int
    link_int <- glm_spec$link_int
  }
  
  # Prepare Z_site BEFORE CV (so CV can use it for weight calculation)
  if (is.null(Z_site)) {
    stop("fit_unified_outcome: Z_site must be provided explicitly.")
  }
  Z_site_mat <- as.matrix(Z_site)
  if (nrow(Z_site_mat) != nrow(W_outcome)) {
    stop(sprintf("fit_unified_outcome: row mismatch between W_outcome (%d) and Z_site (%d).",
                 nrow(W_outcome), nrow(Z_site_mat)))
  }
  
  # Handle lambda = NULL case with cross-validation
  if (is.null(lambda)) {
    lmax <- compute_lambda_max_outcome(W_outcome, Y, A, gamma_s,
                                        A_val = A_val,
                                        family_int = family_int, link_int = link_int,
                                        Z_site = Z_site,
                                        calibrated = calibrated, M_tau = M_tau)
    lambda_min_ratio <- if (nrow(W_outcome) > ncol(W_outcome)) 1e-4 else 0.01
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio)
    
    # Use correct CV function based on calibrated flag
    # - Calibrated (calibrated=TRUE) uses select_lambda_cv_calibrated_outcome_cpp
    # - Refined (calibrated=FALSE) uses select_lambda_cv_general_refined_outcome_cpp
    if (calibrated) {
      # Calibrated outcome CV - uses truncation in loss function
      cv_result <- select_lambda_cv_calibrated_outcome_cpp(
        W_outcome, Y, A, gamma_s, lambda_grid, N_CV_FOLDS_LAMBDA, max_iter, tol, A_val, M_tau, Z_site_mat,
        family_int, link_int
      )
    } else {
      # Refined outcome CV - standard weighted loss without truncation
    cv_result <- select_lambda_cv_general_refined_outcome_cpp(
      W_outcome, Y, A, gamma_s, family_int, link_int, lambda_grid, N_CV_FOLDS_LAMBDA, max_iter, tol, A_val, Z_site_mat
    )
    }
    lambda <- cv_result$best_lambda
  }
  
  # Use unified C++ function with all parameters
  # C++ param name: use_truncation (positional)
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  cpp_result <- fit_unified_outcome_cpp(W_outcome, Y, A, gamma_s, family_int, link_int,
                                        lambda, max_iter, tol, A_val, 
                                        calibrated, M_tau, Z_site_mat, ws)
  # C++ returns 'alpha' matching paper notation (main.tex convention)
  result <- cpp_result$alpha
  attr(result, "lambda_used") <- lambda
  return(result)
}

# ============================================================================
# Shared helpers for weight optimization and variance aggregation
# ============================================================================

#' Convert C_cross to a matrix suitable for C++ functions
#'
#' Validates and converts a cross-site covariance specification to a K x K
#' matrix, or returns an empty matrix if invalid or NULL.
#'
#' @param C_cross Optional K x K cross-site covariance matrix, or NULL.
#' @param K Integer, expected dimension.
#' @return A K x K matrix (valid input) or a 0 x 0 matrix (NULL/invalid).
#' @keywords internal
prepare_cross_matrix <- function(C_cross, K) {
  if (!is.null(C_cross) && is.matrix(C_cross) &&
      nrow(C_cross) == K && ncol(C_cross) == K) {
    # Defensive normalization:
    # - Cross-site covariance is an off-diagonal object by design; diagonal is
    #   accounted for by V_t.
    # - Symmetrize to guard against small numerical asymmetries.
    C_cross <- 0.5 * (C_cross + t(C_cross))
    diag(C_cross) <- 0
    return(C_cross)
  }
  matrix(0, nrow = 0, ncol = 0)
}

#' Optimize aggregation weights with cross-site covariances
#'
#' Solves the convex optimization problem for optimal aggregation weights
#' following eq:final_opt in main.tex.
#'
#' **Optimization Problem (eq:final_opt in main.tex):**
#' $$
#' \min_{\boldsymbol{\eta}} \left\{
#'   \frac{1}{N_t}\left(1 - \sum_{j=1}^K \eta_j\right)^2 \hat{V}_{ot}
#'   + \sum_{j=1}^K \frac{\eta_j^2}{N_t} \hat{V}_{t,s_j}
#'   + \sum_{j=1}^K \frac{\eta_j^2}{N_{s_j}} \hat{V}_{s_j}
#'   + \frac{2}{N_t}\sum_{j=1}^K \eta_j\left(1 - \sum_{k=1}^K \eta_k\right)\hat{C}_{ot,s_j}
#'   + \frac{2}{N_t}\sum_{j<k} \eta_j \eta_k \hat{C}_{t,s_j,s_k}
#'   + \lambda \sum_{j=1}^K |\eta_j|(\hat{\mu}_{ot}^1 - \hat{\mu}_{t,s_j}^1)^2
#' \right\}
#' $$
#'
#' **Key Components:**
#' - First term: Target-only variance contribution
#' - Terms 2-3: Source-specific variance contributions (target and source sides)
#' - Term 4: Covariance between target-only and source-assisted estimators
#' - Term 5: **Cross-site covariances** (eq:cross_site_covariance in main.tex)
#' - Term 6: L1-like penalty for source selection based on estimate discrepancy
#'
#' **Critical Note (eq:variance_components in main.tex):**
#' Ignoring cross-site covariances $\hat{C}_{t,s_j,s_k}$ leads to systematic
#' underestimation of variance and incorrect confidence intervals!
#'
#' @param estimates Vector of source estimates $\hat{\mu}_{t,s_j}^1$ (length K)
#' @param variances List with $\hat{V}_{ot}$, $\hat{V}_t$ (vector), $\hat{V}_s$ (vector)
#' @param C_ot Vector of covariances $\hat{C}_{ot,s_j}$ (length K)
#' @param n_samples List with $N_t$ and $N_{s_j}$ (vector)
#' @param lambda Regularization parameter $\lambda$ for source selection
#' @param mu_ot Target-only estimate $\hat{\mu}_{ot}^1$
#' @param C_cross Cross-site covariance matrix $\hat{C}_{t,s_j,s_k}$ (K × K, optional)
#' @param clip_weights Whether to clip weights to [0,1]. Default FALSE to match
#'        the unconstrained optimization problem in the paper. Set TRUE only as a
#'        pragmatic finite-sample safeguard (it changes the theoretical optimum).
#' @return Optimal weights $\boldsymbol{\eta}^* = (\eta_1^*, \ldots, \eta_K^*)$
#' @export
optimize_weights <- function(estimates, variances, C_ot, n_samples,
                            lambda, mu_ot, C_cross = NULL, clip_weights = FALSE) {
  K <- length(estimates)

  # Defensive: consistent vector lengths
  if (length(variances$V_t) != K || length(variances$V_s) != K) {
    stop(sprintf("optimize_weights: estimates has length %d but V_t has %d, V_s has %d. They must match.",
                 K, length(variances$V_t), length(variances$V_s)))
  }
  if (length(C_ot) != K) {
    stop(sprintf("optimize_weights: C_ot has length %d but K=%d. They must match.", length(C_ot), K))
  }
  if (length(n_samples$n_s) != K) {
    stop(sprintf("optimize_weights: n_samples$n_s has length %d but K=%d. They must match.", length(n_samples$n_s), K))
  }

  # Input validation and clipping with more conservative variance bounds
  estimates <- pmax(pmin(estimates, ESTIMATE_MAX), -ESTIMATE_MAX)
  variances$V_t <- pmax(variances$V_t, VARIANCE_MIN)
  variances$V_s <- pmax(variances$V_s, VARIANCE_MIN)
  variances$V_ot <- max(variances$V_ot, VARIANCE_MIN)
  C_ot <- pmax(pmin(C_ot, ESTIMATE_MAX), -ESTIMATE_MAX)
  lambda <- max(min(lambda, LAMBDA_MAX), LAMBDA_MIN)
  mu_ot <- max(min(mu_ot, ESTIMATE_MAX), -ESTIMATE_MAX)
  
  # Prepare cross-site covariance matrix
  cross_matrix <- prepare_cross_matrix(C_cross, length(estimates))
  
  cpp_result <- optimize_weights_cpp(estimates, variances$V_t, variances$V_s, n_samples$n_s, 
                                    variances$V_ot, n_samples$n_t, C_ot, 
                                    lambda, mu_ot, WEIGHT_OPT_MAX_ITER, WEIGHT_OPT_TOL,
                                    cross_matrix, numeric(0))

  if (!isTRUE(cpp_result$converged)) {
    warning("optimize_weights: optimizer did not converge. Falling back to target-only aggregation (all source weights = 0).")
    return(rep(0, K))
  }
  
  # Validate output
  weights <- cpp_result$weights
  bad_idx <- which(is.na(weights) | is.infinite(weights))
  if (length(bad_idx) > 0) {
    warning(sprintf("optimize_weights: %d/%d weights were NA/Inf (indices: %s). Replaced with 0.",
                    length(bad_idx), length(weights), paste(bad_idx, collapse = ",")))
    weights[bad_idx] <- 0
  }
  
  # NOTE: The optimization problem (eq:final_opt in main.tex) is unconstrained.
  # Clipping to [0,1] is a numerical safeguard but changes the theoretical optimum.
  if (isTRUE(clip_weights)) {
    weights <- pmax(pmin(weights, 1.0), 0.0)  # Clip weights to [0, 1]
  }
  
  return(weights)
}

#' Aggregated variance calculation with cross-site covariances
#'
#' Computes the total variance of the aggregated estimator, accounting for
#' within-site variances and cross-site covariances.
#'
#' @param eta Numeric vector of aggregation weights (length K).
#' @param variances Numeric vector of within-site variance estimates (length K).
#' @param C_ot Numeric vector of target-site covariance terms (length K).
#' @param n_samples Integer vector of per-site sample sizes (length K).
#' @param C_cross Optional K x K matrix of cross-site covariances (default NULL).
#' @param lambda Not used in variance estimation; kept for API parity with the
#'   aggregation objective (default 0).
#' @param mu_ot Not used in variance estimation; kept for API parity with the
#'   aggregation objective (default 0).
#' @param mu_ts Numeric vector of source-site estimates (length K; if non-empty,
#'   ignored by design).
#' @return Scalar aggregated variance.
#' @export
calculate_aggregated_variance <- function(eta, variances, C_ot, n_samples, C_cross = NULL,
                                        lambda = 0.0, mu_ot = 0.0, mu_ts = numeric(0)) {
  
  # Validate and convert cross-site covariance matrix
  cross_matrix <- prepare_cross_matrix(C_cross, length(eta))
  
  # Penalty term mu_ts is excluded by design — it belongs in the optimization
  # objective (eq:final_opt in main.tex) but NOT in variance estimation.
  return(calculate_aggregated_variance_cpp(eta, variances$V_t, variances$V_s,
                                          n_samples$n_s, variances$V_ot,
                                          n_samples$n_t, C_ot, cross_matrix,
                                          lambda, mu_ot))
}
