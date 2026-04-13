# numerical_utils.R - Pure numerical and mathematical utility functions
#
# This file contains helper functions for numerical operations used throughout
# the FACE-C algorithm. All functions are pure (no side effects) and focused
# on mathematical transformations.
#
# Contents:
#   1. Vector/Matrix Transformations
#   2. Activation/Link Functions
#   3. Regularization Helpers
#   4. Safe Numerical Operations
#   5. Statistical Utilities
#   6. Logging Utilities
#   7. RNG State Management

# =============================================================================
# 1. VECTOR/MATRIX TRANSFORMATIONS
# =============================================================================

#' Null-coalescing operator (backport for R < 4.1.0)
#'
#' Returns x if it is not NULL, otherwise returns y.
#' @param x Value to test.
#' @param y Default value when x is NULL.
#' @return x if not NULL, else y.
#' @name op-null-default
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Scale and center a vector or matrix to zero mean and unit variance
#'
#' @param x Numeric vector or matrix.
#' @return Scaled and centered version of x.
#' @note When standard deviation is zero (all values identical), only centers
#'   the data without scaling to avoid division by zero.
#' @examples
#' scale_center(1:5)
#' scale_center(matrix(1:12, nrow = 3))
scale_center <- function(x) {
  if (is.matrix(x)) {
    return(apply(x, 2, function(col) {
      s <- sd(col)
      if (is.na(s) || s == 0) return(col - mean(col))
      return((col - mean(col)) / s)
    }))
  } else {
    s <- sd(x)
    if (is.na(s) || s == 0) return(x - mean(x))
    return((x - mean(x)) / s)
  }
}

#' Clip (truncate) values to a specified range
#'
#' @param x Numeric vector.
#' @param lower Lower bound (default -2.5).
#' @param upper Upper bound (default 2.5).
#' @return Vector with values clipped to [lower, upper].
#' @examples
#' clip_to_range(c(-5, 0, 5), -2, 2)  # Returns c(-2, 0, 2)
clip_to_range <- function(x, lower = -2.5, upper = 2.5) {
  pmax(pmin(x, upper), lower)
}

#' Normalize a vector to unit L2 norm
#'
#' @param x Numeric vector.
#' @return Normalized vector with ||x||_2 = 1, or original if norm is zero.
#' @examples
#' normalize_to_unit(c(3, 4))  # Returns c(0.6, 0.8)
normalize_to_unit <- function(x) {
  norm_x <- sqrt(sum(x^2))
  if (norm_x > 0) {
    x / norm_x
  } else {
    x
  }
}

#' Clip propensity scores to [PROP_SCORE_LOWER, PROP_SCORE_UPPER]
#'
#' Prevents extreme inverse probability weights by bounding propensity scores
#' away from 0 and 1. Uses constants from R/constants.R.
#'
#' @param ps Numeric vector of propensity scores.
#' @return Vector clipped to [PROP_SCORE_LOWER, PROP_SCORE_UPPER].
#' @examples
#' clip_propensity(c(0, 0.5, 1))  # Returns c(0.01, 0.5, 0.99)
clip_propensity <- function(ps) {
  pmax(pmin(ps, PROP_SCORE_UPPER), PROP_SCORE_LOWER)
}

#' Clip outcome model predictions based on GLM family
#'
#' Family-aware clipping to prevent numerical issues:
#'   - binomial: clip to [OUTCOME_PRED_LOWER, OUTCOME_PRED_UPPER] (prevents log(0))
#'   - gaussian: no clipping (continuous predictions can take any real value)
#'
#' @param pred Numeric vector of predictions (probabilities for binomial, means otherwise).
#' @param family GLM family string. Default "binomial".
#' @return Vector with family-appropriate clipping applied.
#' @examples
#' clip_outcome_pred(c(0, 0.5, 1))              # binomial: c(0.001, 0.5, 0.999)
#' clip_outcome_pred(c(-2, 0, 5), "gaussian")   # unchanged: c(-2, 0, 5)
clip_outcome_pred <- function(pred, family = "binomial") {
  switch(tolower(family),
    "binomial" = pmax(pmin(pred, OUTCOME_PRED_UPPER), OUTCOME_PRED_LOWER),
    "gaussian" = pred,
    stop(sprintf("Unsupported family for clip_outcome_pred(): '%s'", family))
  )
}

# =============================================================================
# 2. ACTIVATION/LINK FUNCTIONS
# =============================================================================

#' Logistic (sigmoid) function with numerical stability
#'
#' R reference implementation. For performance-critical code paths, prefer the
#' C++ version \code{logistic_cpp()} / \code{logistic_vec_cpp()} in utils.hpp.
#'
#' @param x Numeric value or vector.
#' @return Logistic transformation: 1 / (1 + exp(-x)).
#' @note Values are clipped to [-LOGISTIC_CLIP, LOGISTIC_CLIP] before computation
#'   to prevent overflow. LOGISTIC_CLIP is defined in R/constants.R and must match
#'   C++ NumericalConstants::ETA_CLIP_MAX in numerical_constants.hpp.
#' @examples
#' logistic(0)    # Returns 0.5
#' logistic(c(-10, 0, 10))
logistic <- function(x) {
  # Clip to prevent numerical overflow
  x <- pmax(pmin(x, LOGISTIC_CLIP), -LOGISTIC_CLIP)
  return(1 / (1 + exp(-x)))
}

# =============================================================================
# 3. SAFE NUMERICAL OPERATIONS
# =============================================================================

#' Safe logarithm with floor to prevent -Inf
#'
#' @param x Numeric value or vector.
#' @param epsilon Small positive floor value (default EPSILON_DEFAULT).
#' @return log(max(x, epsilon)).
#' @examples
#' safe_log(0)       # Returns log(EPSILON_DEFAULT) instead of -Inf
#' safe_log(c(0, 1, exp(1)))
safe_log <- function(x, epsilon = EPSILON_DEFAULT) {
  return(log(pmax(x, epsilon)))
}

#' Safe variance calculation with minimum floor
#'
#' Returns a minimum variance value for single observations or
#' when computed variance is too small/invalid.
#'
#' @param x Numeric vector.
#' @param min_var Minimum allowed variance (default EPSILON_DEFAULT).
#' @return Variance of x, floored at min_var.
#' @examples
#' safe_var(5)          # Returns EPSILON_DEFAULT (single value)
#' safe_var(c(1, 2, 3)) # Returns var(c(1, 2, 3))
safe_var <- function(x, min_var = EPSILON_DEFAULT) {
  if (length(x) <= 1) {
    return(min_var)
  }
  
  var_x <- var(x, na.rm = TRUE)
  if (is.na(var_x) || var_x <= 0) {
    return(min_var)
  }
  
  return(max(var_x, min_var))
}

#' Numerically stable sum via compensated summation (Neumaier)
#'
#' Useful for long reductions where naive summation may accumulate rounding error.
#' Returns the same mathematical quantity as \\code{sum(x)} but with improved
#' floating-point stability.
#'
#' @param x Numeric vector.
#' @return Stable sum of all elements in \\code{x}.
#' @examples
#' stable_sum_kahan(c(1e16, 1, -1e16))
stable_sum_kahan <- function(x) {
  if (!length(x)) return(0)
  s <- 0
  c <- 0
  for (value in x) {
    t <- s + value
    if (abs(s) >= abs(value)) {
      c <- c + ((s - t) + value)
    } else {
      c <- c + ((value - t) + s)
    }
    s <- t
  }
  s + c
}

#' Compute adaptive number of CV folds based on sample size
#'
#' Ensures each fold has at least \code{min_per_fold} observations
#' while respecting minimum and maximum fold count constraints.
#'
#' @param n Sample size.
#' @param min_per_fold Minimum observations per fold (default 5).
#' @param max_folds Maximum number of folds (default 5).
#' @param min_folds Minimum number of folds (default 4 for cv.glmnet validity).
#' @return Integer number of CV folds.
#' @examples
#' get_cv_fold_count(100)         # Returns 5
#' get_cv_fold_count(8)           # Returns 4
#' get_cv_fold_count(100, min_per_fold = 10) # Returns 5
get_cv_fold_count <- function(n, min_per_fold = 5, max_folds = N_CV_FOLDS_LAMBDA, min_folds = 4) {
  n <- as.integer(n)
  if (is.na(n) || n <= 0) return(as.integer(min_folds))
  candidate <- max(min_folds, min(max_folds, n %/% min_per_fold))
  as.integer(min(candidate, n))
}

#' Build a lambda grid for cross-validation (glmnet-style adaptive)
#'
#' Constructs a logarithmically spaced grid of lambda values for CV selection,
#' following the same strategy as \pkg{glmnet}.  The upper endpoint is
#' \code{lambda_max} (the smallest \eqn{\lambda} that zeros out all penalised
#' coefficients, computed from the KKT gradient at zero) and the lower
#' endpoint is \code{lambda_max * lambda_min_ratio}.
#'
#' @param lambda_max Data-adaptive upper bound.  Typically computed from the
#'   un-penalised gradient at \eqn{\beta = 0} via one of the
#'   \code{compute_lambda_max_*()} helpers.
#' @param lambda_min_ratio Ratio \eqn{\lambda_{\min}/\lambda_{\max}}.
#'   Default \code{1e-4} (matching \pkg{glmnet} for \eqn{n > p}).
#'   For \eqn{n \le p}, callers typically pass \code{0.01}.
#' @param nlambda Number of grid points (default
#'   \code{LAMBDA_GRID_SIZE_STANDARD}, i.e.\sspace 100).
#' @return Numeric vector of lambda values on a log-linear scale, sorted
#'   ascending (smallest first).
build_lambda_grid <- function(lambda_max, lambda_min_ratio = 1e-4,
                              nlambda = LAMBDA_GRID_SIZE_STANDARD) {
  stopifnot(is.numeric(lambda_max), length(lambda_max) == 1,
            is.finite(lambda_max), lambda_max > 0)
  lmin <- lambda_max * max(lambda_min_ratio, .Machine$double.eps)
  exp(seq(log(lmin), log(lambda_max), length.out = nlambda))
}

# =============================================================================
# Lambda-max computation helpers (glmnet-style KKT gradient at zero)
# =============================================================================

#' Compute lambda_max for INITIAL density ratio model
#'
#' At \eqn{\gamma = 0} the gradient of the un-penalised initial density
#' ratio loss (treating all source units as \eqn{A=1}) is
#' \deqn{\frac{\partial \ell}{\partial \gamma_j}\bigg|_{\gamma=0}
#'       = \bar\phi_j - \frac{1}{n_{\text{total}}} \sum_{i \in S} X_{ij},}
#' where \eqn{\bar\phi = \widetilde E_t[\phi(X)]} and the sum runs over the
#' source treated units.
#'
#' \eqn{\lambda_{\max}} is the smallest regularisation for which the entire
#' solution is zero (the KKT condition for L1 penalty, excluding the intercept).
#'
#' @param Z_site Covariate matrix (\eqn{n \times p}).
#' @param A Treatment indicator (length \eqn{n}).
#' @param mean_phi Target covariate mean \eqn{\widetilde E_t[\phi(X)]}
#'   (length \eqn{p+1}, first element = intercept = 1).
#' @param A_val Treatment value to select (default 1).
#' @return Scalar \eqn{\lambda_{\max}}.
compute_lambda_max_initial_dr <- function(Z_site, A, mean_phi, A_val = 1L) {
  n_total <- nrow(Z_site)
  treated_idx <- which(A == A_val)
  n_treated <- length(treated_idx)
  if (n_treated == 0) return(1)  # fallback

  # mean_phi already includes intercept at position 1
  # X_treated with intercept prepended, shape (n_treated, p+1)
  X_treated <- cbind(1, Z_site[treated_idx, , drop = FALSE])

  # Gradient at gamma=0 (exp(-0) = 1, psi' = 1 for initial DR):
  # grad_j = mean_phi[j] - sum(X_treated[, j]) / n_total
  grad <- as.numeric(mean_phi) - colSums(X_treated) / n_total

  # Exclude intercept (j=1): L1 penalty does not apply to intercept
  lambda_max <- max(abs(grad[-1]))
  max(lambda_max, 1e-6)  # floor to prevent degenerate grids
}

#' Compute lambda_max for REFINED / CALIBRATED density ratio model
#'
#' Same logic but the gradient depends on \eqn{\psi'(\phi^T \alpha_{\text{init}})}.
#' At \eqn{\gamma = 0}:
#' \deqn{
#'   \nabla_j \ell |_{\gamma=0} = \bar g_j
#'     - \frac{1}{n_{\text{total}}} \sum_{i \in S} X_{ij}\,\psi'(\phi_i^T\alpha)
#' }
#' where \eqn{\bar g = \widetilde E_t[\nabla_\alpha \psi]}.
#'
#' @param Z_site Covariate matrix (\eqn{n \times p}).
#' @param A Treatment indicator (length \eqn{n}).
#' @param mean_grad_psi Mean gradient \eqn{\widetilde E_t[\nabla_\alpha\psi]}
#'   (length \eqn{p+1}).
#' @param alpha_init Initial outcome parameters (length \eqn{p+1}).
#' @param A_val Treatment value to select (default 1).
#' @param family_int GLM family integer.
#' @param link_int GLM link integer.
#' @param W_outcome Outcome features (optional, defaults to Z_site).
#' @param calibrated Whether truncation is used.
#' @param M_tau Truncation parameter.
#' @return Scalar \eqn{\lambda_{\max}}.
compute_lambda_max_refined_dr <- function(Z_site, A, mean_grad_psi, alpha_init,
                                           A_val = 1L,
                                           family_int = 1L, link_int = 1L,
                                           W_outcome = NULL,
                                           calibrated = FALSE,
                                           M_tau = M_TAU_DEFAULT) {
  n_total <- nrow(Z_site)
  treated_idx <- which(A == A_val)
  n_treated <- length(treated_idx)
  if (n_treated == 0) return(1)

  # X_treated with intercept (density ratio features)
  X_treated <- cbind(1, Z_site[treated_idx, , drop = FALSE])

  # Determine outcome features for ψ'
  if (!is.null(W_outcome) && nrow(W_outcome) == n_total) {
    W_tr <- cbind(1, W_outcome[treated_idx, , drop = FALSE])
  } else {
    W_tr <- X_treated
  }

  # Compute ψ'(ϕ^T α_init) for each treated unit
  eta_alpha <- as.numeric(W_tr %*% alpha_init)
  if (calibrated) {
    eta_alpha <- pmin(pmax(eta_alpha, -M_tau), M_tau)
  }
  # ψ'(η): in density-ratio loss this corresponds to b''(θ) and must be non-negative.
  # Keep R lambda_max logic aligned with C++ implementation, which uses |h'(η)|.
  psi_prime <- switch(as.character(link_int),
    "0" = rep(1.0, length(eta_alpha)),             # identity: h'(η) = 1
    "1" = { h <- logistic(eta_alpha); h * (1 - h) }, # logit: h'(η) = h(1-h)
    "2" = exp(eta_alpha),                            # log: h'(η) = exp(η)
    "3" = {
      eta_safe <- ifelse(abs(eta_alpha) < sqrt(.Machine$double.eps),
                         sign(eta_alpha + (eta_alpha == 0)) * sqrt(.Machine$double.eps),
                         eta_alpha)
      1.0 / (eta_safe^2)                              # inverse: |h'(η)| = 1/η²
    },
    rep(1.0, length(eta_alpha))                       # fallback: identity
  )

  # grad_j = mean_grad_psi[j] - sum(X_treated[, j] * psi_prime) / n_total
  grad <- as.numeric(mean_grad_psi) - colSums(X_treated * psi_prime) / n_total

  lambda_max <- max(abs(grad[-1]))
  max(lambda_max, 1e-6)
}

#' Compute lambda_max for WEIGHTED outcome GLM model
#'
#' At \eqn{\beta = 0} the (un-penalised) gradient of the weighted GLM loss is
#' \deqn{
#'   \nabla_j \ell|_{\beta=0} = \frac{1}{n_{\text{total}}}
#'     \sum_{i \in S} w_i \bigl(\psi'(0) - Y_i\bigr) X_{ij}
#' }
#' where \eqn{w_i = \exp(-Z_i^T\gamma)} and \eqn{\psi'(0)} depends on the
#' GLM family (e.g.\ 0.5 for logistic).
#'
#' @param W_outcome Covariate matrix for outcome model (\eqn{n \times p}).
#' @param Y Outcome vector.
#' @param A Treatment indicator.
#' @param gamma_s Density ratio parameters (length \eqn{p_{\text{site}} + 1}).
#' @param A_val Treatment value to select (default 1).
#' @param family_int GLM family integer (0=gaussian, 1=binomial).
#' @param link_int Link integer (0=identity, 1=logit).
#' @param Z_site Site features for gamma (if NULL, uses W_outcome).
#' @param calibrated Whether truncation is used.
#' @param M_tau Truncation parameter.
#' @return Scalar \eqn{\lambda_{\max}}.
compute_lambda_max_outcome <- function(W_outcome, Y, A, gamma_s,
                                        A_val = 1L,
                                        family_int = 1L, link_int = 1L,
                                        Z_site = NULL,
                                        calibrated = FALSE, M_tau = M_TAU_DEFAULT) {
  n_total <- nrow(W_outcome)
  treated_idx <- which(A == A_val)
  n_treated <- length(treated_idx)
  if (n_treated == 0) return(1)

  X_tr <- cbind(1, W_outcome[treated_idx, , drop = FALSE])
  Y_tr <- Y[treated_idx]

  # Compute weights from density ratio
  if (!is.null(Z_site) && nrow(Z_site) == n_total) {
    Z_tr <- cbind(1, Z_site[treated_idx, , drop = FALSE])
  } else {
    Z_tr <- cbind(1, W_outcome[treated_idx, , drop = FALSE])
  }

  g_val <- as.numeric(Z_tr %*% gamma_s[1:ncol(Z_tr)])
  if (calibrated) {
    g_val <- pmin(pmax(g_val, -M_tau), M_tau)
  }
  weights <- exp(-g_val)
  weights <- pmin(pmax(weights, 1e-10), 1e10)

  # h(0) = ψ'(0) for the GLM at β = 0: the response (mean) function, not its derivative.
  # The GLM loss gradient is ∇_j ℓ|_{β=0} = (1/n) Σ w_i (h(0) - Y_i) X_{ij}
  # where h(η) = H'(η) is the inverse link (response) function.
  psi_prime_0 <- if (family_int == 1L) {
    0.5           # binomial/logit: h(0) = logistic(0) = 0.5
  } else {
    0.0           # gaussian/identity: h(0) = 0
  }

  # grad_j = (1/n_total) * sum(w_i * (ψ'(0) - Y_i) * X_ij)
  residuals <- weights * (psi_prime_0 - Y_tr)
  grad <- colSums(X_tr * residuals) / n_total

  lambda_max <- max(abs(grad[-1]))
  max(lambda_max, 1e-6)
}

# =============================================================================
# 4. STATISTICAL UTILITIES
# =============================================================================

#' Compute empirical expectation (weighted mean)
#'
#' @param x Data values.
#' @param weights Optional weights (default NULL for unweighted).
#' @return Weighted mean of x.
#' @examples
#' empirical_expectation(1:5)           # Returns 3
#' empirical_expectation(1:3, c(1,2,1)) # Returns 2
empirical_expectation <- function(x, weights = NULL) {
  if (is.null(weights)) {
    mean(x)
  } else {
    sum(x * weights) / sum(weights)
  }
}

# =============================================================================
# 5. LOGGING UTILITIES
# =============================================================================

#' Conditional logging helper
#'
#' Prints a formatted message only when verbose is TRUE.
#' Replaces scattered `if (verbose) cat(sprintf(...))` blocks throughout
#' the codebase, keeping algorithm logic free of presentation concerns.
#'
#' @param verbose Logical; if FALSE, message is suppressed.
#' @param fmt Format string (passed to \code{sprintf}).
#' @param ... Values interpolated into \code{fmt}.
#' @return Invisible NULL.
#' @examples
#' log_info(TRUE, "Processing fold %d/%d\n", 1, 5)
#' log_info(FALSE, "This is suppressed\n")
log_info <- function(verbose, fmt, ...) {
  if (isTRUE(verbose)) cat(sprintf(fmt, ...))
  invisible(NULL)
}

# =============================================================================
# 6. RNG STATE MANAGEMENT
# =============================================================================

#' Execute an expression with a temporary RNG seed
#'
#' Saves the current RNG state, sets the given seed, evaluates \code{expr},
#' then restores the original RNG state. This ensures deterministic generation
#' inside \code{expr} without polluting the caller's random stream.
#'
#' @param seed Integer seed to use.
#' @param expr Expression to evaluate under the temporary seed.
#' @return The result of evaluating \code{expr}.
#' @examples
#' set.seed(42); runif(1)  # 0.9148...
#' with_seed(99, rnorm(3)) # deterministic under seed 99
#' runif(1)                 # continues from seed 42 stream
with_seed <- function(seed, expr) {
  # Save current RNG state
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv)
  if (old_seed_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  eval.parent(substitute(expr))
}

# =============================================================================
# 7. BLOCK MATRIX UTILITIES
# =============================================================================
