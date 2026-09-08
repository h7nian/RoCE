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
# By default, we estimate μ¹_t = E_t[Y(1)] (treated potential outcome mean).
# Set A_val = 0 to estimate μ⁰_t = E_t[Y(0)] (control potential outcome mean).
# The primary run_tate_crossfit() wrapper runs both arms on common folds and
# contrasts their influence components before estimating TATE-level weights.
#
# ============================================================================
# ARM-SPECIFIC BUILDING BLOCK: μ¹ = E_t[Y(1)]
# ============================================================================
# These nuisance functions estimate one potential-outcome arm at a time. The
# TATE wrapper calls them for A_val=1 and A_val=0 and then uses a common
# aggregation layer; callers should not infer TATE optimality by subtracting
# two separately optimized arm-wise aggregators.
#
# ============================================================================
# NOTATION MAPPING (Code Variable → Paper Symbol in main.tex)
# ============================================================================
# - alpha / alpha_init → α (outcome model parameters)
# - gamma / gamma_s → γ (site/treatment model / density ratio parameters)
# - alpha_ot → α_{ot} (target-only outcome model parameters)
# - X / W_outcome → W(X) (basis expansion for outcome model)
# - Z_site → Z(X) (basis expansion for site/treatment / density-ratio model)
#   main.tex uses a single φ(X); the code allows W(X) and Z(X) to differ.
#   Calibrated density-ratio target summaries therefore use
#   E_t[h'(W(X)^T α) * Z̃(X)] so the moment has the same dimension as γ.
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
#    - optimize_weights(): Optimal aggregation weights (eq:agg_penalized_objective in main.tex)
#    - calculate_aggregated_variance(): Aggregated variance with cross-site cov
#
# ============================================================================

# Dependencies: constants.R, numerical_utils.R
# All loaded automatically by the R package system (see DESCRIPTION Collate field).

.validate_A_val <- function(A_val, caller) {
  if (length(A_val) != 1L || is.na(A_val) || !is.numeric(A_val) || !(A_val %in% c(0, 1))) {
    stop(sprintf("%s: A_val must be a scalar 0 or 1.", caller))
  }
  as.integer(A_val)
}

.validate_lambda_scalar <- function(lambda, caller, arg = "lambda", allow_zero = TRUE) {
  if (length(lambda) != 1L || is.na(lambda) || !is.numeric(lambda) || !is.finite(lambda)) {
    stop(sprintf("%s: %s must be a finite numeric scalar.", caller, arg))
  }
  lambda <- as.numeric(lambda)
  if (allow_zero) {
    if (lambda < 0) stop(sprintf("%s: %s must be non-negative.", caller, arg))
  } else {
    if (lambda <= 0) stop(sprintf("%s: %s must be positive.", caller, arg))
  }
  lambda
}

.match_nuisance_lambda_rule <- function(lambda_rule, caller, arg = "lambda_rule") {
  if (missing(lambda_rule) || is.null(lambda_rule)) {
    return("min")
  }
  tryCatch(
    match.arg(lambda_rule, c("min", "1se")),
    error = function(e) {
      stop(sprintf("%s: %s must be either 'min' or '1se'.", caller, arg),
           call. = FALSE)
    }
  )
}

.nuisance_cv_fold_count <- function(A, A_val, caller) {
  A_val <- .validate_A_val(A_val, caller)
  n_arm <- sum(A == A_val, na.rm = TRUE)
  if (n_arm < 2L) {
    stop(sprintf("%s: at least two observations with A_val=%d are required for CV lambda selection; found %d.",
                 caller, A_val, n_arm), call. = FALSE)
  }
  get_cv_fold_count(n_arm)
}

.select_nuisance_cv_lambda <- function(cv_result, lambda_rule, caller) {
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, caller)
  if (!is.list(cv_result)) {
    stop(sprintf("%s: CV result must be a list with glmnet-style lambda metadata.",
                 caller), call. = FALSE)
  }
  required <- c("lambda_min", "lambda_1se")
  missing_required <- required[!vapply(required, function(field) {
    !is.null(cv_result[[field]]) &&
      length(cv_result[[field]]) == 1L &&
      is.finite(cv_result[[field]])
  }, logical(1L))]
  if (length(missing_required) > 0L) {
    stop(sprintf(
      "%s: CV result is missing finite glmnet-style field(s): %s.",
      caller, paste(missing_required, collapse = ", ")
    ), call. = FALSE)
  }
  field <- if (identical(lambda_rule, "1se")) "lambda_1se" else "lambda_min"
  lambda <- as.numeric(cv_result[[field]])
  attr(lambda, "lambda_min") <- as.numeric(cv_result$lambda_min)
  attr(lambda, "lambda_1se") <- as.numeric(cv_result$lambda_1se)
  attr(lambda, "lambda_rule") <- lambda_rule
  for (diagnostic_name in c(
    "invalid_fold_fits",
    "invalid_lambdas",
    "path_tail_skipped_fold_fits"
  )) {
    diagnostic_value <- suppressWarnings(as.integer(cv_result[[diagnostic_name]]))
    if (length(diagnostic_value) != 1L || is.na(diagnostic_value) ||
        diagnostic_value < 0L) {
      stop(sprintf(
        "%s: CV result has invalid non-negative count field '%s'.",
        caller, diagnostic_name
      ), call. = FALSE)
    }
    attr(lambda, paste0("cv_", diagnostic_name)) <- diagnostic_value
  }
  lambda
}

.attach_nuisance_lambda_attrs <- function(coefs, lambda) {
  attr(coefs, "lambda_used") <- as.numeric(lambda)
  attr(coefs, "lambda_min") <- if (!is.null(attr(lambda, "lambda_min"))) attr(lambda, "lambda_min") else as.numeric(lambda)
  attr(coefs, "lambda_1se") <- if (!is.null(attr(lambda, "lambda_1se"))) attr(lambda, "lambda_1se") else as.numeric(lambda)
  attr(coefs, "lambda_rule") <- if (!is.null(attr(lambda, "lambda_rule"))) attr(lambda, "lambda_rule") else "fixed"
  for (attribute_name in c(
    "lambda_selected_on_path",
    "lambda_selected_before_support_floor",
    "cv_invalid_fold_fits",
    "cv_invalid_lambdas",
    "cv_path_tail_skipped_fold_fits",
    "support_penalty_floor",
    "support_penalty_floor_applied",
    "support_path_floor_applied",
    "support_constant_feature_count"
  )) {
    attr(coefs, attribute_name) <- attr(lambda, attribute_name)
  }
  coefs
}

# A source-arm feature that is constant creates an intercept/feature null
# direction in an exponential-tilting loss.  Along that direction the source
# exponential term is unchanged.  If the target-side linear moment gap is
# larger than the L1 slope penalty, the finite-sample objective is unbounded.
# Return the smallest penalty (plus a tiny strict margin) that makes every such
# coordinate direction coercive.  This is a numerical safeguard for sparse
# arm/site cells; it does not assert population overlap when empirical support
# is absent.
.density_ratio_support_penalty_floor <- function(
    Z_site, A, linear_moment, A_val, caller,
    constant_tolerance = sqrt(.Machine$double.eps),
    strict_margin = 1e-6) {
  Z_site <- as.matrix(Z_site)
  A <- as.numeric(A)
  linear_moment <- as.numeric(linear_moment)
  if (!is.numeric(Z_site) || nrow(Z_site) != length(A) ||
      length(linear_moment) != ncol(Z_site) + 1L ||
      any(!is.finite(Z_site)) || any(!is.finite(A)) ||
      any(!is.finite(linear_moment))) {
    stop(
      caller,
      ": support-floor inputs must be finite with matching rows and moment dimension.",
      call. = FALSE
    )
  }
  arm_rows <- A == A_val
  if (!any(arm_rows)) {
    stop(caller, ": no source observations are available for A_val.",
         call. = FALSE)
  }
  if (ncol(Z_site) == 0L) {
    return(list(floor = 0, constant_feature_count = 0L))
  }

  arm_design <- Z_site[arm_rows, , drop = FALSE]
  column_min <- apply(arm_design, 2L, min)
  column_max <- apply(arm_design, 2L, max)
  column_scale <- pmax(1, abs(column_min), abs(column_max))
  constant <- (column_max - column_min) <=
    constant_tolerance * column_scale
  if (!any(constant)) {
    return(list(floor = 0, constant_feature_count = 0L))
  }

  source_constant <- (column_min[constant] + column_max[constant]) / 2
  moment_gap <- linear_moment[-1L][constant] -
    source_constant * linear_moment[[1L]]
  raw_floor <- max(abs(moment_gap), 0)
  safe_floor <- if (raw_floor > 0) {
    raw_floor + strict_margin * max(1, raw_floor)
  } else {
    0
  }
  list(
    floor = safe_floor,
    constant_feature_count = as.integer(sum(constant))
  )
}

.apply_density_ratio_support_floor <- function(
    lambda, Z_site, A, linear_moment, A_val, caller,
    support = NULL, path_floor_applied = FALSE) {
  if (is.null(support)) {
    support <- .density_ratio_support_penalty_floor(
      Z_site, A, linear_moment, A_val, caller
    )
  }
  selected <- as.numeric(lambda)
  used <- max(selected, support$floor)
  result <- used
  attr(result, "lambda_min") <- attr(lambda, "lambda_min") %||% selected
  attr(result, "lambda_1se") <- attr(lambda, "lambda_1se") %||% selected
  attr(result, "lambda_rule") <- attr(lambda, "lambda_rule") %||% "fixed"
  attr(result, "cv_invalid_fold_fits") <-
    attr(lambda, "cv_invalid_fold_fits") %||% 0L
  attr(result, "cv_invalid_lambdas") <-
    attr(lambda, "cv_invalid_lambdas") %||% 0L
  attr(result, "cv_path_tail_skipped_fold_fits") <-
    attr(lambda, "cv_path_tail_skipped_fold_fits") %||% 0L
  attr(result, "lambda_selected_on_path") <- selected
  attr(result, "lambda_selected_before_support_floor") <- if (
    isTRUE(path_floor_applied)
  ) {
    NA_real_
  } else {
    selected
  }
  attr(result, "support_penalty_floor") <- support$floor
  attr(result, "support_penalty_floor_applied") <-
    isTRUE(path_floor_applied) ||
    support$floor > selected * (1 + 10 * .Machine$double.eps)
  attr(result, "support_path_floor_applied") <- isTRUE(path_floor_applied)
  attr(result, "support_constant_feature_count") <-
    support$constant_feature_count
  result
}

.constrain_density_ratio_lambda_grid <- function(lambda_grid, support) {
  constrained <- pmax(as.numeric(lambda_grid), support$floor)
  attr(constrained, "support_path_floor_applied") <-
    any(as.numeric(lambda_grid) < support$floor)
  constrained
}

.attach_nuisance_fit_attrs <- function(
    coefs, cpp_result, cv_seconds = 0, final_fit_seconds = NA_real_) {
  attr(coefs, "converged") <- isTRUE(cpp_result$converged)
  iterations <- suppressWarnings(as.integer(cpp_result$iterations))
  if (length(iterations) != 1L || is.na(iterations) || iterations < 0L) {
    iterations <- NA_integer_
  }
  attr(coefs, "iterations") <- iterations
  scalar_diagnostic <- function(field) {
    value <- suppressWarnings(as.numeric(cpp_result[[field]]))
    if (length(value) != 1L || !is.finite(value) || value < 0) {
      NA_real_
    } else {
      value
    }
  }
  max_update <- scalar_diagnostic("max_update")
  convergence_threshold <- scalar_diagnostic("convergence_threshold")
  attr(coefs, "max_update") <- max_update
  attr(coefs, "convergence_threshold") <- convergence_threshold
  attr(coefs, "update_to_threshold_ratio") <- if (
      is.finite(max_update) && is.finite(convergence_threshold) &&
      convergence_threshold > 0) {
    max_update / convergence_threshold
  } else {
    NA_real_
  }
  attr(coefs, "max_abs_coefficient") <-
    scalar_diagnostic("max_abs_coefficient")
  line_search_failures <- suppressWarnings(
    as.integer(cpp_result$line_search_failures)
  )
  if (length(line_search_failures) != 1L || is.na(line_search_failures) ||
      line_search_failures < 0L) {
    line_search_failures <- NA_integer_
  }
  attr(coefs, "line_search_failures") <- line_search_failures
  attr(coefs, "cv_seconds") <- as.numeric(cv_seconds)
  attr(coefs, "final_fit_seconds") <- as.numeric(final_fit_seconds)
  coefs
}

# ============================================================================
# R-BASED GLM FITTING (using glmnet)
# ============================================================================

#' Fit initial outcome model using L1-regularized GLM
#'
#' Implements the initial outcome model estimation for the two-level cross-fitting
#' algorithm. Following eq:nuisance_initial_losses in main.tex:
#'
#' \deqn{
#' \ell(\boldsymbol{\alpha}_{t,s_j}^{1}) = \widetilde{E}_{s_j}\left[I(A=1)
#'   \ell_{\text{GLM}}\left(Y, \psi(\phi(\mathbf{X}); \boldsymbol{\alpha}_{t,s_j}^1)\right)\right]
#'   + \lambda_\alpha \|\boldsymbol{\alpha}_{t,s_j}^{1}\|_1.
#' }
#'
#' This function uses `glmnet` for efficient L1-regularized GLM fitting.
#' Lambda is selected via cross-validation if not provided.
#'
#' @param W_outcome Covariate matrix (n x p), the design matrix phi(X)
#' @param Y Outcome vector (n x 1)
#' @param A Treatment indicator vector (n x 1), binary \{0, 1\}
#' @param A_val Treatment value to fit (default 1 for treated group)
#' @param lambda L1 regularization parameter on the full-source empirical
#'        objective scale. Internally this is converted to glmnet's
#'        treatment-arm-only scale after filtering to \code{A == A_val}.
#' @param nlambda Number of lambda values in glmnet path (default 100).
#'        Lower values (e.g., 20) speed up CV with minimal precision loss.
#' @param family GLM family: "gaussian" or "binomial".
#'        Default "binomial".
#' @param lambda_rule CV selection rule when \code{lambda = NULL}: \code{"min"}
#'   selects \code{lambda.min}; \code{"1se"} selects \code{lambda.1se}.
#'
#' @return Vector of outcome model parameters alpha (including intercept as first element)
#'
#' @details
#' The function:
#' 1. Filters data to the specified treatment arm (A = A_val)
#' 2. Uses cv.glmnet for lambda selection if lambda is NULL
#' 3. Converts between glmnet's arm-only loss scale and the paper's
#'    full-source empirical loss scale
#' 4. Fits a GLM with L1 penalty using the specified family
#' 5. Returns coefficients including the intercept
#'
#' @seealso \code{\link{fit_unified_outcome}} for C++ accelerated outcome fitting
#' @export
#' @param cv_group_id Optional positive integer observation-origin IDs aligned
#'   with all input rows. Repeated IDs stay in one nuisance CV fold; NULL or
#'   all-unique IDs preserve the ordinary CV path. Used for resampled data.
fit_initial_outcome <- function(W_outcome, Y, A, A_val = 1L, lambda = NULL, 
                                nlambda = LAMBDA_GRID_SIZE_STANDARD, family = "binomial",
                                lambda_rule = c("min", "1se"), cv_group_id = NULL) {
  
  A_val <- .validate_A_val(A_val, "fit_initial_outcome")
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "fit_initial_outcome")
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, nrow(W_outcome), "fit_initial_outcome"
  )
  # A cached lambda can arrive non-finite (NA) from a degenerate fold that fell
  # back to the constant nuisance; treat it as absent so this fold selects its
  # own lambda by CV (or falls back itself) instead of erroring on the sentinel.
  if (!is.null(lambda) && !(length(lambda) == 1L && is.finite(lambda))) {
    lambda <- NULL
  }
  if (!is.null(lambda)) {
    lambda <- .validate_lambda_scalar(lambda, "fit_initial_outcome")
  }

  # Resolve GLM family for glmnet
  glm_spec <- resolve_glm_family(family)
  glmnet_family <- glm_spec$glmnet_family
  
  # Filter to treatment arm
  arm_idx <- which(A == A_val)
  if (length(arm_idx) == 0) {
    stop(sprintf("fit_initial_outcome: no observations with A_val=%d (n=%d). Outcome model is not identifiable for this fold.",
                 A_val, nrow(W_outcome)))
  }
  
  X_arm <- W_outcome[arm_idx, , drop = FALSE]
  Y_arm <- Y[arm_idx]
  
  # Ensure minimum sample size for CV
  n_arm <- length(Y_arm)
  arm_fraction <- n_arm / nrow(W_outcome)

  # Degenerate (single-class) binomial response: a penalized logistic CV fit is
  # not identifiable (glmnet errors "one ... class has 1 or 0 observations"),
  # e.g. a saturated high-dimensional fold where every treated outcome is 1.
  # Fall back to the constant outcome mean as intercept-only coefficients -- the
  # correct degenerate nuisance -- and count it for the saturation diagnostic.
  degenerate_alpha <- function() {
    .record_or_degenerate_fold()
    m   <- mean(Y_arm)
    eps <- 1e-6
    a   <- c(stats::qlogis(min(max(m, eps), 1 - eps)), rep(0, ncol(X_arm)))
    attr(a, "lambda_used") <- NA_real_
    attr(a, "lambda_fit")  <- NA_real_
    attr(a, "lambda_min")  <- NA_real_
    attr(a, "lambda_1se")  <- NA_real_
    attr(a, "lambda_rule") <- "degenerate_constant"
    a
  }
  # A near-saturated outcome (tiny minority class) at high dimension makes the
  # treated response (near-)perfectly separable, so the penalized logistic CV is
  # unstable: glmnet either errors ("1 or 0 observations") or fails inside
  # cv.glmnet with a non-conformable error from inconsistent per-fold lambda
  # paths. Below glmnet's own <8-per-class threshold the fitted OR is
  # uninformative, so use the constant outcome mean (the correct degenerate
  # nuisance), counted for the saturation diagnostic.
  if (glmnet_family == "binomial" &&
      min(table(factor(Y_arm, levels = c(0, 1)))) < 8L) {
    return(degenerate_alpha())
  }

  # ---------- Standard path: glmnet L1-regularized GLM -----------------------
  if (is.null(lambda)) {
    # Use cross-validation to select lambda
    lambda_fit <- tryCatch({
      # Adaptive number of folds based on sample size
      n_cv_folds <- get_cv_fold_count(n_arm)
      
      cv_args <- list(
        x = X_arm,
        y = Y_arm,
        family = glmnet_family,
        alpha = 1,  # Lasso penalty
        standardize = TRUE,
        nfolds = n_cv_folds,
        nlambda = nlambda,  # Control lambda path length
        maxit = GLMNET_MAX_ITER
      )
      cv_fold_id <- .make_nuisance_cv_fold_id(
        if (is.null(cv_group_id)) NULL else cv_group_id[arm_idx],
        n_cv_folds, "fit_initial_outcome"
      )
      if (!is.null(cv_fold_id)) cv_args$foldid <- cv_fold_id
      cv_fit <- do.call(glmnet::cv.glmnet, cv_args)
      lambda_selected <- switch(lambda_rule,
        min = cv_fit$lambda.min,
        `1se` = cv_fit$lambda.1se
      )
      attr(lambda_selected, "lambda_min") <- cv_fit$lambda.min
      attr(lambda_selected, "lambda_1se") <- cv_fit$lambda.1se
      lambda_selected
    }, error = function(e) {
      msg <- conditionMessage(e)
      # A CV sub-fold can be single-class even when the full arm is not; signal
      # the constant degenerate fallback handled just below.
      # Safety net for any residual saturation-driven cv.glmnet failure that
      # slipped past the pre-check above (a single-class CV sub-fold, or a
      # non-conformable from inconsistent per-fold lambda paths on near-separable
      # data): fall back to the constant degenerate nuisance handled just below.
      if (glmnet_family == "binomial" &&
          grepl("1 or 0 observations|0 or 1 observations|non-conformable", msg)) {
        return(NULL)
      }
      stop(sprintf("fit_initial_outcome: cv.glmnet failed (family=%s, n_arm=%d, nlambda=%d): %s",
                   glmnet_family, n_arm, nlambda, msg))
    })
    if (is.null(lambda_fit)) return(degenerate_alpha())
    lambda_use <- as.numeric(lambda_fit) * arm_fraction
    lambda_min <- as.numeric(attr(lambda_fit, "lambda_min")) * arm_fraction
    lambda_1se <- as.numeric(attr(lambda_fit, "lambda_1se")) * arm_fraction
    lambda_fit <- as.numeric(lambda_fit)
    lambda_rule_used <- lambda_rule
  } else {
    lambda_use <- as.numeric(lambda)
    lambda_fit <- as.numeric(lambda_use / arm_fraction)
    lambda_min <- lambda_use
    lambda_1se <- lambda_use
    lambda_rule_used <- "fixed"
  }
  
  # Fit final model with selected lambda
  fit <- glmnet(
    x = X_arm,
    y = Y_arm,
    family = glmnet_family,
    alpha = 1,
    lambda = lambda_fit,
    standardize = TRUE,
    maxit = GLMNET_MAX_ITER
  )
  
  # Extract coefficients (includes intercept as first element)
  alpha <- as.vector(coef(fit, s = lambda_fit))
  
  # Return the full-source-scale lambda for caching across folds whose arm
  # prevalence may differ. lambda_fit records the glmnet-scale value used here.
  attr(alpha, "lambda_used") <- as.numeric(lambda_use)
  attr(alpha, "lambda_fit") <- as.numeric(lambda_fit)
  attr(alpha, "lambda_min") <- as.numeric(lambda_min)
  attr(alpha, "lambda_1se") <- as.numeric(lambda_1se)
  attr(alpha, "lambda_rule") <- lambda_rule_used
  
  return(alpha)
}

#' Mean GLM gradient expressed on the density-ratio basis
#'
#' Computes the target-site summary used in the calibrated density-ratio
#' moment when the outcome basis \code{W_outcome} and site basis \code{Z_site}
#' differ:
#' \deqn{E_t[h'(W(X)^T\alpha) \tilde Z(X)]}
#' where \eqn{\tilde Z(X)} includes the intercept. This is the dual-basis
#' counterpart of \code{mean_glm_gradient_cpp()}, which returns the same
#' quantity on the outcome basis \eqn{\tilde W(X)}.
#'
#' @param W_outcome Outcome-model features used with \code{alpha}.
#' @param Z_site Site/density-ratio features used with \code{gamma}.
#' @param alpha Outcome-model parameters, including intercept.
#' @param family_int Integer GLM family code, accepted for API symmetry.
#' @param link_int Integer link code.
#' @return Numeric vector of length \code{ncol(Z_site) + 1}.
#' @keywords internal
.mean_glm_gradient_site_basis <- function(W_outcome, Z_site, alpha,
                                          family_int = FAMILY_BINOMIAL,
                                          link_int = LINK_LOGIT) {
  W_mat <- as.matrix(W_outcome)
  Z_mat <- as.matrix(Z_site)
  alpha <- as.numeric(alpha)
  link_int <- as.integer(link_int)

  if (nrow(W_mat) != nrow(Z_mat)) {
    stop(sprintf(".mean_glm_gradient_site_basis: row mismatch between W_outcome (%d) and Z_site (%d).",
                 nrow(W_mat), nrow(Z_mat)))
  }
  if (length(alpha) != ncol(W_mat) + 1L) {
    stop(sprintf(".mean_glm_gradient_site_basis: alpha has length %d, expected %d for W_outcome.",
                 length(alpha), ncol(W_mat) + 1L))
  }
  if (nrow(W_mat) == 0L) {
    stop(".mean_glm_gradient_site_basis: cannot compute a mean gradient on zero observations.")
  }

  eta <- drop(cbind(1, W_mat) %*% alpha)
  h_prime <- switch(
    as.character(link_int),
    "0" = rep(1, length(eta)),
    "1" = {
      eta_clipped <- pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP)
      mu <- 1 / (1 + exp(-eta_clipped))
      mu * (1 - mu)
    },
    stop(sprintf(".mean_glm_gradient_site_basis: unsupported link_int=%s.", link_int))
  )

  colMeans(cbind(1, Z_mat) * as.numeric(abs(h_prime)))
}

# Potential-outcome means are calculated through the fitted outcome models
# below. This section intentionally has no standalone documented function.
# ============================================================================
# C++ ACCELERATED FUNCTIONS
# ============================================================================

# C++ functions (fit_unified_density_ratio_cpp, etc.) are compiled and loaded
# automatically by the package build system via useDynLib(RoCE) in NAMESPACE.
# The src/ C++ files are compiled during R CMD INSTALL or devtools::load_all().

#' Fit INITIAL density ratio function (γ_init parameters)
#'
#' Implements the INITIAL loss function for density ratio estimation.
#' Following eq:gamma_init in main.tex:
#' \deqn{
#' \ell(\boldsymbol{\gamma}_{s_j,1}) = \left(\widetilde{E}_t[\boldsymbol{\phi}(\mathbf{X})]\right)^T \boldsymbol{\gamma}_{s_j,1}
#'   + \widetilde{E}_{s_j}\left[I(A=1) \exp(-\boldsymbol{\phi}(\mathbf{X})^T \boldsymbol{\gamma}_{s_j,1})\right]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}_{s_j,1}\|_1.
#' }
#'
#' NOTE: This does NOT include ψ'(α_init) term - only uses mean_phi from target
#'
#' @param Z_site Covariate matrix (n × p)
#' @param A Treatment indicator vector (n × 1)
#' @param mean_phi Mean of φ(X) from target site: Ẽ_t[φ(X)]
#' @param lambda L1 regularization parameter. \code{NULL} (default) selects
#'   \code{lambda.min} by CV; numeric values use a fixed lambda.
#' @param nlambda Number of candidates in the CV lambda path.
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param lambda_rule CV selection rule when \code{lambda = NULL}: \code{"min"}
#'   selects \code{lambda.min}; \code{"1se"} selects \code{lambda.1se}.
#' @param A_val Treatment arm, either 0 or 1.
#' @param warm_start Optional initial coefficient vector for the final
#'   optimization.
#' @return Vector of initial density ratio parameters
#' @export
#' @inheritParams fit_initial_outcome
fit_initial_density_ratio <- function(Z_site, A, mean_phi, lambda = NULL,
                                       max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                       A_val = 1L, warm_start = NULL,
                                       nlambda = LAMBDA_GRID_SIZE_STANDARD,
                                       lambda_rule = c("min", "1se"), cv_group_id = NULL) {
  A_val <- .validate_A_val(A_val, "fit_initial_density_ratio")
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "fit_initial_density_ratio")
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, nrow(Z_site), "fit_initial_density_ratio"
  )
  cv_started_at <- proc.time()[["elapsed"]]
  support <- .density_ratio_support_penalty_floor(
    Z_site, A, mean_phi, A_val, "fit_initial_density_ratio"
  )
  path_floor_applied <- FALSE

  # Handle lambda = NULL case: use CV selection
  if (is.null(lambda)) {
    n_cv_folds <- .nuisance_cv_fold_count(A, A_val, "fit_initial_density_ratio")
    lmax <- compute_lambda_max_initial_dr(Z_site, A, mean_phi, A_val)
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) LAMBDA_MIN_RATIO_LOW_DIM else LAMBDA_MIN_RATIO_HIGH_DIM
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio,
                                     nlambda = nlambda)
    lambda_grid <- .constrain_density_ratio_lambda_grid(lambda_grid, support)
    path_floor_applied <- isTRUE(attr(
      lambda_grid, "support_path_floor_applied"
    ))
    cv_result <- .call_nuisance_cv_with_groups(
      select_lambda_cv_initial_density_ratio_cpp,
      list(Z_site, A, mean_phi, lambda_grid, n_cv_folds, max_iter, tol, A_val),
      cv_group_id, A, A_val, n_cv_folds, "fit_initial_density_ratio"
    )
    lambda <- .select_nuisance_cv_lambda(cv_result, lambda_rule, "fit_initial_density_ratio")
  } else {
    lambda <- .validate_lambda_scalar(lambda, "fit_initial_density_ratio")
    attr(lambda, "lambda_min") <- lambda
    attr(lambda, "lambda_1se") <- lambda
    attr(lambda, "lambda_rule") <- "fixed"
    attr(lambda, "cv_invalid_fold_fits") <- 0L
    attr(lambda, "cv_invalid_lambdas") <- 0L
    attr(lambda, "cv_path_tail_skipped_fold_fits") <- 0L
  }
  lambda <- .apply_density_ratio_support_floor(
    lambda, Z_site, A, mean_phi, A_val, "fit_initial_density_ratio",
    support = support, path_floor_applied = path_floor_applied
  )
  cv_seconds <- as.numeric(proc.time()[["elapsed"]] - cv_started_at)
  
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  final_fit_started_at <- proc.time()[["elapsed"]]
  cpp_result <- fit_initial_density_ratio_cpp(Z_site, A, mean_phi, lambda, max_iter, tol, A_val, ws)
  final_fit_seconds <- as.numeric(
    proc.time()[["elapsed"]] - final_fit_started_at
  )
  result <- cpp_result$gamma
  result <- .attach_nuisance_lambda_attrs(result, lambda)
  .attach_nuisance_fit_attrs(
    result, cpp_result,
    cv_seconds = cv_seconds,
    final_fit_seconds = final_fit_seconds
  )
}

#' Fit UNIFIED density ratio function (γ parameters)
#'
#' Implements the density ratio estimation loss function.
#' 
#' \strong{Refined loss} (untruncated counterpart of
#' \code{eq:gamma_calibrated_loss}, \code{calibrated = FALSE}):
#' \deqn{
#' \ell(\boldsymbol{\gamma}) = \widetilde{E}_t[\nabla_\alpha \psi]^T \boldsymbol{\gamma}
#'   + \widetilde{E}_{s_j}[I(A=1) \exp(-\phi_{site}^T \boldsymbol{\gamma}) \psi'(\phi_{outcome}^T \boldsymbol{\alpha})]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}\|_1.
#' }
#'
#' \strong{Calibrated loss} (\code{eq:gamma_calibrated_loss} in
#' \code{main.tex}, \code{calibrated = TRUE}):
#' \deqn{
#' \ell_{\gamma,\text{cal}}(\boldsymbol{\gamma}) = \widetilde{E}_t[\nabla_\alpha \psi]^T \boldsymbol{\gamma}
#'   + \widetilde{E}_{s_j}[I(A=1) \exp(-\phi_{site}^T \boldsymbol{\gamma}) \psi'(\mathcal{T}(\phi_{outcome}^T \boldsymbol{\alpha}))]
#'   + \lambda_\gamma \|\boldsymbol{\gamma}\|_1.
#' }
#'
#' @param Z_site Site model features for gamma parameterization (n x p_site)
#' @param A Treatment indicator vector (n x 1)
#' @param mean_grad_psi Mean gradient of ψ from target site, i.e., Ẽ_t[∇_α ψ(ϕ(X); α̂_init)]
#' @param alpha_init Initial outcome model parameters
#' @param lambda L1 regularization parameter. \code{NULL} (default) selects
#'   \code{lambda.min} by CV; numeric values use a fixed lambda.
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @param calibrated Whether to use calibrated loss with truncation
#'   \eqn{\mathcal{T}(\cdot)} (\code{TRUE} for calibrated loss
#'   \code{eq:gamma_calibrated_loss}; \code{FALSE} for its untruncated refined
#'   counterpart). Replaces the previous \code{use_truncation} parameter name.
#' @param M_tau Truncation threshold (only used if calibrated = TRUE)
#' @param W_outcome Outcome model features for psi' computation (n x p_outcome)
#'        Must be provided explicitly and have the same number of rows as \code{Z_site}.
#' @param A_val Treatment arm, either 0 or 1.
#' @param family_int Integer GLM-family code.
#' @param link_int Integer link-function code.
#' @param warm_start Optional initial coefficient vector for the final
#'   optimization.
#' @param lambda_rule CV selection rule when \code{lambda = NULL}: \code{"min"}
#'   selects \code{lambda.min}; \code{"1se"} selects \code{lambda.1se}.
#' @param nlambda Number of candidates in the CV lambda path.
#' @return Vector of density ratio parameters
#' @export
#' @inheritParams fit_initial_outcome
fit_unified_density_ratio <- function(Z_site, A, mean_grad_psi, alpha_init,
                                      lambda = NULL, max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                      calibrated = FALSE, M_tau = M_TAU_DEFAULT,
                                      W_outcome = NULL, A_val = 1L,
                                      family_int = 1L, link_int = 1L,
                                      warm_start = NULL,
                                      nlambda = LAMBDA_GRID_SIZE_STANDARD,
                                      lambda_rule = c("min", "1se"), cv_group_id = NULL) {
  A_val <- .validate_A_val(A_val, "fit_unified_density_ratio")
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "fit_unified_density_ratio")
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, nrow(Z_site), "fit_unified_density_ratio"
  )
  cv_started_at <- proc.time()[["elapsed"]]

  if (is.null(W_outcome)) {
    stop("fit_unified_density_ratio: W_outcome must be provided explicitly.")
  }
  W_outcome_use <- as.matrix(W_outcome)
  if (nrow(W_outcome_use) != nrow(Z_site)) {
    stop(sprintf("fit_unified_density_ratio: row mismatch between Z_site (%d) and W_outcome (%d).",
                 nrow(Z_site), nrow(W_outcome_use)))
  }
  if (length(mean_grad_psi) != ncol(as.matrix(Z_site)) + 1L) {
    stop(
      "fit_unified_density_ratio: mean_grad_psi length must match ",
      "ncol(Z_site) + 1.",
      call. = FALSE
    )
  }
  support <- .density_ratio_support_penalty_floor(
    Z_site, A, mean_grad_psi, A_val, "fit_unified_density_ratio"
  )
  path_floor_applied <- FALSE
  
  # Handle lambda = NULL case
  if (is.null(lambda)) {
    n_cv_folds <- .nuisance_cv_fold_count(A, A_val, "fit_unified_density_ratio")
    lmax <- compute_lambda_max_refined_dr(Z_site, A, mean_grad_psi, alpha_init,
                                           A_val = A_val,
                                           family_int = family_int, link_int = link_int,
                                           W_outcome = W_outcome,
                                           calibrated = calibrated, M_tau = M_tau)
    lambda_min_ratio <- if (nrow(Z_site) > ncol(Z_site)) LAMBDA_MIN_RATIO_LOW_DIM else LAMBDA_MIN_RATIO_HIGH_DIM
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio,
                                     nlambda = nlambda)
    lambda_grid <- .constrain_density_ratio_lambda_grid(lambda_grid, support)
    path_floor_applied <- isTRUE(attr(
      lambda_grid, "support_path_floor_applied"
    ))
    if (calibrated) {
      # Forward W_outcome to CV so ψ' uses correct features
      cv_result <- .call_nuisance_cv_with_groups(
        select_lambda_cv_calibrated_density_ratio_cpp,
        list(Z_site, A, mean_grad_psi, alpha_init, lambda_grid, n_cv_folds, max_iter, tol, M_tau,
             W_outcome_use, A_val, family_int, link_int),
        cv_group_id, A, A_val, n_cv_folds, "fit_unified_density_ratio"
      )
    } else {
      # Refined DR should use ψ'(W_outcome^T alpha_init) with no truncation.
      # Reuse the calibrated CV kernel with M_tau = Inf to keep feature usage
      # (Z_site for gamma, W_outcome for ψ') consistent with final fitting.
      cv_result <- .call_nuisance_cv_with_groups(
        select_lambda_cv_calibrated_density_ratio_cpp,
        list(Z_site, A, mean_grad_psi, alpha_init, lambda_grid, n_cv_folds, max_iter, tol, Inf,
             W_outcome_use, A_val, family_int, link_int),
        cv_group_id, A, A_val, n_cv_folds, "fit_unified_density_ratio"
      )
    }
    lambda <- .select_nuisance_cv_lambda(cv_result, lambda_rule, "fit_unified_density_ratio")
  } else {
    lambda <- .validate_lambda_scalar(lambda, "fit_unified_density_ratio")
    attr(lambda, "lambda_min") <- lambda
    attr(lambda, "lambda_1se") <- lambda
    attr(lambda, "lambda_rule") <- "fixed"
    attr(lambda, "cv_invalid_fold_fits") <- 0L
    attr(lambda, "cv_invalid_lambdas") <- 0L
    attr(lambda, "cv_path_tail_skipped_fold_fits") <- 0L
  }
  lambda <- .apply_density_ratio_support_floor(
    lambda, Z_site, A, mean_grad_psi, A_val,
    "fit_unified_density_ratio", support = support,
    path_floor_applied = path_floor_applied
  )
  cv_seconds <- as.numeric(proc.time()[["elapsed"]] - cv_started_at)
  
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  # C++ param name: mean_grad_psi (positional), calibrated (positional)
  final_fit_started_at <- proc.time()[["elapsed"]]
  cpp_result <- fit_unified_density_ratio_cpp(
    Z_site, A, mean_grad_psi, alpha_init, lambda, max_iter, tol, 
    calibrated, M_tau, W_outcome_use, A_val, family_int, link_int, ws
  )
  final_fit_seconds <- as.numeric(
    proc.time()[["elapsed"]] - final_fit_started_at
  )
  result <- cpp_result$gamma
  result <- .attach_nuisance_lambda_attrs(result, lambda)
  .attach_nuisance_fit_attrs(
    result, cpp_result,
    cv_seconds = cv_seconds,
    final_fit_seconds = final_fit_seconds
  )
}

#' Fit UNIFIED outcome model function (α parameters)
#'
#' Implements the weighted GLM loss for outcome regression following main.tex.
#' Supports GLM families gaussian and binomial.
#'
#' \strong{Refined loss} (untruncated counterpart of
#' \code{eq:alpha_calibrated_loss}, \code{calibrated = FALSE}):
#' \deqn{
#' \ell(\boldsymbol{\alpha}) = \widetilde{E}_{s_j}\left[\frac{I(A=1)}{\exp(\boldsymbol{\phi}^T \hat{\boldsymbol{\gamma}})} 
#'   (-Y \cdot \eta + \Psi(\eta))\right] + \lambda \|\boldsymbol{\alpha}\|_1.
#' }
#'
#' \strong{Calibrated loss} (\code{eq:alpha_calibrated_loss} in
#' \code{main.tex}, \code{calibrated = TRUE}):
#' \deqn{
#' \ell(\boldsymbol{\alpha}) = \widetilde{E}_{s_j}\left[\frac{I(A=1)}{\exp(\mathcal{T}(\boldsymbol{\phi}^T \hat{\boldsymbol{\gamma}}))} 
#'   (-Y \cdot \eta + \Psi(\eta))\right] + \lambda \|\boldsymbol{\alpha}\|_1.
#' }
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
#' @param calibrated Whether to use calibrated loss with truncation
#'   \eqn{\mathcal{T}(\cdot)} (\code{TRUE} for calibrated loss
#'   \code{eq:alpha_calibrated_loss}; \code{FALSE} for its untruncated refined
#'   counterpart).
#' @param M_tau Truncation threshold (only used if calibrated = TRUE)
#' @param Z_site Feature matrix for density ratio (Z_site). Must be provided explicitly.
#' @param lambda_rule CV selection rule when \code{lambda = NULL}: \code{"min"}
#'   selects \code{lambda.min}; \code{"1se"} selects \code{lambda.1se}.
#' @param nlambda Number of candidates in the CV lambda path.
#' @param warm_start Optional initial coefficient vector for the final
#'   optimization.
#' @return Vector of outcome model parameters
#' @export
#' @inheritParams fit_initial_outcome
fit_unified_outcome <- function(W_outcome, Y, A, A_val = 1, gamma_s, lambda = NULL,
                                max_iter = MAX_ITER_DEFAULT, tol = TOL_DEFAULT,
                                family = "binomial", link = NULL,
                                family_int = NULL, link_int = NULL,
                                calibrated = FALSE, M_tau = M_TAU_DEFAULT, Z_site = NULL,
                                warm_start = NULL,
                                nlambda = LAMBDA_GRID_SIZE_STANDARD,
                                lambda_rule = c("min", "1se"), cv_group_id = NULL) {
  A_val <- .validate_A_val(A_val, "fit_unified_outcome")
  lambda_rule <- .match_nuisance_lambda_rule(lambda_rule, "fit_unified_outcome")
  cv_group_id <- .validate_nuisance_cv_group_id(
    cv_group_id, nrow(W_outcome), "fit_unified_outcome"
  )
  cv_started_at <- proc.time()[["elapsed"]]
  
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
    n_cv_folds <- .nuisance_cv_fold_count(A, A_val, "fit_unified_outcome")
    lmax <- compute_lambda_max_outcome(W_outcome, Y, A, gamma_s,
                                        A_val = A_val,
                                        family_int = family_int, link_int = link_int,
                                        Z_site = Z_site,
                                        calibrated = calibrated, M_tau = M_tau)
    lambda_min_ratio <- if (nrow(W_outcome) > ncol(W_outcome)) LAMBDA_MIN_RATIO_LOW_DIM else LAMBDA_MIN_RATIO_HIGH_DIM
    lambda_grid <- build_lambda_grid(lambda_max = lmax,
                                     lambda_min_ratio = lambda_min_ratio,
                                     nlambda = nlambda)
    
    # Use correct CV function based on calibrated flag
    # - Calibrated (calibrated=TRUE) uses select_lambda_cv_calibrated_outcome_cpp
    # - Refined (calibrated=FALSE) uses select_lambda_cv_general_refined_outcome_cpp
    if (calibrated) {
      # Calibrated outcome CV - uses truncation in loss function
      cv_result <- .call_nuisance_cv_with_groups(
        select_lambda_cv_calibrated_outcome_cpp,
        list(W_outcome, Y, A, gamma_s, lambda_grid, n_cv_folds, max_iter, tol, A_val, M_tau, Z_site_mat,
             family_int, link_int),
        cv_group_id, A, A_val, n_cv_folds, "fit_unified_outcome"
      )
    } else {
      # Refined outcome CV - standard weighted loss without truncation
      cv_result <- .call_nuisance_cv_with_groups(
        select_lambda_cv_general_refined_outcome_cpp,
        list(W_outcome, Y, A, gamma_s, family_int, link_int, lambda_grid, n_cv_folds, max_iter, tol, A_val, Z_site_mat),
        cv_group_id, A, A_val, n_cv_folds, "fit_unified_outcome"
      )
    }
    lambda <- .select_nuisance_cv_lambda(cv_result, lambda_rule, "fit_unified_outcome")
  } else {
    lambda <- .validate_lambda_scalar(lambda, "fit_unified_outcome")
    attr(lambda, "lambda_min") <- lambda
    attr(lambda, "lambda_1se") <- lambda
    attr(lambda, "lambda_rule") <- "fixed"
    attr(lambda, "cv_invalid_fold_fits") <- 0L
    attr(lambda, "cv_invalid_lambdas") <- 0L
    attr(lambda, "cv_path_tail_skipped_fold_fits") <- 0L
  }
  cv_seconds <- as.numeric(proc.time()[["elapsed"]] - cv_started_at)
  
  # Use unified C++ function with all parameters
  # C++ param name: use_truncation (positional)
  ws <- if (!is.null(warm_start)) as.numeric(warm_start) else numeric(0)
  final_fit_started_at <- proc.time()[["elapsed"]]
  cpp_result <- fit_unified_outcome_cpp(W_outcome, Y, A, gamma_s, family_int, link_int,
                                        lambda, max_iter, tol, A_val, 
                                        calibrated, M_tau, Z_site_mat, ws)
  final_fit_seconds <- as.numeric(
    proc.time()[["elapsed"]] - final_fit_started_at
  )
  # C++ returns 'alpha' matching paper notation (main.tex convention)
  result <- cpp_result$alpha
  result <- .attach_nuisance_lambda_attrs(result, lambda)
  .attach_nuisance_fit_attrs(
    result, cpp_result,
    cv_seconds = cv_seconds,
    final_fit_seconds = final_fit_seconds
  )
}

# ============================================================================
# Shared helpers for weight optimization and variance aggregation
# ============================================================================

#' Convert C_cross to a matrix suitable for C++ functions
#'
#' Validates and converts a cross-site covariance specification to a K x K
#' matrix, or returns an empty matrix when it is explicitly \code{NULL}.
#' Non-null malformed inputs fail fast rather than silently dropping a
#' covariance term from the aggregation objective.
#'
#' @param C_cross Optional K x K cross-site covariance matrix, or NULL.
#' @param K Integer, expected dimension.
#' @return A K x K matrix (valid input) or a 0 x 0 matrix (\code{NULL}).
#' @keywords internal
prepare_cross_matrix <- function(C_cross, K) {
  if (is.null(C_cross)) {
    return(matrix(0, nrow = 0, ncol = 0))
  }
  if (!is.matrix(C_cross) || !is.numeric(C_cross)) {
    stop("prepare_cross_matrix: C_cross must be NULL or a numeric matrix.",
         call. = FALSE)
  }
  if (!identical(dim(C_cross), c(as.integer(K), as.integer(K)))) {
    stop(sprintf(
      "prepare_cross_matrix: C_cross must have dimension %d x %d; found %s.",
      K, K, paste(dim(C_cross), collapse = " x ")
    ), call. = FALSE)
  }
  if (any(!is.finite(C_cross))) {
    stop("prepare_cross_matrix: C_cross must contain only finite values.",
         call. = FALSE)
  }

  # Defensive normalization:
  # - Cross-site covariance is an off-diagonal object by design; diagonal is
  #   accounted for by V_t.
  # - Symmetrize to guard against small numerical asymmetries.
  C_cross <- 0.5 * (C_cross + t(C_cross))
  diag(C_cross) <- 0
  C_cross
}

#' Optimize aggregation weights with cross-site covariances
#'
#' Solves the convex optimization problem for optimal aggregation weights
#' following eq:agg_penalized_objective in main.tex.
#'
#' The optimization problem in \code{main.tex} is
#' \deqn{
#'   \min_{\boldsymbol{\eta}} \left\{
#'     N_{\mathrm{all}}\widehat{\mathrm{Var}}(\boldsymbol{\eta})
#'     + \sum_{j=1}^K
#'       \bigl(\lambda\widehat t_j-1\bigr)_+|\eta_j|
#'   \right\}.
#' }
#' where \eqn{\widehat{\mathrm{Var}}(\boldsymbol{\eta})} contains the
#' target-only, target-side source, source-side, target--source, and
#' cross-source covariance terms displayed in \code{main.tex}, and
#' \eqn{\widehat t_j} is the standardized target--source discrepancy. For the
#' primary analysis, \code{estimates} and \code{mu_ot} are TATE
#' estimates, so the treated--control covariance is already included in every
#' variance component.
#'
#' Its components are the target-only variance, source-specific target- and
#' source-side variances, target--source covariances, cross-source covariances,
#' and the truncated-Wald weighted L1 penalty. Ignoring the cross-source terms
#' \eqn{\widehat C_{t,s_j,s_k}} generally understates uncertainty.
#'
#' @param estimates Vector of source-assisted estimates (TATE in the
#'   primary analysis; length K).
#' @param variances List with \eqn{\widehat V_{ot}}, \eqn{\widehat V_t}
#'   (vector), and \eqn{\widehat V_s} (vector).
#' @param C_ot Vector of covariances \eqn{\widehat C_{ot,s_j}} (length K).
#' @param n_samples List with \eqn{N_t} and \eqn{N_{s_j}} (vector).
#' @param lambda Wald cutoff factor \eqn{\lambda} for source selection. The
#'   optimized objective is
#'   \eqn{N_{\mathrm{all}}\widehat{\mathrm{Var}}(\eta)
#'   +\sum_j(\lambda\widehat t_j-1)_+|\eta_j|}, matching
#'   \code{main.tex}.
#' @param mu_ot Target-only estimate on the same estimand scale as
#'   \code{estimates}.
#' @param C_cross Optional K by K cross-site covariance matrix
#'   \eqn{\widehat C_{t,s_j,s_k}}.
#' @param clip_weights Whether to clip weights to [0,1]. Default FALSE to match
#'        the unconstrained optimization problem in the paper. Set TRUE only as a
#'        pragmatic finite-sample safeguard (it changes the theoretical optimum).
#' @return Optimal weights
#'   \eqn{\boldsymbol{\eta}^* = (\eta_1^*, \ldots, \eta_K^*)}.
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
  numeric_inputs <- c(
    estimates, variances$V_t, variances$V_s, variances$V_ot,
    C_ot, n_samples$n_t, n_samples$n_s, mu_ot
  )
  if (any(!is.finite(numeric_inputs))) {
    stop(
      "optimize_weights: estimates, variance components, covariances, sample sizes, and mu_ot must be finite.",
      call. = FALSE
    )
  }
  if (length(n_samples$n_t) != 1L || n_samples$n_t <= 0 ||
      any(n_samples$n_s <= 0)) {
    stop("optimize_weights: target and source sample sizes must be positive.",
         call. = FALSE)
  }
  if (!is.null(C_cross) &&
      (!is.numeric(C_cross) || any(!is.finite(C_cross)))) {
    stop("optimize_weights: C_cross must contain only finite numeric values.",
         call. = FALSE)
  }

  # Input validation and clipping with more conservative variance bounds
  estimates <- pmax(pmin(estimates, ESTIMATE_MAX), -ESTIMATE_MAX)
  variances$V_t <- pmax(variances$V_t, VARIANCE_MIN)
  variances$V_s <- pmax(variances$V_s, VARIANCE_MIN)
  variances$V_ot <- max(variances$V_ot, VARIANCE_MIN)
  C_ot <- pmax(pmin(C_ot, ESTIMATE_MAX), -ESTIMATE_MAX)
  lambda <- .validate_lambda_scalar(lambda, "optimize_weights", allow_zero = TRUE)
  if (lambda > LAMBDA_MAX) {
    stop(sprintf("optimize_weights: lambda must be <= %g.", LAMBDA_MAX),
         call. = FALSE)
  }
  mu_ot <- max(min(mu_ot, ESTIMATE_MAX), -ESTIMATE_MAX)
  
  # Prepare cross-site covariance matrix
  cross_matrix <- prepare_cross_matrix(C_cross, length(estimates))
  
  cpp_result <- optimize_weights_cpp(estimates, variances$V_t, variances$V_s, n_samples$n_s, 
                                    variances$V_ot, n_samples$n_t, C_ot, 
                                    lambda, mu_ot, WEIGHT_OPT_MAX_ITER, WEIGHT_OPT_TOL,
                                    cross_matrix, numeric(0))

  if (!isTRUE(cpp_result$converged)) {
    stop("optimize_weights: optimizer did not converge; cannot solve main.tex eq:agg_penalized_objective for aggregation weights.")
  }
  
  # Validate output
  weights <- cpp_result$weights
  bad_idx <- which(is.na(weights) | is.infinite(weights))
  if (length(bad_idx) > 0) {
    stop(sprintf(
      "optimize_weights: %d/%d optimized weights were NA/Inf (indices: %s); refusing to replace them silently.",
      length(bad_idx), length(weights), paste(bad_idx, collapse = ",")
    ))
  }
  
  # NOTE: The optimization problem (eq:agg_penalized_objective in main.tex) is unconstrained.
  # Clipping to [0,1] is a numerical safeguard but changes the theoretical optimum.
  if (isTRUE(clip_weights)) {
    n_below <- sum(weights < 0.0)
    n_above <- sum(weights > 1.0)
    weights <- pmax(pmin(weights, 1.0), 0.0)
    if (n_below + n_above > 0L) {
      warning(sprintf(
        "optimize_weights: clipped %d of %d aggregation weights to [0, 1] (below 0: %d, above 1: %d). The unconstrained eq:agg_penalized_objective optimum is outside the simplex; downstream variance based on the clipped weights may not match the theoretical optimum.",
        n_below + n_above, length(weights), n_below, n_above
      ), call. = FALSE)
    }
  }

  attr(weights, "optimizer_iterations") <- as.integer(cpp_result$iterations)
  attr(weights, "psd_ridge") <- as.numeric(cpp_result$psd_ridge)
  
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
  # objective (eq:agg_penalized_objective in main.tex) but NOT in variance estimation.
  return(calculate_aggregated_variance_cpp(eta, variances$V_t, variances$V_s,
                                          n_samples$n_s, variances$V_ot,
                                          n_samples$n_t, C_ot, cross_matrix,
                                          lambda, mu_ot))
}
