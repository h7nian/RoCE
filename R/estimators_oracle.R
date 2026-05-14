# estimators_oracle.R - Oracle DR estimator for benchmarking
#
# Uses true (known) outcome and density ratio parameters to form a
# DR estimator. Provides the semiparametric efficiency lower bound.

# =============================================================================
# ORACLE DR ESTIMATOR
# =============================================================================
# Uses true (known) outcome parameters and density ratio parameters
# to form a DR estimator.  Aggregation follows the same optimal-weight
# framework as the cross-fitting estimators (eq:final_opt in main.tex),
# with the full quadratic variance formula that accounts for:
#   (a) target-side variance V_{t,s_j} of each source-assisted estimate,
#   (b) covariance C_{ot,s_j} between target-only and source-assisted,
#   (c) cross-site covariance C_{t,s_j,s_k} from the shared target mean.
# =============================================================================

#' Oracle Doubly Robust Estimator
#'
#' Constructs a DR estimator using the TRUE generative parameters (alpha, gamma)
#' rather than estimated nuisance models. This gives the best possible DR
#' performance and serves as an efficiency benchmark.
#'
#' @param data_split Split data by site
#' @param alpha1_true True treatment outcome parameters (intercept + p)
#' @param gamma_params True site model parameters (list from DGP)
#' @param outcome_type "binary" or "continuous"
#' @param A_val Treatment value to estimate (default 1)
#' @param alpha0_true True control outcome parameters (required when A_val=0)
#' @param target_propensity_true Optional numeric vector of true target-site
#'   treatment propensities P(A=1|X) for target observations. When provided,
#'   oracle DR uses this vector (converted to P(A=A_val|X)); otherwise it falls
#'   back to the estimated target-site propensity model.
#' @param lambda_selection Lambda choice for aggregation weights:
#'   - "cv" (default): data-adaptive selection via \code{select_lambda_cv_crossfit}
#'   - numeric: fixed lambda value
#' @param lambda_grid Optional lambda grid passed to
#'   \code{select_lambda_cv_crossfit} when \code{lambda_selection = "cv"}.
#' @param lambda_rule Rule used when \code{lambda_selection = "cv"}:
#'   \code{"min"} selects the variance minimizer and \code{"1se"} selects the
#'   largest lambda within 5\% of the minimum.
#' @return List with estimate, variance, se, method
#' @export
estimate_oracle_dr <- function(data_split, alpha1_true, gamma_params,
                                outcome_type = "binary", A_val = 1L,
                                alpha0_true = NULL,
                                target_propensity_true = NULL,
                                lambda_selection = "cv",
                                lambda_grid = NULL,
                                lambda_rule = c("min", "1se")) {
  lambda_rule <- match.arg(lambda_rule)
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  K <- length(source_sites)

  # Select the correct alpha based on A_val
  alpha_true <- if (A_val == 1) alpha1_true else alpha0_true
  if (is.null(alpha_true)) {
    stop("estimate_oracle_dr: alpha0_true must be provided when A_val=0")
  }

  # ---------- Target predictions: m(X_i; alpha_true) ----------
  W_target <- as.matrix(target_data$W_outcome)
  W_target_design <- cbind(1, W_target)
  n_t <- target_data$n

  # Defensive: alpha dimension must match design matrix
  if (length(alpha_true) != ncol(W_target_design)) {
    stop(sprintf(
      "estimate_oracle_dr: alpha_true has length %d but W_outcome has %d columns (need %d = p+1 with intercept).",
      length(alpha_true), ncol(W_target), ncol(W_target_design)
    ))
  }

  if (outcome_type == "binary") {
    m_target <- logistic(as.numeric(W_target_design %*% alpha_true))
  } else {
    m_target <- as.numeric(W_target_design %*% alpha_true)
  }
  M_t <- mean(m_target)

  # ---------- Target-only AIPW influence function ----------
  ind_a <- as.numeric(target_data$A == A_val)
  used_true_propensity <- FALSE
  pi_a <- if (!is.null(target_propensity_true)) {
    if (length(target_propensity_true) != n_t) {
      stop(sprintf(
        "estimate_oracle_dr: target_propensity_true has length %d but target sample size is %d.",
        length(target_propensity_true), n_t
      ))
    }
    if (any(!is.finite(target_propensity_true))) {
      stop("estimate_oracle_dr: target_propensity_true contains non-finite values.")
    }
    pi1 <- pmax(pmin(as.numeric(target_propensity_true), POSITIVITY_UPPER), POSITIVITY_LOWER)
    used_true_propensity <- TRUE
    if (A_val == 1L) pi1 else pmax(pmin(1 - pi1, POSITIVITY_UPPER), POSITIVITY_LOWER)
  } else {
    tryCatch({
      Z_target <- as.matrix(target_data$Z_site)
      ps_fit <- stats::glm(ind_a ~ Z_target, family = stats::binomial())
      pmax(pmin(ps_fit$fitted.values, POSITIVITY_UPPER), POSITIVITY_LOWER)
    }, error = function(e) {
      stop(sprintf("estimate_oracle_dr: target propensity GLM failed: %s",
                   conditionMessage(e)), call. = FALSE)
    })
  }
  # varphi_ot_i = I(A=a)/pi_a * (Y - m(X)) + m(X) - M_t  (per-obs, length n_t)
  varphi_ot <- ind_a / pi_a * (target_data$Y - m_target) + m_target - M_t
  V_ot <- mean(varphi_ot^2)  # un-normalised; divided by n_t in the quadratic form

  # ---------- Target-side IF pieces used for variance/covariance ----------
  # For each source site j, define zeta_{t,s_j,i} on the target sample.
  # In the default oracle setting alpha_true is shared, so these columns are
  # often identical; we still estimate C_cross from the sample covariance
  # matrix for consistency with the cross-fitting implementation.
  varphi_ot_centered <- varphi_ot - mean(varphi_ot)

  # ---------- Source site DR corrections ----------
  source_delta <- numeric(0)   # delta_s for each source
  source_V_s <- numeric(0)     # V_{s_j} for each source
  source_n_s <- numeric(0)     # n_{s_j}
  source_V_t <- numeric(0)     # V_{t,s_j} for each source
  C_ot_vec <- numeric(0)       # C_{ot,s_j} for each source
  zeta_target_list <- list()   # zeta_{t,s_j} vectors on target sample
  valid_sites <- character(0)

  for (site in source_sites) {
    source_data <- data_split[[site]]
    Z_source <- as.matrix(source_data$Z_site)
    W_source <- as.matrix(source_data$W_outcome)
    W_source_design <- cbind(1, W_source)
    Z_source_int <- cbind(1, Z_source)

    y_source <- source_data$Y
    a_source <- source_data$A
    n_source <- source_data$n

    # True outcome model predictions on source
    if (outcome_type == "binary") {
      m_source <- logistic(as.numeric(W_source_design %*% alpha_true))
    } else {
      m_source <- as.numeric(W_source_design %*% alpha_true)
    }

    # True density ratio weights
    site_num <- as.integer(gsub("s", "", site))
    gamma_key <- paste0("s", site_num, "_", A_val)
    gamma_true <- gamma_params[[gamma_key]]

    if (is.null(gamma_true)) {
      stop(sprintf("estimate_oracle_dr: gamma for %s not found; cannot compute oracle source-assisted estimate.",
                   gamma_key), call. = FALSE)
    }

    eta_dr <- as.numeric(Z_source_int %*% gamma_true)
    w_i <- exp(-eta_dr)

    # DR correction
    treated_mask <- as.numeric(a_source == A_val)
    delta_s <- mean(treated_mask * w_i * (y_source - m_source))

    # Source-side IF variance: V_{s_j} = Var_s(xi_i)
    xi_i <- treated_mask * w_i * (y_source - m_source) - delta_s
    V_s <- mean(xi_i^2)

    # Target-side component for this site: zeta_{t,s_j,i} = m_{t,s_j}(X_i) - M_{t,s_j}
    # Oracle uses true alpha; this is typically common across sites, but we
    # keep a per-site vector so C_cross is always sample-estimated.
    m_target_site <- if (outcome_type == "binary") {
      logistic(as.numeric(W_target_design %*% alpha_true))
    } else {
      as.numeric(W_target_design %*% alpha_true)
    }
    M_t_site <- mean(m_target_site)
    zeta_site <- m_target_site - M_t_site
    V_t_site <- mean(zeta_site^2)
    C_ot_site <- mean(varphi_ot_centered * (zeta_site - mean(zeta_site)))

    valid_sites <- c(valid_sites, site)
    source_delta <- c(source_delta, delta_s)
    source_V_s   <- c(source_V_s, max(V_s, VARIANCE_MIN))
    source_n_s   <- c(source_n_s, n_source)
    source_V_t   <- c(source_V_t, max(V_t_site, VARIANCE_MIN))
    C_ot_vec     <- c(C_ot_vec, C_ot_site)
    zeta_target_list[[length(zeta_target_list) + 1L]] <- zeta_site
  }

  K_valid <- length(valid_sites)

  if (K_valid == 0) {
    stop("estimate_oracle_dr: no valid source sites remained; refusing to silently return target-only oracle estimate.",
         call. = FALSE)
  }

  # ---------- Construct variance components (same interface as cross-fit) ---
  # Source-assisted estimates: mu_{t,s_j} = M_t + delta_s_j
  mu_ts <- M_t + source_delta  # length K_valid

  variances <- list(
    V_ot = max(V_ot, VARIANCE_MIN),
    V_t  = source_V_t,
    V_s  = source_V_s
  )
  n_samples <- list(n_t = n_t, n_s = source_n_s)

  # Cross-site covariance matrix C_cross estimated from target-sample zeta columns
  if (K_valid > 1) {
    zeta_matrix <- do.call(cbind, zeta_target_list)
    zeta_centered <- scale(zeta_matrix, center = TRUE, scale = FALSE)
    C_cross <- crossprod(zeta_centered) / n_t
    diag(C_cross) <- 0  # diagonal handled by V_t
  } else {
    C_cross <- matrix(0, nrow = K_valid, ncol = K_valid)
  }

  lambda_reg <- if (is.character(lambda_selection)) {
    if (!identical(lambda_selection, "cv")) {
      stop("estimate_oracle_dr: lambda_selection must be 'cv' or numeric.")
    }
    select_lambda_cv_crossfit(
      n_folds = 1,
      fold_target_estimate = M_t,
      fold_source_estimates = mu_ts,
      variances_k1 = variances,
      C_ot_k1 = C_ot_vec,
      n_samples_k1 = n_samples,
      C_cross_k1 = C_cross,
      verbose = FALSE,
      crossfit_type = "oracle",
      lambda_grid = lambda_grid,
      lambda_rule = lambda_rule
    )
  } else if (is.numeric(lambda_selection) && length(lambda_selection) == 1L && is.finite(lambda_selection)) {
    lambda_selection
  } else {
    stop("estimate_oracle_dr: lambda_selection must be 'cv' or a finite numeric scalar.")
  }

  # ---------- Optimal weights and aggregated variance ----------
  eta <- optimize_weights(
    estimates    = mu_ts,
    variances    = variances,
    C_ot         = C_ot_vec,
    n_samples    = n_samples,
    lambda       = lambda_reg,
    mu_ot     = M_t,
    C_cross      = C_cross,
    clip_weights = FALSE
  )

  final_estimate <- calculate_aggregated_estimate_cpp(M_t, mu_ts, eta)
  final_variance <- calculate_aggregated_variance(
    eta, variances, C_ot_vec, n_samples, C_cross
  )

  return(list(
    estimate = as.numeric(final_estimate),
    variance = as.numeric(final_variance),
    se       = as.numeric(sqrt(final_variance)),
    method   = "oracle_dr",
    oracle_propensity = if (used_true_propensity) "true" else "estimated",
    lambda_used = as.numeric(lambda_reg),
    aggregation_lambda_rule = lambda_rule,
    weights  = as.numeric(eta),
    n        = n_t + sum(source_n_s)
  ))
}
