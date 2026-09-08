# constants.R - Centralized numerical constants for RoCE algorithms
#
# This file provides a single source of truth for all numerical constants
# used across the codebase. All other files should reference these constants
# instead of defining their own.
#
# =============================================================================
# NAMING CONVENTION NOTE
# =============================================================================
# This project follows snake_case for functions and variables. However,
# variables that directly correspond to mathematical notation in the
# accompanying paper (main.tex / supplemental.tex) preserve their capitalization
# for traceability:
#
#   mu_hat_ot -> target-only estimate hat{mu}^1_{ot} in paper (was Delta_ot)
#   V_ot      -> target-only variance
#   V_t, V_s  -> target / source variance components
#   C_ot      -> covariance between target-only and source-assisted
#   C_cross   -> cross-site covariance matrix
#   mu_pred_ts -> E_t[psi(phi(X); alpha_{t,s_j})] outcome model prediction on target (was M_ts)
#   M_tau     -> truncation threshold
#   K         -> number of source sites
#   varphi_ot    -> CENTERED influence function: hat{varphi}_{ot,i} - hat{mu}^1_{ot}
#               (eq:if_target_only in main.tex, after subtracting the estimate)
#
# Single-letter data variables (A, Y, X, Z, W, R, K, p, n) follow standard
# statistical notation and are acceptable in this context.
#
# =============================================================================
# OPTIMIZATION DEFAULTS
# =============================================================================

#' Default truncation radius M_tau for the calibrated estimation. Applied as a
#' SINGLE radius to BOTH the calibrated fitting losses (eq:gamma_calibrated_loss,
#' eq:alpha_calibrated_loss in main.tex) AND the inference weight (see
#' M_TAU_INFERENCE_DEFAULT), mirroring SMMAL (hou2025efficient), which truncates
#' the linear predictor at 2M with M the logit-scale positivity bound (M = 2.2 for
#' probabilities in [0.1, 0.9], i.e. radius 2M ~= 4.4). We use 5 as a
#' finite-sample safeguard on the FACE DGP: observed fitted logits are usually
#' well inside this radius, while rare excursions on the high-dimensional
#' [X, X^2] design are capped (weight cap e^5 ~= 148, versus e^10 ~= 22000).
#' A fixed radius is not, by itself, an asymptotic-inactivity guarantee;
#' production summaries therefore report clipping frequencies and include an
#' M_tau sensitivity analysis. The theory instead permits a radius sequence
#' whose clipping effect is o_p(N^{-1/2}). See the truncation discussion in
#' main.tex and supplemental.tex.
M_TAU_DEFAULT <- 5.0

#' Truncation radius for inference (correction term and source variance). Set
#' equal to M_TAU_DEFAULT so ONE radius governs fitting and inference, as in SMMAL
#' (hou2025efficient): the truncated weight w_T = exp(-T(phi^T gamma)) feeds the
#' final estimator AND its variance (SMMAL eq. 13 -> eq. 6/7). On the high-dim
#' [X, X^2] density-ratio design a few source units can otherwise receive an
#' extreme untruncated weight exp(-phi^T gamma), producing heavy-tailed estimates
#' (the per-sim SE co-moves, so coverage holds, but the marginal empirical SE is
#' inflated ~2x). Truncating the inference linear predictor at |g| <= 5 collapses
#' the empirical SE onto the model SE (verified) while leaving the bulk -- and
#' hence the estimand -- unchanged.
M_TAU_INFERENCE_DEFAULT <- 5.0

#' Legacy arm-wise truncation parameter for the RHC application.
#' \code{run_rhc_experiment()} retains this historical radius for reproducible
#' secondary arm-specific summaries. The primary TATE manuscript call
#' uses \code{M_tau_inference = 5}, matching \code{M_TAU_INFERENCE_DEFAULT} and
#' the value stated in \code{main.tex}.
M_TAU_INFERENCE_RHC <- 3.0

#' Maximum iterations for coordinate descent optimization
MAX_ITER_DEFAULT <- 10000L

#' Convergence tolerance for optimization
TOL_DEFAULT <- 1e-6

#' Lambda grid size for standard CV selection
LAMBDA_GRID_SIZE_STANDARD <- 100

#' Maximum iterations for glmnet paths used by comparison and initial-outcome
#' nuisances. High-dimensional binomial folds can require more than glmnet's
#' default 100,000 iterations at the least-regularized tail even when the
#' selected lambda is well behaved. Warnings are intentionally not suppressed;
#' the larger cap reduces avoidable tail failures while preserving an explicit
#' signal if a path still does not converge.
GLMNET_MAX_ITER <- 1000000L

#' Number of consecutive nonconverged candidates tolerated within one fold of
#' a descending nuisance-parameter lambda path before the remaining, less
#' regularized tail is marked ineligible. Keep synchronized with the C++
#' `NumericalConstants::CV_FAILURE_PATIENCE` value.
NUISANCE_CV_FAILURE_PATIENCE <- 3L

#' Ratio of \code{lambda_min} to \code{lambda_max} in the CV grid when the
#' design is low-dimensional (n > p). Matches glmnet's default in that regime.
LAMBDA_MIN_RATIO_LOW_DIM <- 1e-4

#' Ratio of \code{lambda_min} to \code{lambda_max} in the CV grid when the
#' design is high-dimensional (n <= p). Matches glmnet's default in that regime.
LAMBDA_MIN_RATIO_HIGH_DIM <- 0.01

# =============================================================================
# VARIANCE / WEIGHT OPTIMIZATION BOUNDS
# =============================================================================

#' Minimum variance floor (prevent division by zero)
VARIANCE_MIN <- 1e-6

#' Maximum absolute value for estimates (clipping bound)
ESTIMATE_MAX <- 1e6

#' Minimum lambda for weight optimization
LAMBDA_MIN <- 1e-6

#' Maximum lambda for weight optimization.
#' This is a numerical sanity cap for user-supplied sensitivity grids; the
#' production aggregation path uses the fixed \code{AGG_WALD_LAMBDA}.
LAMBDA_MAX <- 1e6

#' Aggregation Wald cutoff constants
#'
#' The FACE truncated-Wald penalty factor uses lambda as the reciprocal of the
#' penalty-activation cutoff: a source is unpenalized when its discrepancy
#' statistic is no larger than \eqn{1/\lambda}. Crossing the cutoff initiates
#' soft shrinkage and does not guarantee an exactly zero finite-sample weight.
#' The coverage-blind disjoint pilot documented in the manuscript selected
#' cutoff \eqn{c=1}, hence multiplier \eqn{\lambda=1}.
#'
#' \code{AGG_WALD_CUTOFF_DEFAULT} records that locked production cutoff.
#' \code{AGG_WALD_LAMBDA} reads the optional build-time environment override
#' \code{ROCE_AGG_WALD_LAMBDA}, while \code{AGG_WALD_CUTOFF} always records
#' its exact reciprocal for the installed build. A near-zero positive override
#' approximately disables adaptive source penalization and is intended only for
#' explicit sensitivity analysis.
#'
#' @name aggregation_wald_defaults
#' @aliases AGG_WALD_CUTOFF_DEFAULT AGG_WALD_LAMBDA AGG_WALD_CUTOFF
#' @return A positive numeric scalar.
NULL

#' @rdname aggregation_wald_defaults
AGG_WALD_CUTOFF_DEFAULT <- 1.0

#' @rdname aggregation_wald_defaults
AGG_WALD_LAMBDA <- local({
  .val <- suppressWarnings(as.numeric(Sys.getenv(
    "ROCE_AGG_WALD_LAMBDA",
    format(1 / AGG_WALD_CUTOFF_DEFAULT, scientific = FALSE)
  )))
  if (length(.val) != 1L || is.na(.val) || !is.finite(.val) || .val <= 0) {
    stop("AGG_WALD_LAMBDA: ROCE_AGG_WALD_LAMBDA must be a single positive finite number.",
         call. = FALSE)
  }
  .val
})

#' @rdname aggregation_wald_defaults
AGG_WALD_CUTOFF <- 1 / AGG_WALD_LAMBDA

#' Maximum iterations for weight optimization (optimize_weights)
WEIGHT_OPT_MAX_ITER <- 10000L

#' Convergence tolerance for weight optimization (optimize_weights)
WEIGHT_OPT_TOL <- 1e-8

#' Default number of multiplier (wild) bootstrap replicates for comparison-method
#' standard errors (\code{.multiplier_bootstrap_se}). Large enough to keep the
#' Monte-Carlo error of the bootstrap SE small relative to its sampling variability.
BOOTSTRAP_REPLICATES_DEFAULT <- 5000L

# =============================================================================
# SITE BALANCE
# =============================================================================

#' Base intercept for balanced site allocation.
#' In the multinomial logistic site model, target is the baseline category
#' with all-zero coefficients.  Setting a negative intercept on each source
#' category's gamma reduces the source probabilities so that
#' P(target) ≈ 1/(K+1) under the baseline setting.
#'
#' NOTE: shift_strength is intended to control the *covariate-dependent* part
#' of the site model (the ε components), not the intercept. Scaling both the
#' intercept and covariate terms would be cancelled by the unit-norm
#' normalization in data generation.
#'
#' Derived analytically (unit-norm + Gaussian covariates) and calibrated
#' empirically; gives approximately 25\% target for K=3, 20\% for K=4, etc.
GAMMA_BALANCE_INTERCEPT <- -0.15

#' Default covariate-shift strength for the RoCE DGP.
#' A moderate default keeps the main simulation away from the strong-transport
#' stress-test regime while preserving nontrivial source/target covariate shift.
ROCE_SHIFT_STRENGTH_DEFAULT <- 0.5

# =============================================================================
# NUMERICAL STABILITY
# =============================================================================

#' Default epsilon for numerical stability operations
EPSILON_DEFAULT <- 1e-8

#' Division floor to prevent divide-by-zero (denominator clipping)
#' @export
DIVISION_FLOOR <- 1e-10

#' Clipping bounds for logistic/eta.
#' Alias for C++ NumericalConstants::ETA_CLIP_MIN/MAX (== 50).
#' Used in R-side logistic/GLM functions. The C++ side uses the same value
#' via NumericalConstants::ETA_CLIP_MIN/MAX in numerical_constants.h.
LOGISTIC_CLIP <- 50

# =============================================================================
# GLM FAMILY/LINK INTEGER CODES
# =============================================================================
# Named constants for GLM integer codes, matching C++ enums in utils.h.
# Use these instead of raw integers for type safety and readability.
# See resolve_glm_family() for string ↔ integer mapping.

#' Integer codes for supported GLM families and links
#'
#' These constants mirror the enumerations used by the compiled nuisance-model
#' kernels.
#' @name glm_integer_codes
#' @aliases FAMILY_GAUSSIAN FAMILY_BINOMIAL LINK_IDENTITY LINK_LOGIT
#' @return Integer scalar.
NULL

#' @rdname glm_integer_codes
#' @export
FAMILY_GAUSSIAN <- 0L
#' @rdname glm_integer_codes
#' @export
FAMILY_BINOMIAL <- 1L

#' @rdname glm_integer_codes
#' @export
LINK_IDENTITY <- 0L
#' @rdname glm_integer_codes
#' @export
LINK_LOGIT    <- 1L

# =============================================================================
# C++ MATCHING CONSTANTS
# =============================================================================
# These constants mirror C++ NumericalConstants in numerical_constants.h.
# Keep both files in sync when changing values.

#' Default number of folds for the outer cross-fitting loop (K_f in paper)
N_FOLDS_DEFAULT <- 10

#' Number of CV folds for internal lambda selection in C++ nuisance models.
#' Separate from N_FOLDS_DEFAULT (outer cross-fitting) to allow independent tuning.
N_CV_FOLDS_LAMBDA <- 5L

#' Density ratio weight bounds (matches C++ WEIGHT_MIN/MAX)
WEIGHT_MIN <- 1e-10
WEIGHT_MAX <- 1e10

#' Parameter clipping bound (matches C++ PARAM_MAX)
PARAM_MAX <- 100.0

#' Active set threshold (matches C++ ACTIVE_SET_THRESHOLD)
ACTIVE_SET_THRESHOLD <- 1e-10

#' Hessian floor (matches C++ HESSIAN_FLOOR)
HESSIAN_FLOOR <- 1e-8

# =============================================================================
# COMPARISON METHODS CONSTANTS
# =============================================================================

#' Propensity score clipping bounds (prevent extreme IPW weights)
PROP_SCORE_LOWER <- 0.01
PROP_SCORE_UPPER <- 0.99

#' Outcome prediction clipping bounds (prevent log(0) in binomial deviance)
OUTCOME_PRED_LOWER <- 0.001
OUTCOME_PRED_UPPER <- 0.999

#' Ridge regularization for matrix inversion
RIDGE_DEFAULT <- 1e-6

#' Density ratio weight clipping bounds
DR_WEIGHT_LOWER <- 0.1
DR_WEIGHT_UPPER <- 10.0

# =============================================================================
# CONFIDENCE INTERVAL
# =============================================================================

#' Z-score for 95\% confidence intervals
#' @export
Z_ALPHA_05 <- qnorm(0.975)  # 1.959964...

# =============================================================================
# DATA GENERATION / POSITIVITY
# =============================================================================

#' Positivity lower bound for propensity score clipping in data generation
POSITIVITY_LOWER <- 0.1

#' Positivity upper bound for propensity score clipping in data generation
POSITIVITY_UPPER <- 0.9

#' Default noise standard deviation for continuous outcomes in RoCE DGP
ROCE_NOISE_SD_DEFAULT <- 0.1

#' Shift SD multiplier for independent site allocation covariate shifts
#' Controls magnitude of site-specific mean shifts: delta ~ N(0, INDEPENDENT_SHIFT_SD * shift_strength)
INDEPENDENT_SHIFT_SD <- 0.3

#' Minimum number of treated units required to fit an outcome model
#' Used across estimator functions; falls back to marginal mean when below threshold
MIN_TREATED_FOR_MODEL <- 5L

#' Checkpoint save interval: save every N parameter settings
CHECKPOINT_SETTING_INTERVAL <- 5L

#' Checkpoint save interval: save every N simulations within a setting
CHECKPOINT_SIM_INTERVAL <- 20L

# =============================================================================
# TRUE OUTCOME PARAMETERS (single source of truth)
# =============================================================================

#' Get true outcome model parameters for simulation
#'
#' Returns the (intercept + p covariate) parameter vectors for the treated
#' (alpha1) and control (alpha0) outcome models used in data generation.
#' This is the single source of truth — both \code{generate_simulation_data}
#' and \code{calculate_superpopulation_truth} call this function.
#'
#' @param p Number of covariates (excluding intercept).
#' @return Named list with \code{alpha1} and \code{alpha0} vectors of length p+1.
#' @export
get_true_outcome_parameters <- function(p) {
  # Treatment outcome: intercept + 2 non-zero covariate coefficients (X1, X2)
  # Control outcome: intercept + 2 non-zero covariate coefficients (X1, X2)
  # Sparsity matches the site model gamma (also 2 non-zero covariates)
  # Coefficients (0.5, 0.75) provide a difficulty gradient: the stronger
  # coefficient (0.75) is easier for LASSO to detect, while 0.5 requires
  # pooling across sites for reliable recovery (binary outcome, eff_n~60).
  alpha1 <- c(0.5, 0.5, 0.75, rep(0, max(0, p - 2)))
  alpha0 <- c(-0.5, 0.75, 0, rep(0, max(0, p - 2)))
  list(
    alpha1 = alpha1[1:(p + 1)],
    alpha0 = alpha0[1:(p + 1)]
  )
}

# =============================================================================
# GLM FAMILY UTILITIES
# =============================================================================

#' Valid GLM family names (glmnet-style API)
VALID_GLM_FAMILIES <- c("gaussian", "binomial")

# =============================================================================
# VALID OPTION VECTORS (single source of truth for string validation)
# =============================================================================
# All parameter-validation helpers should reference these vectors instead
# of hard-coding option lists.

#' Valid configuration identifiers
VALID_CONFIGS <- c("C1", "C2", "C3", "C4")

#' Valid estimand types
VALID_ESTIMAND_TYPES <- c("sample", "superpopulation")

#' Valid site allocation methods
VALID_SITE_ALLOCATIONS <- c("model", "uniform", "balanced", "target_heavy", "source_heavy", "independent")

#' Valid covariate transformation types
VALID_TRANSFORM_TYPES <- c("strong", "mild", "none")

#' Valid outcome types
VALID_OUTCOME_TYPES <- c("binary", "continuous")

#' Valid heterogeneity types
VALID_HETEROGENEITY_TYPES <- c("none", "mild", "strong", "partial")

#' Valid data generating process (DGP) types
#'
#' "roce"      – current RoCE DGP (multinomial logistic site model,
#'                nonlinear covariate transforms, binary or continuous outcome)
#' "face" – DGP adapted from Han et al. (JASA 2023, Section 5.1):
#'                site-specific skew-normal covariates, linear + squared
#'                outcome predictors, and site-specific constant treatment
#'                shifts for either continuous or binary outcomes
VALID_DGP_TYPES <- c("roce", "face")

# =============================================================================
# FACE PAPER DGP CONSTANTS (Han et al., JASA 2023, Section 5.1)
# =============================================================================

#' True ATE for the target population in FACE paper DGP, continuous outcome
#' (Δ_T = 3.0, the conditional mean shift between arms).
FACE_ATE_TARGET <- 3.0

#' True ATE for the target population in FACE paper DGP, binary outcome
#' (Δ_bin = 1.0, the conditional log-odds shift between arms). The continuous
#' Δ_T = 3.0 is a mean shift and saturates expit(); the binary outcome instead
#' uses a moderate log-odds shift so neither arm is pushed against 0/1. The
#' resulting target-population risk difference is ≈ 0.21 (not 1.0).
FACE_BINARY_ATE_TARGET <- 1.0

#' Target logit-scale standard deviation of the covariate signal for the binary
#' FACE paper DGP. The raw linear predictor (X−κ)ᵀβ_lin + (X²)ᵀβ_sq has a
#' target-law sd of ≈ 3.2, which saturates expit() (treated-arm prevalence
#' ≈ 0.98 → near-degenerate outcome regression). For the binary outcome the
#' signal is standardized to this controlled spread; see
#' \code{get_face_binary_calibration()}.
FACE_BINARY_SIGNAL_SD <- 1.0

#' Noise standard deviation for FACE paper outcomes: 2√5 ≈ 4.47
FACE_NOISE_SD <- 2 * sqrt(5)

#' Shared location parameter κ for the skewed-normal covariate distribution.
#' All sites use the same κ (midpoint of the (0.10, 0.15) range in the paper).
FACE_KAPPA <- 0.125

#' Maximum skewness parameter ν for source-site covariates (paper: ν ∈ [0, 0.2]).
#' Source sites receive equally-spaced ν values between ν_max/K and ν_max.
FACE_NU_SOURCE_MAX <- 0.2

#' Get true outcome model parameters for the FACE paper DGP
#'
#' Returns linear (β_lin) and squared (β_sq) covariate coefficients.
#' The first \code{min(4, p)} covariates have non-zero coefficients
#' equally spaced from 0.4 to 1.2 (consistent with the paper's p = 10 spacing);
#' all remaining covariates have zero coefficients.
#'
#' Both β_lin and β_sq use the same values because the paper sets
#' β_{1a} = β_{2a} (same coefficients for linear and squared blocks,
#' same coefficients for a = 0 and a = 1). The treatment effect is
#' captured solely by the site-specific ATE Δ_k.
#'
#' @param p Number of covariates (excluding intercept).
#' @return Named list with \code{beta_linear} and \code{beta_squared},
#'   each a numeric vector of length p.
#' @export
get_face_outcome_parameters <- function(p) {
  n_nonzero <- min(4L, p)
  # Equally-spaced from 0.4 to 1.2 (same spacing as the paper's 10-covariate setting)
  beta_nonzero <- seq(0.4, 1.2, length.out = n_nonzero)
  list(
    beta_linear  = c(beta_nonzero, rep(0, p - n_nonzero)),
    beta_squared = c(beta_nonzero, rep(0, p - n_nonzero))
  )
}

#' Get propensity-score model parameters for the FACE paper DGP
#'
#' Returns linear (α_1) and squared (α_2) covariate coefficients used in:
#'   π_k = expit(X α_1 + X² α_2)
#'
#' The first \code{min(4, p)} covariates have non-zero coefficients:
#' - α_1: equally-spaced from 0.5 to -0.5 (decreasing, as in the paper)
#' - α_2: (-0.5, 0, …, 0) for the same block
#'
#' @param p Number of covariates.
#' @return Named list with \code{alpha1} and \code{alpha2},
#'   each a numeric vector of length p.
#' @export
get_face_ps_parameters <- function(p) {
  n_nonzero <- min(4L, p)
  alpha1_nonzero <- seq(0.5, -0.5, length.out = n_nonzero)  # equally-spaced decrements
  alpha2_nonzero <- c(-0.5, rep(0, n_nonzero - 1L))
  list(
    alpha1 = c(alpha1_nonzero, rep(0, p - n_nonzero)),
    alpha2 = c(alpha2_nonzero, rep(0, p - n_nonzero))
  )
}

#' Resolve GLM family specification into all derived parameters
#'
#' Takes a glmnet-style family string and returns a list containing the
#' canonical link, integer codes for C++, and glmnet family name.
#' This centralizes the mapping so it doesn't need to be repeated.
#'
#' @param family Character string: one of "gaussian", "binomial".
#' @return Named list with elements:
#'   \itemize{
#'     \item \code{family}: normalized family name
#'     \item \code{link}: canonical link function name ("identity", "logit")
#'     \item \code{family_int}: integer code for C++ (0/1)
#'     \item \code{link_int}: integer code for C++ (0/1)
#'     \item \code{glmnet_family}: family string for \code{glmnet::cv.glmnet}
#'   }
#' @export
resolve_glm_family <- function(family) {
  family <- tolower(family)

  if (!family %in% VALID_GLM_FAMILIES) {
    stop(sprintf("Unsupported GLM family: '%s'. Must be one of: %s",
                 family, paste(VALID_GLM_FAMILIES, collapse = ", ")))
  }
  
  link <- switch(family,
    "gaussian" = "identity",
    "binomial" = "logit"
  )
  family_int <- switch(family,
    "gaussian" = FAMILY_GAUSSIAN,
    "binomial" = FAMILY_BINOMIAL
  )
  link_int <- switch(link,
    "identity" = LINK_IDENTITY,
    "logit" = LINK_LOGIT
  )
  glmnet_family <- switch(family,
    "gaussian" = "gaussian",
    "binomial" = "binomial"
  )
  
  list(
    family      = family,
    link        = link,
    family_int  = family_int,
    link_int    = link_int,
    glmnet_family = glmnet_family
  )
}
