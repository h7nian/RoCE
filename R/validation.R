# validation.R - Input validation functions for FACE-C algorithms
#
# This file contains functions to validate inputs before algorithm execution.
# Proper validation helps catch errors early and provides informative messages.
#
# Contents:
#   1. Data Structure Validation
#   2. Variance Component Validation

# =============================================================================
# 1. DATA STRUCTURE VALIDATION
# =============================================================================

#' Check a numeric vector or matrix for NA, NaN, or Inf values.
#'
#' @param x numeric vector or matrix
#' @param field_name name of the field (for error message)
#' @param site_label site identifier (for error message)
#' @return invisible(TRUE); stops with informative error on bad values
#' @keywords internal
check_finite <- function(x, field_name, site_label) {
  if (anyNA(x)) {
    n_na <- sum(is.na(x))
    stop(sprintf(
      "%s contains %d NA/NaN value(s) in site '%s'. Remove or impute missing values before calling the algorithm.",
      field_name, n_na, site_label
    ))
  }
  if (any(is.infinite(x))) {
    n_inf <- sum(is.infinite(x))
    stop(sprintf(
      "%s contains %d Inf/-Inf value(s) in site '%s'. All values must be finite.",
      field_name, n_inf, site_label
    ))
  }
  invisible(TRUE)
}

.is_non_default_numeric <- function(value, default, tol = sqrt(.Machine$double.eps)) {
  is.na(value) || !isTRUE(all.equal(as.numeric(value), as.numeric(default),
                                    tolerance = tol, check.attributes = FALSE))
}

.warn_ignored_param <- function(param, value, dgp_type, reason) {
  warning(sprintf(
    "Parameter '%s' = %s is ignored when dgp_type = '%s': %s",
    param, paste(value, collapse = ","), dgp_type, reason
  ), call. = FALSE)
}

#' Validate algorithm input data structure
#'
#' Checks that data_split contains properly formatted data for target and
#' source sites, including required fields, dimensions, and value ranges.
#'
#' @param data_split List of site data, split by site identifier.
#' @param lambda_selection Lambda selection method ("cv" or numeric value).
#'        Optional (default NULL). When NULL, the lambda_selection check is
#'        skipped; this lets comparison estimators that do not use a lambda
#'        (e.g. \code{estimate_tilted_aipw}, \code{estimate_federated_dr},
#'        \code{estimate_pooled_dr}) reuse this validator for data-shape
#'        checks without passing a placeholder.
#' @param family GLM family: "gaussian" or "binomial".
#'        Y values are validated based on the family: binomial requires Y in \{0,1\}.
#' @return TRUE if valid; throws informative error otherwise.
#'
#' @details
#' Required fields in each site's data:
#' \itemize{
#'   \item \code{W_outcome}: Covariate matrix for outcome model (n x p)
#'   \item \code{Z_site}: Covariate matrix for site model (n x p)
#'   \item \code{A}: Treatment indicator vector (binary: 0/1)
#'   \item \code{Y}: Outcome vector (family-specific validation)
#'   \item \code{n}: Sample size
#' }
validate_algorithm_inputs <- function(data_split, lambda_selection = NULL,
                                      family = "binomial", A_val = 1L) {
  # Check for NULL
  if (is.null(data_split)) {
    stop("data_split must not be NULL")
  }
  
  # Check if data_split is a list
  if (!is.list(data_split)) {
    stop("data_split must be a list")
  }
  
  # Check if data_split is empty

  if (length(data_split) == 0) {
    stop("data_split must not be empty")
  }
  
  # Check if target site exists
  if (!"t" %in% names(data_split)) {
    stop("Target site 't' not found in data_split")
  }
  
  # Validate target site data
  target_data <- data_split[["t"]]
  required_fields <- c("W_outcome", "Z_site", "A", "Y", "n")
  
  for (field in required_fields) {
    if (!field %in% names(target_data)) {
      stop(paste("Required field", field, "not found in target data"))
    }
  }
  
  # Check data dimensions
  n_t <- target_data$n
  if (nrow(target_data$W_outcome) != n_t ||
      nrow(target_data$Z_site) != n_t ||
      length(target_data$A) != n_t ||
      length(target_data$Y) != n_t) {
    stop("Inconsistent data dimensions in target site")
  }
  
  # Check for NA/NaN/Inf in target data
  check_finite(target_data$W_outcome, "W_outcome", "t")
  check_finite(target_data$Z_site,    "Z_site",    "t")
  check_finite(target_data$A,         "A",         "t")
  check_finite(target_data$Y,         "Y",         "t")
  
  # Check treatment values
  if (!all(target_data$A %in% c(0, 1))) {
    stop("Treatment variable A must be binary (0/1)")
  }
  
  # Check outcome values based on GLM family
  if (family == "binomial") {
    if (!all(target_data$Y %in% c(0, 1))) {
      stop("Outcome variable Y must be binary (0/1) for family='binomial'. ",
           "Use family='gaussian' for continuous outcomes.")
    }
  }
  # gaussian: no constraint on Y
  
  # Check if there are units with the target treatment value
  if (sum(target_data$A == A_val) == 0) {
    stop(sprintf("No units with A=%d found in target site", A_val))
  }
  
  # Validate source sites
  source_sites <- setdiff(names(data_split), "t")
  if (length(source_sites) == 0) {
    stop("No source sites found in data_split")
  }
  
  for (site in source_sites) {
    source_data <- data_split[[site]]
    
    # Check required fields
    for (field in required_fields) {
      if (!field %in% names(source_data)) {
        stop(paste("Required field", field, "not found in source site", site))
      }
    }
    
    # Check data dimensions
    n_s <- source_data$n
    if (nrow(source_data$W_outcome) != n_s ||
        nrow(source_data$Z_site) != n_s ||
        length(source_data$A) != n_s ||
        length(source_data$Y) != n_s) {
      stop(paste("Inconsistent data dimensions in source site", site))
    }
    
    # Check for NA/NaN/Inf in source data
    check_finite(source_data$W_outcome, "W_outcome", site)
    check_finite(source_data$Z_site,    "Z_site",    site)
    check_finite(source_data$A,         "A",         site)
    check_finite(source_data$Y,         "Y",         site)
    
    # Check treatment and outcome values
    if (!all(source_data$A %in% c(0, 1))) {
      stop(paste("Treatment variable A must be binary (0/1) in source site", site))
    }
    
    if (family == "binomial" && !all(source_data$Y %in% c(0, 1))) {
      stop(paste("Outcome variable Y must be binary (0/1) in source site", site,
                 "for family='binomial'."))
    }
    
    # Check if there are units with the target treatment value
    if (sum(source_data$A == A_val) == 0) {
      stop(sprintf("No units with A=%d found in source site %s", A_val, site))
    }
  }
  
  # Validate lambda_selection only when callers actually use one. Comparison
  # methods that pass lambda_selection = NULL skip this block.
  if (!is.null(lambda_selection)) {
    if (!is.character(lambda_selection) && !is.numeric(lambda_selection)) {
      stop("lambda_selection must be either 'cv' or a numeric value")
    }

    if (is.character(lambda_selection) && lambda_selection != "cv") {
      stop("lambda_selection must be 'cv' if character")
    }

    if (is.numeric(lambda_selection) && (lambda_selection <= 0 || lambda_selection > 1)) {
      stop("lambda_selection must be between 0 and 1 if numeric")
    }
  }

  return(TRUE)
}

# =============================================================================
# 2. CROSS-FITTING SAMPLE SIZE VALIDATION
# =============================================================================

#' Validate sample sizes for cross-fitting algorithms
#'
#' Checks that each site has sufficient sample size for the specified number
#' of cross-fitting folds. With K_f folds, each fold should have at least
#' MIN_FOLD_SIZE observations for reliable estimation.
#'
#' @param data_split List of site data, split by site identifier.
#' @param n_folds Number of cross-fitting folds (K_f).
#' @param min_fold_size Minimum observations per fold (default 5).
#' @return TRUE if valid; throws informative error or warning otherwise.
#'
#' @details
#' For two-level cross-fitting with K_f folds:
#' - Each site needs at least K_f * min_fold_size observations
#' - Target site needs at least K_f * min_fold_size treated units
#' - Source sites need at least K_f * min_fold_size treated units
#' 
#' If sample sizes are insufficient, the function will:
#' - Throw an error if critically insufficient (< 3 per fold)
#' - Issue a warning if marginal (< min_fold_size per fold)
validate_crossfit_sample_sizes <- function(data_split, n_folds, min_fold_size = 5,
                                           A_val = 1L) {
  # Validate n_folds itself
  if (n_folds < 2) {
    stop("n_folds must be at least 2 for cross-fitting.")
  }
  
  if (n_folds == 2) {
    warning(sprintf(
      "validate_crossfit_sample_sizes: n_folds = %d may lead to unstable variance estimation; consider n_folds >= 3 for more reliable inference.",
      n_folds
    ), call. = FALSE)
  }
  
  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  
  # Critical minimum: at least 3 observations per fold for any estimation
  CRITICAL_MIN_PER_FOLD <- 3
  
  # Check target site
  n_t <- target_data$n
  n_t_treated <- sum(target_data$A == A_val)
  
  if (n_t < n_folds * CRITICAL_MIN_PER_FOLD) {
    stop(sprintf(
      "Target site has insufficient sample size for %d-fold cross-fitting. Need at least %d observations, but only have %d.",
      n_folds, n_folds * CRITICAL_MIN_PER_FOLD, n_t
    ))
  }
  
  if (n_t_treated < n_folds * CRITICAL_MIN_PER_FOLD) {
    stop(sprintf(
      "Target site has insufficient A=%d units for %d-fold cross-fitting. Need at least %d, but only have %d.",
      A_val, n_folds, n_folds * CRITICAL_MIN_PER_FOLD, n_t_treated
    ))
  }
  
  # Warning if marginal
  if (n_t < n_folds * min_fold_size) {
    warning(sprintf(
      "Target site sample size (%d) may be insufficient for stable %d-fold cross-fitting. Recommend at least %d observations.",
      n_t, n_folds, n_folds * min_fold_size
    ))
  }
  
  # Check source sites
  for (site in source_sites) {
    source_data <- data_split[[site]]
    n_s <- source_data$n
    n_s_treated <- sum(source_data$A == A_val)
    
    if (n_s < n_folds * CRITICAL_MIN_PER_FOLD) {
      stop(sprintf(
        "Source site '%s' has insufficient sample size for %d-fold cross-fitting. Need at least %d observations, but only have %d.",
        site, n_folds, n_folds * CRITICAL_MIN_PER_FOLD, n_s
      ))
    }
    
    if (n_s_treated < n_folds * CRITICAL_MIN_PER_FOLD) {
      stop(sprintf(
        "Source site '%s' has insufficient A=%d units for %d-fold cross-fitting. Need at least %d, but only have %d.",
        site, A_val, n_folds, n_folds * CRITICAL_MIN_PER_FOLD, n_s_treated
      ))
    }
    
    # Warning if marginal
    if (n_s < n_folds * min_fold_size) {
      warning(sprintf(
        "Source site '%s' sample size (%d) may be insufficient for stable %d-fold cross-fitting.",
        site, n_s, n_folds
      ))
    }
  }
  
  return(TRUE)
}

#' Validate truncation parameters used in FACE-C
#'
#' @param M_tau Truncation bound for calibrated nuisance optimization.
#' @param M_tau_inference Truncation bound for inference-stage DR terms.
#' @return TRUE invisibly; throws error on invalid inputs.
validate_truncation_parameters <- function(M_tau, M_tau_inference) {
  if (!is.numeric(M_tau) || length(M_tau) != 1L || is.na(M_tau) || M_tau <= 0) {
    stop("M_tau must be a single positive numeric value.")
  }

  valid_inference <- is.numeric(M_tau_inference) &&
    length(M_tau_inference) == 1L &&
    !is.na(M_tau_inference) &&
    (is.infinite(M_tau_inference) || M_tau_inference > 0)

  if (!valid_inference) {
    stop("M_tau_inference must be a single positive numeric value or Inf.")
  }

  invisible(TRUE)
}

# =============================================================================
# 3. VARIANCE COMPONENT VALIDATION
# =============================================================================

# =============================================================================
# 4. SIMULATION PARAMETER VALIDATION
# =============================================================================

#' Validate simulation configuration parameters (single source of truth)
#'
#' Centralized validation for all simulation parameters. Called once at the
#' entry point (main.R) and by generate_simulation_data(). Replaces the
#' duplicated validation blocks that previously existed in multiple files.
#'
#' @param estimand_type "sample" or "superpopulation"
#' @param site_allocation "model", "uniform", "balanced", "target_heavy", "source_heavy", or "independent"
#' @param transform_type "strong", "mild", or "none"
#' @param outcome_type "binary" or "continuous"
#' @param heterogeneity_type "none", "mild", "strong", or "partial"
#' @param shift_strength Positive numeric multiplier
#' @param n_folds Integer >= 3
#' @param n_sims Positive integer (NULL to skip this check)
#' @param n_total Positive integer total sample size (NULL to skip this check)
#' @param K Positive integer number of source sites (NULL to skip this check)
#' @param p Positive integer number of covariates (NULL to skip this check)
#' @param config Configuration string (NULL to skip, else must be "C1"-"C4")
#' @param dgp_type "facec" (default) or "face"
#' @param ate_deviation Numeric ATE deviation for FACE paper non-informative sites (>= 0)
#' @param n_deviated_sites Non-negative integer: how many source sites deviate (FACE paper only)
#' @param warn_ignored Logical. If TRUE, warn when parameters are accepted for
#'   API compatibility but ignored by the selected DGP.
#' @return TRUE invisibly; throws informative error on invalid input
#' @export
validate_simulation_params <- function(estimand_type = "superpopulation",
                                        site_allocation = "model",
                                        transform_type = "mild",
                                        outcome_type = "binary",
                                        heterogeneity_type = "none",
                                        shift_strength = 1.0,
                                        n_folds = N_FOLDS_DEFAULT,
                                        n_sims = NULL,
                                        n_total = NULL,
                                        K = NULL,
                                        p = NULL,
                                        config = NULL,
                                        dgp_type = "facec",
                                        ate_deviation = 0.0,
                                        n_deviated_sites = 0L,
                                        warn_ignored = TRUE) {
  if (!(dgp_type %in% VALID_DGP_TYPES)) {
    stop(sprintf("Invalid dgp_type: '%s'. Must be one of: %s",
                 dgp_type, paste(VALID_DGP_TYPES, collapse = ", ")))
  }

  if (!(estimand_type %in% VALID_ESTIMAND_TYPES)) {
    stop(sprintf("Invalid estimand_type: '%s'. Must be one of: %s",
                 estimand_type, paste(VALID_ESTIMAND_TYPES, collapse = ", ")))
  }

  # site_allocation, transform_type, outcome_type, heterogeneity_type, and
  # shift_strength are specific to the FACE-C DGP.  For the FACE paper DGP
  # these parameters are not used, so we skip their validation (they keep
  # their defaults and generate a warning when non-default values are passed).
  if (dgp_type == "facec") {
    if (!(site_allocation %in% VALID_SITE_ALLOCATIONS)) {
      stop(sprintf("Invalid site_allocation: '%s'. Must be one of: %s",
                   site_allocation, paste(VALID_SITE_ALLOCATIONS, collapse = ", ")))
    }

    if (!(transform_type %in% VALID_TRANSFORM_TYPES)) {
      stop(sprintf("Invalid transform_type: '%s'. Must be one of: %s",
                   transform_type, paste(VALID_TRANSFORM_TYPES, collapse = ", ")))
    }

    if (!(outcome_type %in% VALID_OUTCOME_TYPES)) {
      stop(sprintf("Invalid outcome_type: '%s'. Must be one of: %s",
                   outcome_type, paste(VALID_OUTCOME_TYPES, collapse = ", ")))
    }

    if (!(heterogeneity_type %in% VALID_HETEROGENEITY_TYPES)) {
      stop(sprintf("Invalid heterogeneity_type: '%s'. Must be one of: %s",
                   heterogeneity_type, paste(VALID_HETEROGENEITY_TYPES, collapse = ", ")))
    }

    if (is.na(shift_strength) || shift_strength <= 0) {
      stop("shift_strength must be a positive number.")
    }

    if (isTRUE(warn_ignored)) {
      if (.is_non_default_numeric(ate_deviation, 0.0)) {
        .warn_ignored_param("ate_deviation", ate_deviation, dgp_type,
                            "only the FACE paper DGP uses source-site ATE deviations.")
      }
      if (.is_non_default_numeric(n_deviated_sites, 0L)) {
        .warn_ignored_param("n_deviated_sites", n_deviated_sites, dgp_type,
                            "only the FACE paper DGP uses deviated source-site counts.")
      }
    }

  } else if (dgp_type == "face") {
    # FACE paper DGP-specific validation
    if (is.na(ate_deviation) || ate_deviation < 0) {
      stop("ate_deviation must be a non-negative number.")
    }
    if (is.na(n_deviated_sites) || n_deviated_sites < 0) {
      stop("n_deviated_sites must be a non-negative integer.")
    }
    if (!(outcome_type %in% VALID_OUTCOME_TYPES)) {
      stop(sprintf("Invalid outcome_type: '%s'. Must be one of: %s",
                   outcome_type, paste(VALID_OUTCOME_TYPES, collapse = ", ")))
    }

    if (isTRUE(warn_ignored)) {
      if (!identical(site_allocation, "model")) {
        .warn_ignored_param("site_allocation", site_allocation, dgp_type,
                            "FACE paper DGP has its own site-specific covariate generator.")
      }
      if (!identical(transform_type, "mild")) {
        .warn_ignored_param("transform_type", transform_type, dgp_type,
                            "FACE paper DGP builds its own linear/quadratic design matrices.")
      }
      if (!identical(heterogeneity_type, "none")) {
        .warn_ignored_param("heterogeneity_type", heterogeneity_type, dgp_type,
                            "FACE paper DGP uses ate_deviation/n_deviated_sites instead.")
      }
      if (.is_non_default_numeric(shift_strength, 1.0)) {
        .warn_ignored_param("shift_strength", shift_strength, dgp_type,
                            "FACE paper DGP uses fixed skewed-normal covariate shifts.")
      }
    }
  }

  if (!is.null(n_total)) {
    if (any(is.na(n_total)) || any(n_total < 1) ||
        any(abs(n_total - round(n_total)) > sqrt(.Machine$double.eps))) {
      stop("n_total must contain only positive integers.")
    }
  }

  if (!is.null(K)) {
    if (any(is.na(K)) || any(K < 1) ||
        any(abs(K - round(K)) > sqrt(.Machine$double.eps))) {
      stop("K must contain only positive integer source-site counts.")
    }
  }

  if (!is.null(p)) {
    if (any(is.na(p)) || any(p < 1) ||
        any(abs(p - round(p)) > sqrt(.Machine$double.eps))) {
      stop("p must contain only positive integer covariate counts.")
    }
  }

  if (!is.null(n_total) && !is.null(K) && min(n_total) < max(K) + 1L) {
    stop(sprintf(
      "n_total contains value(s) too small for K: min(n_total)=%d, max(K)=%d, so at least max(K) + 1 observations are required.",
      as.integer(min(n_total)), as.integer(max(K))
    ))
  }

  if (is.na(n_folds) || n_folds < 3 ||
      abs(n_folds - round(n_folds)) > sqrt(.Machine$double.eps)) {
    stop("n_folds must be an integer >= 3.")
  }

  if (!is.null(n_sims) && (is.na(n_sims) || n_sims < 1 ||
                           abs(n_sims - round(n_sims)) > sqrt(.Machine$double.eps))) {
    stop("n_sims must be a positive integer.")
  }

  if (!is.null(config)) {
    if (any(!(config %in% VALID_CONFIGS))) {
      stop(sprintf("Invalid config: '%s'. Must be one of: %s",
                   paste(config[!(config %in% VALID_CONFIGS)], collapse = ","),
                   paste(VALID_CONFIGS, collapse = ", ")))
    }
  }

  invisible(TRUE)
}
