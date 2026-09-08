# data_generation.R - Simulation data generation for RoCE algorithms
#
# This file contains functions to generate synthetic data for Monte Carlo
# simulations following the data generating process described in main.tex.
#
# Contents:
#   1. Covariate Generation
#   2. Site/Treatment Assignment
#   3. Outcome Generation (RoCE DGP)
#   4. Full Simulation Data Generation (RoCE DGP)
#   4b. FACE Paper DGP → extracted to R/data_generation_face.R
#   5. Data Splitting Utilities

#' Generate base covariates
#' @param n number of observations
#' @param p number of covariates (default 4)
#' @return matrix of base covariates
#' @export
generate_base_covariates <- function(n, p = 4) {
  # Generate from multivariate normal
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  
  # Clip to [-2.5, 2.5]
  X <- apply(X, 2, clip_to_range)
  
  # Standardize
  X <- scale_center(X)
  
  colnames(X) <- paste0("X", 1:p)
  return(X)
}

#' Apply nonlinear transformations to covariates
#' @param X matrix of base covariates
#' @param transform_type type of transformation:
#'   - "strong": original aggressive transformations (exp, cube, square)
#'   - "mild": gentle transformations (square, interaction) that preserve IF orthogonality
#'   - "none": no transformation (X_dagger = X)
#' @return matrix of transformed covariates
#' @details
#' The "strong" transformation uses highly nonlinear functions that can break
#' the Neyman orthogonality condition required for valid IF-based variance estimation.
#' The "mild" transformation uses only quadratic terms and simple interactions,
#' which are easier to approximate with linear models and preserve orthogonality.
#' @export
transform_covariates <- function(X, transform_type = "strong") {
  n <- nrow(X)
  p <- ncol(X)
  
  # Validate transform_type
  if (!(transform_type %in% c("strong", "mild", "none"))) {
    stop(sprintf("Invalid transform_type: '%s'. Must be one of: 'strong', 'mild', 'none'.", transform_type))
  }
  
  # Initialize transformed matrix
  X_dagger <- matrix(0, nrow = n, ncol = p)
  
  if (transform_type == "none") {
    # No transformation: X_dagger = X
    X_dagger <- X
    
  } else if (transform_type == "mild") {
    # Mild transformations: preserve IF orthogonality
    # Use only quadratic terms and simple interactions
    X_dagger[, 1] <- scale_center(X[, 1]^2)                    # Simple square
    X_dagger[, 2] <- scale_center(X[, 1] * X[, 2])             # Interaction
    X_dagger[, 3] <- scale_center(X[, 3]^2)                    # Simple square
    X_dagger[, 4] <- scale_center(X[, 2] + X[, 4])             # Linear combination
    
    # If more than 4 covariates, use mild transformations
    if (p > 4) {
      for (j in 5:p) {
        X_dagger[, j] <- scale_center(X[, j]^2)  # Simple square
      }
    }
    
  } else {
    # Strong transformations (original): can break IF orthogonality
    # Apply transformations as specified in the document
    X_dagger[, 1] <- scale_center(exp(0.5 * X[, 1]))
    X_dagger[, 2] <- scale_center(10 + X[, 2] / (1 + exp(X[, 1])))
    X_dagger[, 3] <- scale_center((0.04 * X[, 1] * X[, 3] + 0.5)^3)
    X_dagger[, 4] <- scale_center((X[, 2] + X[, 4] + 20)^2)
    
    # If more than 4 covariates, keep the rest unchanged
    if (p > 4) {
      X_dagger[, 5:p] <- X[, 5:p]
    }
  }
  
  colnames(X_dagger) <- paste0("X", 1:p, "_dagger")
  return(X_dagger)
}

#' Generate site-specific γ (site/treatment model) parameters
#'
#' Following main.tex: Multinomial logistic model parameters with unit norm
#' constraint and a balancing intercept.
#'
#' \strong{Parameter structure:}
#' \deqn{\boldsymbol{\gamma}_{s_j,a}
#'   = (c_0, \epsilon_{j,a,1}, \epsilon_{j,a,2}, \mathbf{0}_{p-2}).}
#' Here \eqn{\epsilon_{j,a,k} \sim \mathcal{N}(0, \sigma_\epsilon^2)}, with
#' \eqn{\sigma_\epsilon = \mathtt{sigma\_epsilon} \times
#' \mathtt{shift\_strength}}, and
#' \eqn{c_0 = \mathtt{GAMMA\_BALANCE\_INTERCEPT}} is a negative intercept chosen so
#' that \eqn{P(R=t) \approx 1/(K+1)} under the unit-strength setting.
#'
#' \strong{Balancing derivation:} In the multinomial logistic model the target is
#' the baseline category (all-zero γ).  Each of the 2K source categories has
#' \eqn{\mathbb{E}[\exp(\gamma^T Z)]
#' = \exp(c_{\rm eff} + (1 - c_{\rm eff}^2)/2)} under unit norm. Setting this
#' equal to \eqn{1/2} yields
#' \eqn{P(R=t) = 1/(1 + 2K \cdot 1/2) = 1/(K+1)} regardless of K.
#'
#' \strong{Unit norm constraint:}
#' \deqn{\|\boldsymbol{\gamma}_{s_j,a}\|_2 = 1 \quad \forall j, a.}
#' The full vector (including intercept) is normalized to unit L2 norm.
#'
#' @param K number of source sites
#' @param p number of covariates (excluding intercept)
#' @param sigma_epsilon base standard deviation for variations (default 0.1)
#' @param shift_strength multiplier for sigma_epsilon controlling covariate shift
#'        intensity (default \code{ROCE_SHIFT_STRENGTH_DEFAULT})
#' @return list with γ parameters for each site and treatment combination
generate_site_model_parameters <- function(K, p, sigma_epsilon = 0.1,
                                           shift_strength = ROCE_SHIFT_STRENGTH_DEFAULT) {
  gamma_list <- list()

  for (j in 1:K) {
    # Base parameters with zeros, length p+1 (including intercept)
    base_gamma_0 <- rep(0, p + 1)
    base_gamma_1 <- rep(0, p + 1)

    # γ_{s_j,a} = (c0, ε_{j,a,1}, ε_{j,a,2}, 0_{p-2})
    # Negative intercept c0 balances target vs source site sizes:
    # P(target) ≈ 1/(K+1) regardless of K.
    # IMPORTANT: Do NOT scale the intercept by shift_strength.
    # We normalize the full vector to unit norm below; scaling both intercept
    # and covariate components would cancel out and make shift_strength a no-op.
    gamma_intercept <- GAMMA_BALANCE_INTERCEPT
    base_gamma_0[1] <- gamma_intercept
    base_gamma_1[1] <- gamma_intercept

    # Add small variations to first 2 covariate parameters
    if (p >= 2) {
      # shift_strength scales the variation: larger values produce stronger covariate shift.
      epsilon_0 <- rnorm(2, 0, sigma_epsilon * shift_strength)
      epsilon_1 <- rnorm(2, 0, sigma_epsilon * shift_strength)

      base_gamma_0[2:3] <- epsilon_0
      base_gamma_1[2:3] <- epsilon_1
      # Components 4:(p+1) remain 0 as specified
    } else if (p > 0) {
      # For p < 2, add variations to all available covariate parameters
      epsilon_0 <- rnorm(p, 0, sigma_epsilon * shift_strength)
      epsilon_1 <- rnorm(p, 0, sigma_epsilon * shift_strength)

      base_gamma_0[2:(p+1)] <- epsilon_0
      base_gamma_1[2:(p+1)] <- epsilon_1
    }

    # Normalize FULL vector to unit norm (including intercept)
    # This matches main.tex lines 284-287: ||γ_{s_j,a}||_2 = 1
    full_norm_0 <- sqrt(sum(base_gamma_0^2))
    full_norm_1 <- sqrt(sum(base_gamma_1^2))

    if (full_norm_0 > 0) {
      base_gamma_0 <- base_gamma_0 / full_norm_0
    }
    if (full_norm_1 > 0) {
      base_gamma_1 <- base_gamma_1 / full_norm_1
    }

    # Store the parameters
    gamma_list[[paste0("s", j, "_0")]] <- base_gamma_0
    gamma_list[[paste0("s", j, "_1")]] <- base_gamma_1
  }

  return(gamma_list)
}

#' Calculate site membership probabilities
#' Following eq:site_membership in main.tex: P(R=t|X) and P(R=s_j,A=a|X) using multinomial logistic model
#' @param Z covariate matrix (X or X_dagger)
#' @param gamma_params list of γ (site/treatment model) parameters
#' @param K number of source sites
#' @return list with probabilities for each site and treatment combination
calculate_site_probabilities <- function(Z, gamma_params, K) {
  n <- nrow(Z)
  
  # Add intercept column to Z matrix
  Z_with_intercept <- cbind(1, Z)
  
  # Calculate linear predictors for all sites and treatments
  g_values <- matrix(0, nrow = n, ncol = 2 * K)
  
  for (j in 1:K) {
    g_values[, 2*j - 1] <- Z_with_intercept %*% gamma_params[[paste0("s", j, "_0")]]
    g_values[, 2*j] <- Z_with_intercept %*% gamma_params[[paste0("s", j, "_1")]]
  }
  
  # Calculate denominator
  exp_g <- exp(g_values)
  denominator <- 1 + rowSums(exp_g)
  
  # Calculate probabilities
  probs <- list()
  
  # Target site probability
  probs$p_target <- 1 / denominator
  
  # Source site probabilities
  for (j in 1:K) {
    probs[[paste0("p_s", j, "_0")]] <- exp_g[, 2*j - 1] / denominator
    probs[[paste0("p_s", j, "_1")]] <- exp_g[, 2*j] / denominator
  }
  
  return(probs)
}

#' Compute marginal treatment propensity from multinomial site probabilities
#'
#' Derives P(A=1|X) from the multinomial logistic model probabilities.
#' This is used consistently across all site allocation methods to maintain
#' the same confounding structure regardless of how subjects are assigned to sites.
#'
#' Formula: P(A=1|X) = Σ_j P(R=s_j, A=1|X) / Σ_j P(R=s_j|X)
#'
#' @param site_probs Site probability list from calculate_site_probabilities()
#' @param K Number of source sites
#' @param n Number of observations
#' @return Numeric vector of treatment propensities, length n,
#'   clipped to [POSITIVITY_LOWER, POSITIVITY_UPPER]
calculate_treatment_propensity <- function(site_probs, K, n) {
  n_cat <- 1 + 2 * K
  pm <- matrix(0, nrow = n, ncol = n_cat)
  pm[, 1] <- site_probs$p_target
  for (j in 1:K) {
    pm[, 2 * j]     <- site_probs[[paste0("p_s", j, "_0")]]
    pm[, 2 * j + 1] <- site_probs[[paste0("p_s", j, "_1")]]
  }
  p_num <- rowSums(pm[, seq(3, n_cat, by = 2), drop = FALSE])
  p_den <- rowSums(pm[, 2:n_cat, drop = FALSE])
  p_treat <- pmax(pmin(p_num / pmax(p_den, DIVISION_FLOOR),
                       POSITIVITY_UPPER), POSITIVITY_LOWER)
  return(p_treat)
}

#' Generate site membership and treatment assignments (vectorized)
#' @param n number of observations
#' @param site_probs site probability list
#' @param K number of source sites
#' @return list with R (site membership) and A (treatment) vectors
generate_site_treatment_assignments <- function(n, site_probs, K) {
  # Build probability matrix: each row is an observation, each column is a (site, treatment) combination
  # Column order: target, s1_0, s1_1, s2_0, s2_1, ...
  n_categories <- 1 + 2 * K  # target + (control + treated) for each source
  prob_matrix <- matrix(0, nrow = n, ncol = n_categories)
  
  # First column: target site probability
  prob_matrix[, 1] <- site_probs$p_target
  
  # Source site probabilities
    for (j in 1:K) {
    prob_matrix[, 2 * j] <- site_probs[[paste0("p_s", j, "_0")]]      # control
    prob_matrix[, 2 * j + 1] <- site_probs[[paste0("p_s", j, "_1")]]  # treated
  }
  
  # Vectorized multinomial sampling using cumulative probabilities
  # For each row, sample one category based on probabilities
  cumprob <- t(apply(prob_matrix, 1, cumsum))
  u <- runif(n)
  outcomes <- rowSums(u > cumprob) + 1  # Which category was selected
  
  # Initialize output vectors
  R <- character(n)
  A <- integer(n)
  
  # Target site (outcome == 1)
  target_idx <- outcomes == 1
  R[target_idx] <- "t"
  
  # For target site, calculate treatment probability based on source site proportions
  if (any(target_idx)) {
    # Use the same propensity formula as non-model allocations
    # Only compute for target-site observations
    target_probs <- list()
    target_probs$p_target <- site_probs$p_target[target_idx]
    for (j_inner in 1:K) {
      target_probs[[paste0("p_s", j_inner, "_0")]] <- site_probs[[paste0("p_s", j_inner, "_0")]][target_idx]
      target_probs[[paste0("p_s", j_inner, "_1")]] <- site_probs[[paste0("p_s", j_inner, "_1")]][target_idx]
    }
    p_treat <- calculate_treatment_propensity(target_probs, K, sum(target_idx))
    A[target_idx] <- rbinom(sum(target_idx), 1, p_treat)
  }
  
  # Source sites (outcome > 1)
  source_idx <- outcomes > 1
  if (any(source_idx)) {
    source_outcomes <- outcomes[source_idx]
    # Site index: outcome 2,3 -> s1; 4,5 -> s2; etc.
    site_nums <- ceiling((source_outcomes - 1) / 2)
    R[source_idx] <- paste0("s", site_nums)
    # Treatment: even outcomes (2,4,6,...) -> 0; odd outcomes (3,5,7,...) -> 1
    A[source_idx] <- source_outcomes %% 2
  }
  
  return(list(R = R, A = as.integer(A)))
}

#' Generate outcomes using GLM framework (vectorized)
#' Following eq:linear_phi_assumption in main.tex: GLM outcome models
#' E[Y(a)|W,R=r] = h(W; α_r^a) where h depends on outcome_type
#' @param W covariate matrix for outcome model
#' @param A treatment vector
#' @param R site membership vector
#' @param alpha1 treatment outcome parameters
#' @param alpha0 control outcome parameters
#' @param noise_sd standard deviation of noise (used for continuous outcomes)
#' @param outcome_type type of outcome: "binary" (logistic link) or "continuous"
#'   (identity link with Gaussian noise)
#' @param heterogeneity_type type of outcome heterogeneity across sites:
#'   - "none" (default): all sites use same parameters
#'   - "mild": approximately 40\% change in 2 non-zero coefficients for all source sites
#'   - "strong": approximately 80\% change in 2 non-zero coefficients for all source sites
#'   - "partial": only first half of source sites have mild heterogeneity
#' @return observed outcomes (binary: 0/1 for binary; continuous for continuous)
generate_outcomes <- function(W, A, R, alpha1, alpha0, noise_sd = ROCE_NOISE_SD_DEFAULT,
                              outcome_type = "binary",
                              heterogeneity_type = "none") {
  n <- length(A)
  
  # Create design matrix with intercept
  W_design <- cbind(1, W)
  
  # Compute source-site-specific parameters based on heterogeneity_type
  if (heterogeneity_type != "none") {
    # Determine which source sites are heterogeneous
    unique_sources <- sort(unique(R[R != "t"]))
    K_actual <- length(unique_sources)
    
    if (heterogeneity_type == "partial") {
      # Only first half of source sites get heterogeneous parameters
      n_hetero <- max(1, floor(K_actual / 2))
      hetero_sites <- unique_sources[1:n_hetero]
    } else {
      # All source sites are heterogeneous
      hetero_sites <- unique_sources
    }
    
    # Create source-specific parameter modifications
    alpha1_source <- alpha1
    alpha0_source <- alpha0
    
    if (heterogeneity_type %in% c("mild", "partial")) {
      # Mild heterogeneity: ~40% change in 2 non-zero coefficients
      # Original: alpha1[2] = 0.5, alpha1[3] = 0.75
      # Changed to: alpha1_source[2] = 0.3 (-40%), alpha1_source[3] = 0.45 (-40%)
      alpha1_source[2] <- 0.3
      alpha1_source[3] <- 0.45
    } else if (heterogeneity_type == "strong") {
      # Strong heterogeneity: ~80% change in 2 non-zero coefficients
      # Makes source sites substantially different from target
      alpha1_source[2] <- 0.1   # was 0.5 (-80%)
      alpha1_source[3] <- 0.15  # was 0.75 (-80%)
      # Also modify control outcome for asymmetric bias
      alpha0_source[2] <- 1.5   # was 0.75 (+100%)
    }
    
    # Identify target and source observations
    is_target <- R == "t"
    is_hetero_source <- R %in% hetero_sites
    is_homo_source <- (!is_target) & (!is_hetero_source)
    
    # Compute linear predictors vectorized
    eta1 <- numeric(n)
    eta0 <- numeric(n)
    
    # Target site: use original parameters
    if (any(is_target)) {
      eta1[is_target] <- as.numeric(W_design[is_target, , drop = FALSE] %*% alpha1)
      eta0[is_target] <- as.numeric(W_design[is_target, , drop = FALSE] %*% alpha0)
    }
    
    # Heterogeneous source sites: use modified parameters
    if (any(is_hetero_source)) {
      eta1[is_hetero_source] <- as.numeric(W_design[is_hetero_source, , drop = FALSE] %*% alpha1_source)
      eta0[is_hetero_source] <- as.numeric(W_design[is_hetero_source, , drop = FALSE] %*% alpha0_source)
    }
    
    # Homogeneous source sites (for partial): use original parameters
    if (any(is_homo_source)) {
      eta1[is_homo_source] <- as.numeric(W_design[is_homo_source, , drop = FALSE] %*% alpha1)
      eta0[is_homo_source] <- as.numeric(W_design[is_homo_source, , drop = FALSE] %*% alpha0)
    }
  } else {
    # Homogeneous case: all sites use same parameters (fully vectorized)
    eta1 <- as.numeric(W_design %*% alpha1)
    eta0 <- as.numeric(W_design %*% alpha0)
  }
  
  if (outcome_type == "binary") {
    # Apply logistic link function (vectorized)
    mu1 <- logistic(eta1)
    mu0 <- logistic(eta0)
    
    # Generate outcomes based on treatment assignment (vectorized)
    mu <- ifelse(A == 1, mu1, mu0)
    Y <- rbinom(n, 1, mu)
  } else if (outcome_type == "continuous") {
    # Identity link with Gaussian noise
    mu1 <- eta1
    mu0 <- eta0
    
    mu <- ifelse(A == 1, mu1, mu0)
    Y <- mu + rnorm(n, 0, noise_sd)
  } else {
    stop(paste0("Invalid outcome_type: '", outcome_type, "'. Must be 'binary' or 'continuous'."))
  }
  
  return(Y)
}

#' Compute superpopulation true potential outcome mean (TARGET POPULATION SPECIFIC)
#' 
#' This function computes a fixed "superpopulation" true value E_t[Y(1)] that
#' is specific to the TARGET population. It accounts for covariate shift between
#' the overall population and the target site using importance weighting.
#'
#' @param p number of base covariates
#' @param K number of source sites (needed for site probability calculation)
#' @param config configuration ("C1", "C2", "C3", or "C4"). The RoCE
#'   superpopulation truth is fixed across configurations; \code{config}
#'   is validated here for interface consistency with \code{generate_simulation_data()}.
#' @param n_ref size of reference population for Monte Carlo integration (default 100000)
#' @param ref_seed fixed seed for reference population (default 99999)
#' @param transform_type type of covariate transformation ("strong", "mild", or "none")
#' @param site_allocation Site-allocation mechanism. For \code{"model"}, the
#'   target truth is weighted by the covariate-dependent target probability;
#'   otherwise allocation is independent of the target covariate law.
#' @param outcome_type Outcome family, \code{"binary"} or \code{"continuous"}.
#' @param shift_strength Positive multiplier controlling covariate shift in the
#'   site-allocation model.
#' @return list with mu1_superpop (E_t[Y(1)]) and mu0_superpop (E_t[Y(0)])
#' @details
#' The TARGET-SPECIFIC superpopulation parameter is defined as:
#' \deqn{E_t[Y(1)] = E[Y(1)\mid R=t]
#' = \frac{E_X[\operatorname{logistic}((X^\dagger)^T\alpha_1)
#' P(R=t\mid X^\dagger)]}{E_X[P(R=t\mid X^\dagger)]}.}
#' 
#' This is computed via importance-weighted Monte Carlo integration:
#' 1. Generate large reference population with fixed seed
#' 2. Generate fixed alpha parameters (for consistent site probabilities)
#' 3. Compute P(R=t|X_dagger) for each observation
#' 4. Weight potential outcomes by P(R=t|X_dagger) to get target-specific expectation
#'
#' The same seeds ensure all simulations use the same true value.
#' 
#' \strong{Why target-specific?}
#' RoCE's primary estimand is the TATE μ^1_t - μ^0_t in the TARGET
#' population, with each arm also available as a secondary estimand. Due to
#' covariate shift, these target-arm means differ from their overall-population
#' counterparts. Target-specific superpopulation truths therefore provide the
#' correct reference for both arms and their contrast.
calculate_superpopulation_truth <- function(p, K = 3, config,
                                            n_ref = 100000, ref_seed = 99999,
                                            transform_type = "mild",
                                            site_allocation = "model",
                                            outcome_type = "binary",
                                            shift_strength = ROCE_SHIFT_STRENGTH_DEFAULT) {
  if (!(config %in% VALID_CONFIGS)) {
    stop(sprintf("Invalid config: '%s'. Must be one of: %s",
                 config, paste(VALID_CONFIGS, collapse = ", ")))
  }

  with_seed(ref_seed, {
    # IMPORTANT: generate gamma_params FIRST before any other RNG calls.
    # This ensures the same RNG state as in generate_simulation_data(),
    # which also calls with_seed(99999, generate_site_model_parameters(...)).
    gamma_params <- generate_site_model_parameters(K, p, shift_strength = shift_strength)

    # Generate large reference population (RNG state continues after gamma generation)
    #
    # For independent allocation, R is deterministic and the estimand is
    # E[Y(1) | R=t] under the TARGET site's covariate distribution.
    # To match generate_simulation_data(..., site_allocation="independent"),
    # we therefore integrate over the target covariate distribution used there.
    if (site_allocation == "independent") {
      X_ref <- matrix(rnorm(n_ref * p), nrow = n_ref, ncol = p)
      X_ref <- apply(X_ref, 2, clip_to_range)
      colnames(X_ref) <- paste0("X", 1:p)
    } else {
      X_ref <- generate_base_covariates(n_ref, p)
    }
    X_dagger_ref <- transform_covariates(X_ref, transform_type = transform_type)

    # The RoCE DGP is fixed across C1-C4. Configurations only change the
    # fitted working bases exposed to estimators in generate_simulation_data().
    Z_site_ref <- X_dagger_ref
    W_ref <- X_dagger_ref

    # Calculate site probabilities P(R=t|X_dagger), P(R=s_j,A=a|X_dagger)
    site_probs <- calculate_site_probabilities(Z_site_ref, gamma_params, K)

    # Target probability weights depend on allocation method.
    # For model allocation, P(R=t|X_dagger) varies with X (covariate shift).
    # For non-model allocations where R is independent of X, E_t[.] = E[.].
    # For independent allocation, we already generated X_ref from the target
    # distribution, so uniform weights are appropriate.
    if (site_allocation == "model") {
      p_target <- site_probs$p_target
    } else {
      p_target <- rep(1, n_ref)
    }

    # True outcome parameters (single source of truth)
    alphas <- get_true_outcome_parameters(p)
    alpha1 <- alphas$alpha1
    alpha0 <- alphas$alpha0

    # Design matrix with intercept
    W_design <- cbind(1, W_ref)

    # Compute conditional means E[Y(a)|X_dagger] = logistic(X_dagger^T alpha_a)
    eta1 <- as.numeric(W_design %*% alpha1)
    eta0 <- as.numeric(W_design %*% alpha0)

    mu1_given_x <- logistic(eta1)
    mu0_given_x <- logistic(eta0)

    # For continuous outcomes, truth is identity link (no logistic)
    if (outcome_type == "continuous") {
      mu1_given_x <- eta1
      mu0_given_x <- eta0
    }

    # Target-specific superpopulation parameters via importance weighting.
    normalizing_constant <- mean(p_target)
    mu1_superpop <- sum(mu1_given_x * p_target) / sum(p_target)
    mu0_superpop <- sum(mu0_given_x * p_target) / sum(p_target)

    # Also compute overall (marginal) superpopulation for reference.
    mu1_overall <- mean(mu1_given_x)
    mu0_overall <- mean(mu0_given_x)

    return(list(
      # Target-specific superpopulation (what RoCE estimates)
      mu1_superpop = mu1_superpop,
      mu0_superpop = mu0_superpop,
      ate_superpop = mu1_superpop - mu0_superpop,
      # Overall (marginal) superpopulation for reference
      mu1_overall = mu1_overall,
      mu0_overall = mu0_overall,
      ate_overall = mu1_overall - mu0_overall,
      # Diagnostics
      mean_p_target = normalizing_constant,
      n_ref = n_ref,
      ref_seed = ref_seed
    ))
  })  # end with_seed
}

# =============================================================================
# 4b. FACE Paper DGP (Han et al., JASA 2023, Section 5.1)
# =============================================================================
# Extracted to R/data_generation_face.R for better organization.
# Functions: generate_skewed_normal, generate_face_covariates,
#   calculate_face_propensity, build_face_ate_map, generate_face_outcomes,
#   calculate_face_truth, generate_face_data
# =============================================================================

#' Generate complete dataset for simulation
#' @param n_total total number of observations
#' @param K number of source sites
#' @param p number of base covariates
#' @param config configuration ("C1", "C2", "C3", or "C4"). Under
#'   \code{dgp_type = "roce"}, C1 uses the exact transformed basis for both
#'   nuisances, C2 reduces only the outcome basis, C3 reduces only the site
#'   basis, and C4 reduces both. Under \code{dgp_type = "face"}, C1--C2 use a
#'   rich quadratic calibration basis and C1--C3 use the correct quadratic
#'   outcome basis; see \code{generate_face_data()} for the important
#'   skew-normal density-ratio qualification.
#' @param estimand_type type of estimand to use for true value calculation:
#'   - "sample": sample-specific true value E_n[Y(1)] based on realized
#'     covariates. This varies across simulations and is appropriate for 
#'     conditional inference. SE should be compared with SD(bias).
#'   - "superpopulation" (default): fixed superpopulation parameter
#'     E[Y(1)] = E_X[logistic(X_dagger^T α)] computed via Monte Carlo integration.
#'     This is the same across all simulations and appropriate for
#'     unconditional/marginal inference. SE should be compared with SD(estimate).
#' @param site_allocation method for site allocation:
#'   - "model" (default): multinomial logistic model (joint site & treatment)
#'   - "uniform": equal probability per site, propensity-based treatment
#'   - "balanced": deterministic equal site sizes, propensity-based treatment
#'   - "target_heavy"/"source_heavy": skewed allocation proportions
#'   - "independent": each site generates its own covariates from its own
#'     distribution. Target: X ~ N(0, I_p). Source j: X ~ N(mu_j, I_p) with
#'     deterministic shift mu_j controlled by shift_strength.
#'   Non-model and non-independent methods use the same covariate-dependent
#'   propensity P(A=1|X) as the model allocation to maintain consistent
#'   confounding structure.
#' @param transform_type type of covariate transformation:
#'   - "strong": aggressive nonlinear transformations (original)
#'   - "mild": gentle transformations that preserve IF orthogonality
#'   - "none": no transformation (X_dagger = X)
#' @param outcome_type Outcome family: \code{"binary"} or \code{"continuous"}.
#' @param heterogeneity_type Source outcome-model heterogeneity setting.
#' @param shift_strength Positive multiplier controlling covariate shift.
#' @param dgp_type Data-generating process: \code{"face"} or \code{"roce"}.
#' @param ate_deviation Non-negative additive source treatment-shift deviation
#'   under the FACE DGP. It is a mean shift for Gaussian outcomes and a
#'   log-odds shift for binary outcomes, not generally the induced marginal
#'   TATE difference.
#' @param n_deviated_sites Number of leading deviated source sites under the
#'   FACE DGP.
#' @param effect_mod_strength Source-only treatment-effect-modification strength
#'   under the FACE DGP; zero recovers the standard DGP.
#' @param n_target Optional target-site size for explicit FACE-DGP allocation.
#' @param n_source_sizes Optional vector of source-site sizes for explicit
#'   FACE-DGP allocation.
#' @param warn_ignored Whether to warn when a valid argument is irrelevant to
#'   the selected DGP.
#' @return list with all generated data
#' @details
#' The choice of estimand_type affects how to validate SE estimates:
#' 
#' \strong{sample:}
#' - True value = mean(logistic(X_dagger,target^T α)) for realized sample target covariates
#' - Different simulations have different true values
#' - Correct comparison: SE vs SD(bias), where bias = estimate - true_value_i
#' - Coverage: should be approximately 95\% when SE is correctly estimated
#' 
#' \strong{superpopulation (default):}
#' - True value = E_X[logistic(X_dagger^T α)] (fixed across simulations)
#' - All simulations use the same true value
#' - Correct comparison: SE vs SD(estimate)
#' - Coverage: may be below 95\% due to additional variability from true value estimation
#' 
#' For most purposes, "superpopulation" provides clearer diagnostics for SE validation.
#'
#' The choice of site_allocation affects the difficulty of variance estimation:
#' - "model": realistic but can lead to imbalanced site sizes (harder for variance estimation)
#' - "uniform": simpler, more balanced, good for validating variance formulas
#' - "balanced": ensures equal sizes, eliminates site-size variability
#' - "independent": each site generates data independently (federated-realistic),
#'   covariate shift comes from different site-specific distributions
#' @export
generate_simulation_data <- function(n_total = NULL, K = 3, p = 4, config = "C1",
                                     estimand_type = "superpopulation",
                                     site_allocation = "model",
                                     transform_type = "mild",
                                     outcome_type = "binary",
                                     heterogeneity_type = "none",
                                     shift_strength = ROCE_SHIFT_STRENGTH_DEFAULT,
                                   # DGP selector: "face" (FACE negative-transfer DGP, default)
                                   # or "roce" (RoCE DGP). ate_deviation / n_deviated_sites
                                   # apply only when dgp_type = "face".
                                   dgp_type = "face",
                                   ate_deviation   = 0.0,
                                   n_deviated_sites = 0L,
                                   # Source-only effect modification (FACE DGP); 0 = standard DGP.
                                   effect_mod_strength = 0,
                                   # Explicit per-site sample sizes (FACE DGP only)
                                   n_target         = NULL,
                                   n_source_sizes   = NULL,
                                   warn_ignored = TRUE) {
  # Explicit per-site sample sizes are a FACE-DGP feature. Reject them for the
  # roce DGP (which controls site sizes via site_allocation) rather than
  # silently ignoring the request.
  if ((!is.null(n_source_sizes) || !is.null(n_target)) && dgp_type != "face") {
    stop("n_target / n_source_sizes (explicit per-site sizes) are supported only ",
         "for dgp_type = 'face'; the roce DGP sets site sizes via site_allocation.",
         call. = FALSE)
  }
  # In per-site mode, derive a concrete n_total (and K) so the centralized
  # validation and result bookkeeping see values consistent with the allocation.
  if (!is.null(n_source_sizes)) {
    site_sizes <- resolve_face_site_sizes(n_total, n_target, n_source_sizes, K)
    n_total    <- site_sizes$n_total
    K          <- site_sizes$K
  }

  # Validate all parameters via centralized function (single source of truth)
  validate_simulation_params(
    estimand_type    = estimand_type,
    site_allocation  = site_allocation,
    transform_type   = transform_type,
    outcome_type     = outcome_type,
    heterogeneity_type = heterogeneity_type,
    shift_strength   = shift_strength,
    n_total          = n_total,
    K                = K,
    p                = p,
    config           = config,
    dgp_type         = dgp_type,
    ate_deviation    = ate_deviation,
    n_deviated_sites = n_deviated_sites,
    warn_ignored     = warn_ignored
  )

  # ---- FACE paper DGP (Han et al., JASA 2023, Section 5.1) ----
  # Self-contained early return: all data generation for this DGP is handled
  # here so that the existing RoCE code below is not affected in any way.
  if (dgp_type == "face") {
    return(
      generate_face_data(
        n_total          = n_total,
        K                = K,
        p                = p,
        config           = config,
        estimand_type    = estimand_type,
        outcome_type     = outcome_type,
        ate_deviation    = ate_deviation,
        n_deviated_sites = as.integer(n_deviated_sites),
        effect_mod_strength = effect_mod_strength,
        n_target         = n_target,
        n_source_sizes   = n_source_sizes
      )
    )
  }

  # Generate gamma parameters with FIXED seed (must match calculate_superpopulation_truth)
  # Uses with_seed() to isolate RNG so subsequent calls continue from caller's stream.
  gamma_params <- with_seed(99999,
    generate_site_model_parameters(K, p, shift_strength = shift_strength)
  )
  
  # ---- Covariate generation (depends on site_allocation) ----
  if (site_allocation == "independent") {
    # Independent per-site data generation:
    # Each site generates its own covariates from its own distribution
    # Target: X ~ N(0, I_p)   (baseline distribution)
    # Source j: X ~ N(μ_j, I_p)  (shifted distribution)
    # This mimics federated settings where sites have genuinely separate populations.
    
    # Allocate sample sizes (equal split, remainder goes to target)
    n_per_site <- floor(n_total / (K + 1))
    n_target <- n_total - K * n_per_site
    
    # Generate deterministic site-specific shift vectors (fixed across simulations)
    # Uses a different seed (99998) from gamma_params (99999) to avoid correlation.
    # Shift magnitude controlled by shift_strength (INDEPENDENT_SHIFT_SD from constants.R):
    #   shift_strength = 0.5 -> default moderate shift
    #   shift_strength = 1.0 -> stronger shift
    #   shift_strength = 2.0 -> stress-test shift
    site_shifts <- with_seed(99998, {
      shifts_list <- list()
      for (j in 1:K) {
        delta_j <- rep(0, p)
        n_shift_dims <- min(4, p)
        delta_j[1:n_shift_dims] <- rnorm(n_shift_dims, 0, INDEPENDENT_SHIFT_SD * shift_strength)
        shifts_list[[j]] <- delta_j
      }
      shifts_list
    })
    
    # Target site: X_t ~ N(0, I_p), clip to [-2.5, 2.5]
    X_target <- matrix(rnorm(n_target * p), nrow = n_target, ncol = p)
    X_target <- apply(X_target, 2, clip_to_range)
    colnames(X_target) <- paste0("X", 1:p)
    
    # Source sites: X_j ~ N(μ_j, I_p), clip to [-2.5, 2.5]
    X_sources <- list()
    for (j in 1:K) {
      X_j <- matrix(rnorm(n_per_site * p), nrow = n_per_site, ncol = p)
      X_j <- sweep(X_j, 2, site_shifts[[j]], "+")  # Add site-specific shift
      X_j <- apply(X_j, 2, clip_to_range)
      colnames(X_j) <- paste0("X", 1:p)
      X_sources[[j]] <- X_j
    }
    
    # Stack all covariates and assign site labels
    X <- rbind(X_target, do.call(rbind, X_sources))
    R <- c(rep("t", n_target),
           unlist(lapply(1:K, function(j) rep(paste0("s", j), n_per_site))))
    
  } else {
    # Shared-pool generation: all subjects drawn from common distribution
    X <- generate_base_covariates(n_total, p)
  }

  # Generate transformed covariates
  X_dagger <- transform_covariates(X, transform_type = transform_type)

  # True RoCE DGP bases. These stay fixed across C1-C4 so configurations
  # only change the fitted working bases exposed to estimators.
  Z_site_true <- X_dagger
  W_outcome_true <- X_dagger

  # Determine fitted site basis based on config.
  # C1, C2: correctly specified (uses X_dagger)
  # C3, C4: misspecified (uses X)
  if (config %in% c("C1", "C2")) {
    Z_site <- X_dagger
  } else {
    Z_site <- X
  }

  # Determine fitted outcome basis based on config.
  # C1, C3: correctly specified (uses X_dagger)
  # C2, C4: misspecified (uses X)
  if (config %in% c("C1", "C3")) {
    W_outcome <- X_dagger
  } else {
    W_outcome <- X
  }

  # ---- Site membership & treatment assignment ----
  # All allocation methods share the same covariate-dependent treatment propensity
  # P(A=1|X) derived from the multinomial γ parameters. Only the mechanism that
  # assigns subjects to sites differs. This keeps the confounding structure
  # constant so that differences in estimation accuracy are attributable solely
  # to site-size proportions, not to changes in the treatment mechanism.
  site_probs <- calculate_site_probabilities(Z_site_true, gamma_params, K)
  p_treat_true <- calculate_treatment_propensity(site_probs, K, n_total)

  if (site_allocation == "model") {
    # Multinomial logistic: site & treatment jointly sampled from the model
    assignments <- generate_site_treatment_assignments(n_total, site_probs, K)
    R <- assignments$R
    A <- assignments$A
  } else if (site_allocation == "independent") {
    # Independent allocation: R already set above during covariate generation.
    # Treatment assigned via the SAME covariate-dependent propensity formula
    # to maintain a consistent confounding structure across allocation methods.
    A <- rbinom(n_total, 1, p_treat_true)
  } else {
    # --- site assignment by specified method ---
    if (site_allocation == "uniform") {
      site_ids <- sample(0:K, n_total, replace = TRUE)
    } else if (site_allocation == "balanced") {
      site_ids <- sample(rep(0:K, length.out = n_total))
    } else if (site_allocation == "target_heavy") {
      site_ids <- sample(0:K, n_total, replace = TRUE,
                         prob = c(0.40, rep(0.60 / K, K)))
    } else if (site_allocation == "source_heavy") {
      site_ids <- sample(0:K, n_total, replace = TRUE,
                         prob = c(0.10, rep(0.90 / K, K)))
    }
    R <- ifelse(site_ids == 0, "t", paste0("s", site_ids))

    # Treatment from covariate-dependent propensity (same formula as model allocation)
    A <- rbinom(n_total, 1, p_treat_true)
  }

  # True outcome parameters (single source of truth)
  alphas <- get_true_outcome_parameters(p)
  alpha1 <- alphas$alpha1
  alpha0 <- alphas$alpha0

  # Generate potential outcomes for coverage evaluation, then select the
  # observed outcome so consistency holds observation by observation.
  # Y_1: potential outcome under treatment A=1
  # Y_0: potential outcome under treatment A=0
  A_ones <- rep(1, n_total)
  A_zeros <- rep(0, n_total)

  Y_1 <- generate_outcomes(W_outcome_true, A_ones, R, alpha1, alpha0,
                           noise_sd = ROCE_NOISE_SD_DEFAULT,
                           outcome_type = outcome_type,
                           heterogeneity_type = heterogeneity_type)
  Y_0 <- generate_outcomes(W_outcome_true, A_zeros, R, alpha1, alpha0,
                           noise_sd = ROCE_NOISE_SD_DEFAULT,
                           outcome_type = outcome_type,
                           heterogeneity_type = heterogeneity_type)
  Y <- ifelse(A == 1, Y_1, Y_0)
  
  # Calculate true potential outcome means
  # The interpretation depends on estimand_type:
  #   - "sample": sample-specific E_n[Y(1)] based on realized target covariates
  #   - "superpopulation": fixed E_X[Y(1)] = E[logistic(X_dagger^T α)] over covariate distribution
  target_idx <- R == "t"
  
  if (estimand_type == "sample") {
    # Sample-specific true values (original behavior)
    # These vary across simulations due to different realized covariates
    # mu1_sample = (1/n_t) * sum_{i in target} logistic(X_i^T α_1)
    W_target <- W_outcome_true[target_idx, , drop = FALSE]
    W_target_design <- cbind(1, W_target)
    
    eta1_target <- as.numeric(W_target_design %*% alpha1)
    eta0_target <- as.numeric(W_target_design %*% alpha0)
    
    mu1_true <- mean(logistic(eta1_target))  # E_n[logistic(X_dagger^T α_1)]
    mu0_true <- mean(logistic(eta0_target))  # E_n[logistic(X_dagger^T α_0)]
    
    # For continuous outcomes, truth is identity link (no logistic)
    if (outcome_type == "continuous") {
      mu1_true <- mean(eta1_target)
      mu0_true <- mean(eta0_target)
    }
    
    # Also compute realized potential outcomes for reference
    mu1_realized <- mean(Y_1[target_idx])  # Sample mean of realized Y(1)
    mu0_realized <- mean(Y_0[target_idx])  # Sample mean of realized Y(0)
    
  } else {
    # Superpopulation true values (fixed across all simulations)
    # Now computes TARGET-SPECIFIC superpopulation E_t[Y(1)] instead of overall E[Y(1)]
    superpop_truth <- calculate_superpopulation_truth(
      p = p, 
      K = K,
      config = config, 
      transform_type = transform_type,
      site_allocation = site_allocation,
      outcome_type = outcome_type,
      shift_strength = shift_strength
    )
    mu1_true <- superpop_truth$mu1_superpop
    mu0_true <- superpop_truth$mu0_superpop

    # Also compute sample-specific values for reference
    W_target <- W_outcome_true[target_idx, , drop = FALSE]
    W_target_design <- cbind(1, W_target)
    eta1_target <- as.numeric(W_target_design %*% alpha1)
    eta0_target <- as.numeric(W_target_design %*% alpha0)
    mu1_realized <- mean(logistic(eta1_target))
    mu0_realized <- mean(logistic(eta0_target))
    if (outcome_type == "continuous") {
      mu1_realized <- mean(eta1_target)
      mu0_realized <- mean(eta0_target)
    }
  }
  
  # Create data list
  data <- list(
    n = n_total,
    K = K,
    p = p,
    X = X,
    X_dagger = X_dagger,
    R = R,
    A = A,
    Y = Y,
    Y_1 = Y_1,  # Potential outcome under treatment
    Y_0 = Y_0,  # Potential outcome under control
    p_treat_true = p_treat_true,  # True P(A=1|X) used by the DGP
    mu1_true = mu1_true,  # True E[Y(1)] (interpretation depends on estimand_type)
    mu0_true = mu0_true,  # True E[Y(0)]
    mu1_realized = mu1_realized,  # Sample-specific mean (for reference)
    mu0_realized = mu0_realized,  # Sample-specific mean (for reference)
    gamma_params = gamma_params,
    alpha1_true = alpha1,
    alpha0_true = alpha0,
    config = config,
    Z_site_true = Z_site_true,
    W_outcome_true = W_outcome_true,
    Z_site = Z_site,
    W_outcome = W_outcome,
    estimand_type = estimand_type,  # Record which estimand type was used
    outcome_type = outcome_type,    # Record outcome type
    dgp_type = "roce"              # Record DGP type
  )
  
  return(data)
}

#' Split data by site
#' @param data full dataset
#' @return list with data for each site
#' @export
split_data_by_site <- function(data) {
  site_data <- list()
  
  # Get unique sites
  sites <- unique(data$R)
  
  for (site in sites) {
    idx <- which(data$R == site)
    site_data[[site]] <- list(
      n = length(idx),
      X = data$X[idx, , drop = FALSE],
      X_dagger = data$X_dagger[idx, , drop = FALSE],
      A = data$A[idx],
      Y = data$Y[idx],
      Z_site_true = (data$Z_site_true %||% data$Z_site)[idx, , drop = FALSE],
      W_outcome_true = (data$W_outcome_true %||% data$W_outcome)[idx, , drop = FALSE],
      Z_site = data$Z_site[idx, , drop = FALSE],
      W_outcome = data$W_outcome[idx, , drop = FALSE]
    )
  }
  
  return(site_data)
}
