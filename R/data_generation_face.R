# data_generation_face.R - FACE Paper DGP (Han et al., JASA 2023, Section 5.1)
#
# This file implements the data-generating process described in the published
# FACE paper. The key differences from the FACE-C DGP (in data_generation.R) are:
#   - Covariates: site-specific skewed-normal distributions (no pooled draw)
#   - Outcome:    continuous, linear + squared covariate terms, site-specific
#                 constant ATE Δ_k (no within-site effect modification)
#   - PS model:   logistic with both linear and squared covariate terms
#   - Config C1–C4: correctness/misspecification refers to including vs.
#                   omitting the squared covariate block
#
# Functions:
#   - generate_skewed_normal
#   - generate_face_covariates
#   - calculate_face_propensity
#   - build_face_ate_map
#   - generate_face_outcomes
#   - calculate_face_truth
#   - generate_face_data  (main orchestrator)

#' Generate samples from the skewed-normal distribution SN(κ, φ², ν)
#'
#' Uses the standard stochastic representation (Azzalini 1985):
#' \deqn{X = \kappa + \phi \left(\delta |Z_1| + \sqrt{1-\delta^2}\, Z_2\right)}
#' where \eqn{Z_1, Z_2 \sim N(0,1)} independently and \eqn{\delta = \nu/\sqrt{1+\nu^2}}.
#' When \eqn{\nu = 0} (target site) this reduces to \eqn{X \sim N(\kappa, \phi^2)}.
#'
#' @param n   Number of observations.
#' @param kappa Location parameter (scalar).
#' @param phi   Scale parameter (standard deviation, > 0).
#' @param nu    Skewness parameter (\eqn{\nu = 0} gives the symmetric normal).
#' @return Numeric vector of length n.
generate_skewed_normal <- function(n, kappa = 0, phi = 1, nu = 0) {
  delta <- nu / sqrt(1 + nu^2)
  Z1    <- rnorm(n)
  Z2    <- rnorm(n)
  kappa + phi * (delta * abs(Z1) + sqrt(1 - delta^2) * Z2)
}

#' Generate site-specific covariates for the FACE paper DGP
#'
#' Each site draws its own covariates from a skewed-normal distribution:
#' \itemize{
#'   \item Target (\eqn{R = t}):
#'         \eqn{X_{tp} \sim \text{SN}(\kappa, 1, 0) = N(\kappa, 1)}
#'   \item Source \eqn{j} (\eqn{R = s_j}):
#'         \eqn{X_{jp} \sim \text{SN}(\kappa, 1, \nu_j)},
#'         with \eqn{\nu_j} equally spaced on
#'         \eqn{[\nu_{\max}/K,\; \nu_{\max}]}
#' }
#' All sites share the same location \eqn{\kappa} and unit scale, consistent
#' with \eqn{\kappa_{k\cdot} \in (0.10, 0.15)} and \eqn{\phi_{k\cdot} = 1}
#' in the paper.
#'
#' @param n_target          Number of target-site observations.
#' @param n_source_per_site Number of observations per source site.
#' @param K                 Number of source sites (0 for reference-population calls).
#' @param p                 Number of covariates.
#' @param kappa             Location parameter (default \code{FACE_KAPPA}).
#' @param nu_source_max     Maximum skewness for source sites
#'                          (default \code{FACE_NU_SOURCE_MAX}).
#' @return Named list with:
#'   \itemize{
#'     \item \code{X}         – (\eqn{n_{\text{total}} \times p}) covariate matrix.
#'     \item \code{R}         – character vector of site labels.
#'     \item \code{kappa}     – location parameter used.
#'     \item \code{nu_source} – skewness values assigned to each source site.
#'   }
generate_face_covariates <- function(n_target, n_source_per_site, K, p,
                                           kappa         = FACE_KAPPA,
                                           nu_source_max = FACE_NU_SOURCE_MAX) {
  n_total <- n_target + K * n_source_per_site
  X <- matrix(0, nrow = n_total, ncol = p)
  R <- character(n_total)

  # Source-site skewness: equally spaced on [ν_max / K, ν_max]
  # (Avoids ν = 0 for source sites so they are always distinguishable from target.)
  nu_source <- if (K > 0) seq(nu_source_max / K, nu_source_max, length.out = K) else numeric(0)

  # Target site (ν = 0 → symmetric, matches the paper specification)
  if (n_target > 0) {
    for (j in 1:p) {
      X[1:n_target, j] <- generate_skewed_normal(n_target, kappa = kappa, phi = 1, nu = 0)
    }
    R[1:n_target] <- "t"
  }

  # Source sites
  for (k in seq_len(K)) {
    idx_start <- n_target + (k - 1L) * n_source_per_site + 1L
    idx_end   <- n_target + k * n_source_per_site
    for (j in 1:p) {
      X[idx_start:idx_end, j] <- generate_skewed_normal(
        n_source_per_site, kappa = kappa, phi = 1, nu = nu_source[k]
      )
    }
    R[idx_start:idx_end] <- paste0("s", k)
  }

  colnames(X) <- paste0("X", 1:p)
  list(X = X, R = R, kappa = kappa, nu_source = nu_source)
}

#' Compute treatment propensity scores for the FACE paper DGP
#'
#' \deqn{\pi_k = \text{expit}(X \alpha_1 + X^{\circ 2} \alpha_2)}
#' where \eqn{X^{\circ 2}} denotes element-wise squaring.
#' Output is clipped to [\code{POSITIVITY_LOWER}, \code{POSITIVITY_UPPER}].
#'
#' @param X      Covariate matrix (\eqn{n \times p}).
#' @param alpha1 Linear PS coefficients (length p).
#' @param alpha2 Squared PS coefficients (length p).
#' @return Numeric vector of propensity scores (length n).
calculate_face_propensity <- function(X, alpha1, alpha2) {
  X_sq  <- X^2
  eta   <- as.numeric(X %*% alpha1 + X_sq %*% alpha2)
  p_treat <- logistic(eta)
  pmax(pmin(p_treat, POSITIVITY_UPPER), POSITIVITY_LOWER)
}

#' Build a site-to-ATE mapping for the FACE paper DGP
#'
#' The target site always has \eqn{\Delta_T = } \code{FACE_ATE_TARGET}.
#' The first \code{n_deviated} source sites receive
#' \eqn{\Delta_k = \Delta_T + \texttt{ate\_deviation}};
#' the remaining source sites are informative (\eqn{\Delta_k = \Delta_T}).
#' Setting \code{n_deviated = 0} (default) makes all source sites informative
#' (deviation level 1 in Table 1 of the paper).
#'
#' @param K             Number of source sites.
#' @param ate_deviation Additive ATE deviation for non-informative sites
#'                      (default 0.0 → all informative).
#' @param n_deviated    Number of leading source sites that deviate
#'                      (default 0L → none deviate).
#' @return Named numeric vector mapping site label → ATE.
build_face_ate_map <- function(K, ate_deviation = 0.0, n_deviated = 0L) {
  ate_map      <- c("t" = FACE_ATE_TARGET)
  source_ates  <- rep(FACE_ATE_TARGET, K)
  if (n_deviated > 0L) {
    deviated_idx        <- seq_len(min(n_deviated, K))
    source_ates[deviated_idx] <- FACE_ATE_TARGET + ate_deviation
  }
  names(source_ates) <- paste0("s", seq_len(K))
  c(ate_map, source_ates)
}

#' Generate potential outcomes for the FACE paper DGP
#'
#' **Continuous** (\code{outcome_type = "continuous"}, default):
#' \deqn{Y_k(a) = (X_k - \kappa)^\top \beta_{\text{lin}}
#'               + (X_k^{\circ 2})^\top \beta_{\text{sq}}
#'               + \Delta_k \cdot \mathbf{1}(a=1) + \varepsilon_k,
#'               \quad \varepsilon_k \sim N(0, \sigma^2)}
#'
#' **Binary** (\code{outcome_type = "binary"}):
#' \deqn{Y_k(a) \sim \text{Bernoulli}\bigl(
#'       \text{expit}\bigl(
#'         (X_k - \kappa)^\top \beta_{\text{lin}}
#'         + (X_k^{\circ 2})^\top \beta_{\text{sq}}
#'         + \Delta_k \cdot \mathbf{1}(a=1)
#'       \bigr)\bigr)}
#'
#' Both potential-outcome arms share the same covariate coefficients
#' (\eqn{\beta_{1,a} = \beta_{2,a}} for a = 0, 1, matching the paper); the
#' treatment effect is captured entirely by the site-specific constant \eqn{\Delta_k}.
#'
#' @param X            Covariate matrix (\eqn{n \times p}).
#' @param A            Treatment indicator vector (0/1, length n).
#' @param R            Site membership character vector (length n).
#' @param beta_lin     Linear covariate coefficients (length p).
#' @param beta_sq      Squared covariate coefficients (length p).
#' @param ate_map      Named numeric vector (from \code{build_face_ate_map()}).
#' @param kappa        Centering constant (scalar, default \code{FACE_KAPPA}).
#' @param noise_sd     Noise standard deviation (default \code{FACE_NOISE_SD});
#'                     only used when \code{outcome_type = "continuous"}.
#' @param outcome_type "continuous" (default) or "binary".
#' @return Numeric vector of outcomes (length n).
generate_face_outcomes <- function(X, A, R, beta_lin, beta_sq, ate_map,
                                         kappa        = FACE_KAPPA,
                                         noise_sd     = FACE_NOISE_SD,
                                         outcome_type = "continuous") {
  X_centered <- sweep(X, 2L, kappa, "-")
  X_sq       <- X^2
  eta        <- as.numeric(X_centered %*% beta_lin + X_sq %*% beta_sq)
  delta      <- ate_map[R]          # site-specific ATE, same length as R
  linpred    <- eta + delta * A

  if (outcome_type == "binary") {
    prob <- logistic(linpred)
    rbinom(length(A), 1L, prob)
  } else {
    linpred + rnorm(length(A), 0, noise_sd)
  }
}

#' Compute superpopulation truth for the FACE paper DGP (target-population specific)
#'
#' Evaluates \eqn{E_t[Y(a)]} via Monte Carlo integration over the target-site
#' covariate distribution (\eqn{X \sim N(\kappa, I_p)}, \eqn{\nu = 0}).
#'
#' Analytically, for the target site:
#' \deqn{E_t[Y(1)] = \sum_{j \le 4} \beta_{\text{sq},j}(1 + \kappa^2) + \Delta_T}
#' \deqn{E_t[Y(0)] = \sum_{j \le 4} \beta_{\text{sq},j}(1 + \kappa^2)}
#' because \eqn{E[X_j - \kappa] = 0} and \eqn{E[X_j^2] = 1 + \kappa^2} for
#' \eqn{X_j \sim N(\kappa, 1)}.
#' The MC approach is used for consistency with \code{calculate_superpopulation_truth()}.
#'
#' @param p           Number of covariates.
#' @param kappa       Location parameter (default \code{FACE_KAPPA}).
#' @param outcome_type "continuous" or "binary".
#' @param n_ref       Reference population size (default 100000).
#' @param ref_seed    Fixed RNG seed for reproducibility (default 99999).
#' @return Named list with \code{mu1_superpop}, \code{mu0_superpop},
#'   \code{ate_superpop} (always \code{FACE_ATE_TARGET} by construction),
#'   \code{n_ref}, and \code{ref_seed}.
calculate_face_truth <- function(p, kappa = FACE_KAPPA,
                                       outcome_type = "continuous",
                                       n_ref = 100000L, ref_seed = 99999L) {
  with_seed(ref_seed, {
    # Target-site covariates: X ~ N(kappa, 1) (ν = 0 → symmetric)
    ref_covs <- generate_face_covariates(
      n_target          = n_ref,
      n_source_per_site = 0L,
      K                 = 0L,
      p                 = p,
      kappa             = kappa
    )
    X_ref      <- ref_covs$X
    out_params <- get_face_outcome_parameters(p)

    X_centered <- sweep(X_ref, 2L, kappa, "-")
    X_sq       <- X_ref^2
    eta_ref    <- as.numeric(X_centered %*% out_params$beta_linear +
                             X_sq       %*% out_params$beta_squared)

    if (outcome_type == "binary") {
      # Binary: E[Y(a)] = E[expit(η + Δ_T · a)]
      mu0_superpop <- mean(logistic(eta_ref))
      mu1_superpop <- mean(logistic(eta_ref + FACE_ATE_TARGET))
    } else {
      # Continuous: E[Y(a)] = E[η] + Δ_T · a
      mu0_superpop <- mean(eta_ref)
      mu1_superpop <- mu0_superpop + FACE_ATE_TARGET
    }

    list(
      mu1_superpop = mu1_superpop,
      mu0_superpop = mu0_superpop,
      ate_superpop = mu1_superpop - mu0_superpop,
      n_ref        = n_ref,
      ref_seed     = ref_seed
    )
  })
}

#' Orchestrate full data generation for the FACE paper DGP
#'
#' Creates a dataset whose structure (field names and types) is identical to that
#' produced by the FACE-C DGP path of \code{generate_simulation_data()}, so that
#' downstream code (\code{split_data_by_site()}, estimation algorithms, etc.)
#' requires no modification.
#'
#' **Config semantics for the FACE paper DGP:**
#' \itemize{
#'   \item C1 (both correctly specified): PS uses \eqn{[X, X^2]};
#'         OR uses \eqn{[X - \kappa, X^2]}.
#'   \item C2 (OR misspecified): PS uses \eqn{[X, X^2]};
#'         OR uses only \eqn{X} (drops quadratic block).
#'   \item C3 (PS misspecified): PS uses only \eqn{X};
#'         OR uses \eqn{[X - \kappa, X^2]}.
#'   \item C4 (both misspecified): PS uses only \eqn{X};
#'         OR uses only \eqn{X}.
#' }
#'
#' @param n_total          Total observations across all sites.
#' @param K                Number of source sites.
#' @param p                Number of covariates (only first 4 have non-zero coefficients).
#' @param config           Configuration string: "C1", "C2", "C3", or "C4".
#' @param estimand_type    "superpopulation" (fixed truth) or "sample" (realized truth).
#' @param outcome_type     "continuous" (default) or "binary".
#' @param ate_deviation    ATE deviation for non-informative source sites (default 0).
#' @param n_deviated_sites Number of leading source sites with deviated ATE (default 0L).
#' @return Named list with the same fields as the FACE-C DGP output.
generate_face_data <- function(n_total, K = 3, p = 4, config = "C1",
                                     estimand_type    = "superpopulation",
                                     outcome_type     = "continuous",
                                     ate_deviation    = 0.0,
                                     n_deviated_sites = 0L) {
  # ---- 1. Allocate equal sample sizes (remainder goes to target) ----
  n_source_per_site <- floor(n_total / (K + 1L))
  n_target          <- n_total - K * n_source_per_site

  # ---- 2. Generate site-specific covariates from skewed-normal distribution ----
  cov_list <- generate_face_covariates(
    n_target          = n_target,
    n_source_per_site = n_source_per_site,
    K                 = K,
    p                 = p
  )
  X  <- cov_list$X
  R  <- cov_list$R
  X_dagger <- X   # No nonlinear transform in this DGP; keep for structural consistency

  # ---- 3. Build design matrices for PS and Outcome (config-dependent) ----
  X_sq <- X^2
  kappa <- cov_list$kappa

  # Z_site: used to fit the site / propensity model
  #   C1, C2 (correctly specified): include squared block
  #   C3, C4 (misspecified):        linear-only (drop squared block)
  Z_site <- if (config %in% c("C1", "C2")) cbind(X, X_sq) else X

  # W_outcome: used to fit the outcome regression
  #   C1, C3 (correctly specified): [X - kappa, X^2]
  #   C2, C4 (misspecified):        X only (drop squared block)
  X_centered <- sweep(X, 2L, kappa, "-")
  W_outcome  <- if (config %in% c("C1", "C3")) cbind(X_centered, X_sq) else X

  # ---- 4. Propensity scores and treatment assignment ----
  ps_params <- get_face_ps_parameters(p)
  p_treat   <- calculate_face_propensity(X, ps_params$alpha1, ps_params$alpha2)
  A         <- rbinom(n_total, 1L, p_treat)

  # ---- 5. Site-specific ATEs and outcome generation ----
  out_params <- get_face_outcome_parameters(p)
  ate_map    <- build_face_ate_map(K, ate_deviation, n_deviated_sites)

  Y <- generate_face_outcomes(
    X, A, R,
    beta_lin     = out_params$beta_linear,
    beta_sq      = out_params$beta_squared,
    ate_map      = ate_map,
    kappa        = kappa,
    outcome_type = outcome_type
  )

  # Potential outcomes (needed for coverage evaluation and ATE estimation)
  A_ones  <- rep(1L, n_total)
  A_zeros <- rep(0L, n_total)
  Y_1 <- generate_face_outcomes(X, A_ones,  R, out_params$beta_linear,
                                      out_params$beta_squared, ate_map, kappa,
                                      outcome_type = outcome_type)
  Y_0 <- generate_face_outcomes(X, A_zeros, R, out_params$beta_linear,
                                      out_params$beta_squared, ate_map, kappa,
                                      outcome_type = outcome_type)

  # ---- 6. True potential-outcome means ----
  target_idx <- R == "t"

  if (estimand_type == "sample") {
    # Sample-specific truth: computed from the realized target observations
    X_centered_target <- sweep(X[target_idx, , drop = FALSE], 2L, kappa, "-")
    X_sq_target       <- X[target_idx, , drop = FALSE]^2
    eta_target <- as.numeric(
      X_centered_target %*% out_params$beta_linear +
      X_sq_target       %*% out_params$beta_squared
    )
    if (outcome_type == "binary") {
      mu0_true <- mean(logistic(eta_target))
      mu1_true <- mean(logistic(eta_target + FACE_ATE_TARGET))
    } else {
      mu0_true <- mean(eta_target)
      mu1_true <- mu0_true + FACE_ATE_TARGET
    }
    mu1_realized <- mean(Y_1[target_idx])
    mu0_realized <- mean(Y_0[target_idx])

  } else {
    # Superpopulation truth: fixed across simulations (Monte Carlo integration)
    superpop_truth <- calculate_face_truth(p = p, kappa = kappa,
                                                 outcome_type = outcome_type)
    mu1_true     <- superpop_truth$mu1_superpop
    mu0_true     <- superpop_truth$mu0_superpop
    mu1_realized <- mean(Y_1[target_idx])   # sample-specific (for reference logging)
    mu0_realized <- mean(Y_0[target_idx])
  }

  # ---- 7. Assemble output list (matches FACE-C DGP structure) ----
  list(
    n            = n_total,
    K            = K,
    p            = p,
    X            = X,
    X_dagger     = X_dagger,
    R            = R,
    A            = A,
    Y            = Y,
    Y_1          = Y_1,
    Y_0          = Y_0,
    p_treat_true = p_treat,
    mu1_true     = mu1_true,
    mu0_true     = mu0_true,
    mu1_realized = mu1_realized,
    mu0_realized = mu0_realized,
    gamma_params = NULL,           # not applicable for face DGP
    alpha1_true   = NULL,           # not applicable (use out_params instead)
    alpha0_true   = NULL,
    out_params   = out_params,     # face-specific: linear + squared coefficients
    ate_map      = ate_map,        # face-specific: site → ATE lookup
    config       = config,
    Z_site       = Z_site,
    W_outcome    = W_outcome,
    estimand_type = estimand_type,
    outcome_type  = outcome_type,
    dgp_type      = "face"
  )
}
