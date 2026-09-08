# data_generation_face.R - FACE Paper DGP (Han et al., JASA 2023, Section 5.1)
#
# This file implements the data-generating process described in the published
# FACE paper. The key differences from the RoCE DGP (in data_generation.R) are:
#   - Covariates: site-specific skewed-normal distributions (no pooled draw)
#   - Outcome:    continuous or binary, with linear + squared covariate terms
#                 and a site-specific constant treatment shift Δ_k. The shift
#                 equals an ATE only for the continuous identity-link outcome;
#                 for binary outcomes it is a conditional log-odds shift.
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
#   - resolve_face_site_sizes
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
#' @param n_target       Number of target-site observations.
#' @param n_source_sizes Integer vector of per-site source sample sizes; its
#'                       length determines the number of source sites \eqn{K}
#'                       (use \code{integer(0)} for reference-population calls).
#' @param p              Number of covariates.
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
generate_face_covariates <- function(n_target, n_source_sizes, p,
                                           kappa         = FACE_KAPPA,
                                           nu_source_max = FACE_NU_SOURCE_MAX) {
  K       <- length(n_source_sizes)
  n_total <- n_target + sum(n_source_sizes)
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

  # Source sites (sizes may differ across sites; advance a running row offset)
  offset <- n_target
  for (k in seq_len(K)) {
    n_k <- n_source_sizes[k]
    if (n_k <= 0L) next
    idx <- (offset + 1L):(offset + n_k)
    for (j in 1:p) {
      X[idx, j] <- generate_skewed_normal(
        n_k, kappa = kappa, phi = 1, nu = nu_source[k]
      )
    }
    R[idx] <- paste0("s", k)
    offset <- offset + n_k
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

#' Build a site-to-treatment-shift mapping for the FACE paper DGP
#'
#' The historical function name is retained for API compatibility. The target
#' site always has \eqn{\Delta_T = } \code{FACE_ATE_TARGET} for continuous
#' outcomes or \code{FACE_BINARY_ATE_TARGET} on the binary log-odds scale.
#' The first \code{n_deviated} source sites receive
#' \eqn{\Delta_k = \Delta_T + \texttt{ate\_deviation}};
#' the remaining source sites are informative (\eqn{\Delta_k = \Delta_T}).
#' Setting \code{n_deviated = 0} (default) makes all source sites informative
#' (deviation level 1 in Table 1 of the paper).
#'
#' @param K             Number of source sites.
#' @param ate_deviation Additive treatment-shift deviation for non-informative sites
#'                      (default 0.0 → all informative). On the natural scale of
#'                      \code{base_ate} (mean shift for continuous, log-odds shift
#'                      for binary).
#' @param n_deviated    Number of leading source sites that deviate
#'                      (default 0L → none deviate).
#' @param base_ate      Target / informative-source treatment shift on the
#'                      outcome linear-predictor scale: \code{FACE_ATE_TARGET}
#'                      (continuous mean shift, default) or
#'                      \code{FACE_BINARY_ATE_TARGET} (binary log-odds shift).
#' @return Named numeric vector mapping site label to treatment shift.
build_face_ate_map <- function(K, ate_deviation = 0.0, n_deviated = 0L,
                               base_ate = FACE_ATE_TARGET) {
  ate_map      <- c("t" = base_ate)
  source_ates  <- rep(base_ate, K)
  if (n_deviated > 0L) {
    deviated_idx        <- seq_len(min(n_deviated, K))
    source_ates[deviated_idx] <- base_ate + ate_deviation
  }
  names(source_ates) <- paste0("s", seq_len(K))
  c(ate_map, source_ates)
}

#' Logit-scale calibration for the binary FACE paper DGP
#'
#' The raw linear predictor \eqn{\eta = (X-\kappa)^\top\beta_{\text{lin}} +
#' (X^{\circ 2})^\top\beta_{\text{sq}}} has a target-law standard deviation of
#' \eqn{\approx 3.2}, which saturates \code{expit()} (treated-arm prevalence
#' \eqn{\approx 0.98}, a near-degenerate outcome regression). For the binary
#' outcome the covariate signal is standardized to a controlled logit-scale
#' spread (\code{FACE_BINARY_SIGNAL_SD}) via \code{face_binary_logit()}.
#'
#' The centring / scaling moments are computed \strong{once} from the fixed
#' \eqn{\nu = 0} target reference population — the same population (and RNG seed)
#' used by \code{calculate_face_truth()} — so the realized data and both truth
#' branches share one deterministic transform. They must never be recomputed from
#' realized (mixed-site) covariates, which would silently shift the estimand.
#'
#' @inheritParams calculate_face_truth
#' @param signal_sd Target logit-scale sd of the covariate signal
#'   (default \code{FACE_BINARY_SIGNAL_SD}).
#' @return Named list containing the calibration moments and the corresponding
#'   target-population binary potential-outcome means. The latter are reused by
#'   \code{calculate_face_truth()} so a simulation task does not regenerate the
#'   same deterministic reference population.
#' @keywords internal
get_face_binary_calibration <- function(p, kappa = FACE_KAPPA,
                                        signal_sd = FACE_BINARY_SIGNAL_SD,
                                        n_ref = 100000L, ref_seed = 99999L) {
  out_params <- get_face_outcome_parameters(p)
  eta_ref <- with_seed(ref_seed, {
    X_ref <- generate_face_covariates(n_target       = n_ref,
                                      n_source_sizes = integer(0),
                                      p              = p,
                                      kappa          = kappa)$X
    as.numeric(sweep(X_ref, 2L, kappa, "-") %*% out_params$beta_linear +
               (X_ref^2)                    %*% out_params$beta_squared)
  })
  calibration <- list(
    eta_mean = mean(eta_ref),
    eta_sd = max(stats::sd(eta_ref), EPSILON_DEFAULT),
    signal_sd = signal_sd
  )
  g_ref <- face_binary_logit(eta_ref, calibration)
  calibration$mu0_superpop <- mean(logistic(g_ref))
  calibration$mu1_superpop <- mean(
    logistic(g_ref + FACE_BINARY_ATE_TARGET)
  )
  calibration$n_ref <- as.integer(n_ref)
  calibration$ref_seed <- as.integer(ref_seed)
  calibration$p <- as.integer(p)
  calibration$kappa <- as.numeric(kappa)
  calibration
}

#' Map the raw FACE linear predictor onto the standardized binary logit scale
#'
#' \eqn{g(X) = \texttt{signal\_sd} \cdot (\eta - E_t[\eta]) / \mathrm{sd}_t[\eta]}.
#' Because this is affine in \eqn{\eta}, the conditional logit
#' \eqn{g(X) + \Delta\,\mathbf{1}(a=1)} is linear in \eqn{[1, X-\kappa, X^{\circ 2}]},
#' so the correctly-specified outcome regression (C1/C3) stays exactly correct.
#'
#' @param eta   Raw linear predictor (numeric vector).
#' @param calib Calibration list from \code{get_face_binary_calibration()}.
#' @return Standardized logit-scale signal (numeric vector, same length as
#'   \code{eta}).
#' @keywords internal
face_binary_logit <- function(eta, calib) {
  calib$signal_sd * (eta - calib$eta_mean) / calib$eta_sd
}

#' Generate potential outcomes for the FACE paper DGP
#'
#' \strong{Continuous} (\code{outcome_type = "continuous"}, default):
#' \deqn{Y_k(a) = (X_k - \kappa)^\top \beta_{\text{lin}}
#'               + (X_k^{\circ 2})^\top \beta_{\text{sq}}
#'               + \Delta_k \cdot \mathbf{1}(a=1) + \varepsilon_k,
#'               \quad \varepsilon_k \sim N(0, \sigma^2)}
#'
#' \strong{Binary} (\code{outcome_type = "binary"}): the raw linear predictor
#' \eqn{\eta} (which saturates \code{expit()}) is standardized to a controlled
#' logit-scale signal \eqn{g(X)} (see \code{face_binary_logit()}):
#' \deqn{Y_k(a) \sim \text{Bernoulli}\bigl(
#'       \text{expit}\bigl( g(X_k) + \Delta_k \cdot \mathbf{1}(a=1) \bigr)\bigr),
#'       \quad g(X) = \sigma_g \cdot (\eta - E_t[\eta]) / \mathrm{sd}_t[\eta].}
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
#' @param binary_calib Calibration list from \code{get_face_binary_calibration()};
#'                     required when \code{outcome_type = "binary"}, ignored
#'                     otherwise. Standardizes the logit-scale signal.
#' @param effect_mod_strength Source-only effect-modification strength; zero
#'   recovers the standard FACE DGP.
#' @param em_direction Optional covariate direction for effect modification.
#'   When omitted, a deterministic normalized direction is used.
#' @return Numeric vector of outcomes (length n).
generate_face_outcomes <- function(X, A, R, beta_lin, beta_sq, ate_map,
                                         kappa        = FACE_KAPPA,
                                         noise_sd     = FACE_NOISE_SD,
                                         outcome_type = "continuous",
                                         binary_calib = NULL,
                                         effect_mod_strength = 0,
                                         em_direction = NULL) {
  X_centered <- sweep(X, 2L, kappa, "-")
  X_sq       <- X^2
  eta        <- as.numeric(X_centered %*% beta_lin + X_sq %*% beta_sq)
  delta      <- ate_map[R]  # site-specific treatment shift, same length as R

  # Source-only effect modification. A covariate-dependent treatment-effect term
  # is added for source sites only and centered on the target covariate location
  # (kappa), so E_t[em_term] = 0. For the continuous identity-link DGP used by
  # the one-vs-two-round experiment, each source's target-population ATE is
  # unchanged while its conditional outcome *shape* differs from the target's.
  # Mean-zero logit-scale modification would not imply the same marginal risk
  # difference for a binary outcome. This is the regime that separates the one-round
  # (target-initialized) and two-round (source-initialized) estimators; the
  # default of 0 leaves the standard FACE DGP unchanged.
  tau <- delta
  if (effect_mod_strength != 0) {
    if (is.null(em_direction)) {
      em_direction <- numeric(ncol(X)); em_direction[1L] <- 1
    }
    is_source <- as.numeric(R != "t")
    tau <- delta + effect_mod_strength * is_source *
           as.numeric(X_centered %*% em_direction)
  }

  if (outcome_type == "binary") {
    if (is.null(binary_calib)) {
      stop("binary outcomes require 'binary_calib' from get_face_binary_calibration().",
           call. = FALSE)
    }
    prob <- logistic(face_binary_logit(eta, binary_calib) + tau * A)
    rbinom(length(A), 1L, prob)
  } else {
    eta + tau * A + rnorm(length(A), 0, noise_sd)
  }
}

#' Compute superpopulation truth for the FACE paper DGP (target-population specific)
#'
#' Evaluates \eqn{E_t[Y(a)]} via Monte Carlo integration over the target-site
#' covariate distribution (\eqn{X \sim N(\kappa, I_p)}, \eqn{\nu = 0}).
#'
#' \strong{Continuous:} \eqn{E_t[Y(a)] = E_t[\eta] + \Delta_T\,a}, so
#' \eqn{\text{ate\_superpop} = \Delta_T = } \code{FACE_ATE_TARGET}. Analytically
#' \eqn{E_t[\eta] = \sum_{j \le 4} \beta_{\text{sq},j}(1 + \kappa^2)} because
#' \eqn{E[X_j - \kappa] = 0} and \eqn{E[X_j^2] = 1 + \kappa^2}.
#'
#' \strong{Binary:} \eqn{E_t[Y(a)] = E_t[\text{expit}(g(X) + \Delta_{\text{bin}}\,a)]}
#' with \eqn{g} the standardized logit signal (\code{face_binary_logit()}); here
#' \eqn{\Delta_{\text{bin}}} is a log-odds shift, so \code{ate_superpop} is the
#' implied target-population \emph{risk difference} (\eqn{\approx 0.21}), not
#' \eqn{\Delta_{\text{bin}}}.
#' The MC approach is used for consistency with \code{calculate_superpopulation_truth()}.
#'
#' @param p           Number of covariates.
#' @param kappa       Location parameter (default \code{FACE_KAPPA}).
#' @param outcome_type "continuous" or "binary".
#' @param n_ref       Reference population size (default 100000).
#' @param ref_seed    Fixed RNG seed for reproducibility (default 99999).
#' @param binary_calib Optional result from
#'   \code{get_face_binary_calibration()} computed with the same reference
#'   settings. When supplied for a binary outcome, its stored truth moments are
#'   reused instead of drawing the deterministic reference population again.
#' @return Named list with \code{mu1_superpop}, \code{mu0_superpop},
#'   \code{ate_superpop} (continuous: \code{FACE_ATE_TARGET}; binary: the implied
#'   risk difference), \code{n_ref}, and \code{ref_seed}.
calculate_face_truth <- function(p, kappa = FACE_KAPPA,
                                       outcome_type = "continuous",
                                       n_ref = 100000L, ref_seed = 99999L,
                                       binary_calib = NULL) {
  if (outcome_type == "binary" && !is.null(binary_calib)) {
    required <- c(
      "mu0_superpop", "mu1_superpop", "n_ref", "ref_seed", "p", "kappa"
    )
    if (!is.list(binary_calib) ||
        any(!required %in% names(binary_calib)) ||
        !identical(as.integer(binary_calib$n_ref), as.integer(n_ref)) ||
        !identical(as.integer(binary_calib$ref_seed), as.integer(ref_seed)) ||
        !identical(as.integer(binary_calib$p), as.integer(p)) ||
        !isTRUE(all.equal(
          as.numeric(binary_calib$kappa), as.numeric(kappa),
          tolerance = 0, check.attributes = FALSE
        )) ||
        any(!is.finite(c(
          binary_calib$mu0_superpop, binary_calib$mu1_superpop
        )))) {
      stop(
        "calculate_face_truth: binary_calib does not match the requested reference population.",
        call. = FALSE
      )
    }
    return(list(
      mu1_superpop = as.numeric(binary_calib$mu1_superpop),
      mu0_superpop = as.numeric(binary_calib$mu0_superpop),
      ate_superpop = as.numeric(
        binary_calib$mu1_superpop - binary_calib$mu0_superpop
      ),
      n_ref = as.integer(n_ref),
      ref_seed = as.integer(ref_seed)
    ))
  }
  with_seed(ref_seed, {
    # Target-site covariates: X ~ N(kappa, 1) (ν = 0 → symmetric)
    ref_covs <- generate_face_covariates(
      n_target       = n_ref,
      n_source_sizes = integer(0),
      p              = p,
      kappa          = kappa
    )
    X_ref      <- ref_covs$X
    out_params <- get_face_outcome_parameters(p)

    X_centered <- sweep(X_ref, 2L, kappa, "-")
    X_sq       <- X_ref^2
    eta_ref    <- as.numeric(X_centered %*% out_params$beta_linear +
                             X_sq       %*% out_params$beta_squared)

    if (outcome_type == "binary") {
      # Binary: E[Y(a)] = E[expit(g(X) + Δ_bin · a)], g the standardized signal.
      # This reference (same n_ref / ref_seed) reproduces the calibration moments
      # baked into the data via get_face_binary_calibration(), so truth and data
      # share one transform by construction.
      calib <- list(eta_mean  = mean(eta_ref),
                    eta_sd    = max(stats::sd(eta_ref), EPSILON_DEFAULT),
                    signal_sd = FACE_BINARY_SIGNAL_SD)
      g_ref        <- face_binary_logit(eta_ref, calib)
      mu0_superpop <- mean(logistic(g_ref))
      mu1_superpop <- mean(logistic(g_ref + FACE_BINARY_ATE_TARGET))
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

#' Resolve per-site sample sizes for the FACE paper DGP
#'
#' Supports two mutually exclusive ways of specifying how observations are
#' distributed across the target and source sites:
#' \itemize{
#'   \item \strong{Explicit per-site} (federated-realistic): supply both
#'         \code{n_target} and \code{n_source_sizes}; the number of source sites
#'         \eqn{K} is taken from \code{length(n_source_sizes)} and the total is
#'         their sum.
#'   \item \strong{Total with equal split} (default, matches the FACE paper's
#'         equal \eqn{n_k}): supply \code{n_total}; each source receives
#'         \eqn{\lfloor n_{total}/(K+1)\rfloor} observations and the target the
#'         remainder.
#' }
#'
#' @param n_total        Total sample size across all sites (equal-split mode).
#' @param n_target       Target-site sample size (explicit per-site mode).
#' @param n_source_sizes Integer vector of per-site source sizes (explicit mode).
#' @param K              Number of source sites (equal-split mode only; ignored
#'                       in per-site mode, where it is taken from
#'                       \code{length(n_source_sizes)}).
#' @return Named list with \code{n_total}, \code{n_target},
#'   \code{n_source_sizes} (length-\eqn{K} integer vector), and \code{K}.
#' @keywords internal
resolve_face_site_sizes <- function(n_total = NULL, n_target = NULL,
                                    n_source_sizes = NULL, K = NULL) {
  per_site <- !is.null(n_source_sizes) || !is.null(n_target)

  if (per_site) {
    if (is.null(n_source_sizes) || is.null(n_target)) {
      stop("Per-site allocation requires BOTH n_target and n_source_sizes.",
           call. = FALSE)
    }
    n_target       <- as.integer(n_target)
    n_source_sizes <- as.integer(n_source_sizes)
    if (length(n_source_sizes) < 1L || anyNA(n_source_sizes) ||
        any(n_source_sizes <= 0L)) {
      stop("n_source_sizes must be a non-empty vector of positive integers.",
           call. = FALSE)
    }
    if (is.na(n_target) || n_target <= 0L) {
      stop("n_target must be a positive integer.", call. = FALSE)
    }
    # In per-site mode the number of source sites is defined by n_source_sizes;
    # any K argument is redundant and ignored (callers report the resolved K).
    K       <- length(n_source_sizes)
    n_total <- n_target + sum(n_source_sizes)
  } else {
    if (is.null(n_total) || is.null(K)) {
      stop("Equal-split allocation requires both n_total and K.", call. = FALSE)
    }
    n_source_per_site <- floor(n_total / (K + 1L))
    if (n_source_per_site <= 0L) {
      stop(sprintf(
        "n_total (%d) is too small for K=%d sites (each source would get 0).",
        as.integer(n_total), as.integer(K)), call. = FALSE)
    }
    n_target       <- n_total - K * n_source_per_site
    n_source_sizes <- rep(n_source_per_site, K)
  }

  list(n_total        = as.integer(n_total),
       n_target       = as.integer(n_target),
       n_source_sizes = as.integer(n_source_sizes),
       K              = as.integer(K))
}

#' Orchestrate full data generation for the FACE paper DGP
#'
#' Creates a dataset whose structure (field names and types) is identical to that
#' produced by the RoCE DGP path of \code{generate_simulation_data()}, so that
#' downstream code (\code{split_data_by_site()}, estimation algorithms, etc.)
#' requires no modification.
#'
#' \strong{Config semantics for the FACE paper DGP:}
#' \itemize{
#'   \item C1 (rich calibration, correct OR): calibration uses
#'         \eqn{[X, X^2]}; OR uses \eqn{[X - \kappa, X^2]}.
#'   \item C2 (rich calibration, misspecified OR): calibration uses
#'         \eqn{[X, X^2]};
#'         OR uses only \eqn{X} (drops quadratic block).
#'   \item C3 (misspecified calibration, correct OR): calibration uses only
#'         \eqn{X};
#'         OR uses \eqn{[X - \kappa, X^2]}.
#'   \item C4 (both working bases reduced): calibration uses only \eqn{X};
#'         OR uses only \eqn{X}.
#' }
#' The treatment propensity is quadratic before the DGP's finite-sample
#' clipping, but the skew-normal source-to-target density ratio is not exactly
#' log-quadratic. Thus C1--C2 provide a rich calibration basis rather than an
#' exactly specified merged site/treatment model.
#'
#' @param n_total          Total observations across all sites (equal-split
#'                         allocation). Ignored when \code{n_source_sizes} is
#'                         supplied; may be left \code{NULL} in that case.
#' @param K                Number of source sites (equal-split allocation).
#' @param p                Number of covariates (only first 4 have non-zero coefficients).
#' @param config           Configuration string: "C1", "C2", "C3", or "C4".
#' @param estimand_type    "superpopulation" (fixed truth) or "sample" (realized truth).
#' @param outcome_type     "continuous" (default) or "binary".
#' @param ate_deviation    ATE deviation for non-informative source sites (default 0).
#' @param n_deviated_sites Number of leading source sites with deviated ATE (default 0L).
#' @param effect_mod_strength Source-only effect-modification strength; zero
#'   recovers the standard FACE DGP.
#' @param n_target         Optional target-site sample size for explicit per-site
#'                         allocation (paired with \code{n_source_sizes}).
#' @param n_source_sizes   Optional integer vector of per-site source sample
#'                         sizes; when supplied, \code{K} and the total are taken
#'                         from it and \code{n_total} is ignored.
#' @return Named list with the same fields as the RoCE DGP output.
generate_face_data <- function(n_total = NULL, K = 3, p = 4, config = "C1",
                                     estimand_type    = "superpopulation",
                                     outcome_type     = "continuous",
                                     ate_deviation    = 0.0,
                                     n_deviated_sites = 0L,
                                     effect_mod_strength = 0,
                                     n_target         = NULL,
                                     n_source_sizes   = NULL) {
  # ---- 1. Resolve per-site sample sizes. Equal split across K + 1 sites by
  #         default (FACE-paper convention); explicit per-site sizes when both
  #         n_target and n_source_sizes are supplied. ----
  site_sizes     <- resolve_face_site_sizes(n_total, n_target, n_source_sizes, K)
  n_total        <- site_sizes$n_total
  K              <- site_sizes$K
  n_target       <- site_sizes$n_target
  n_source_sizes <- site_sizes$n_source_sizes

  # ---- 2. Generate site-specific covariates from skewed-normal distribution ----
  cov_list <- generate_face_covariates(
    n_target       = n_target,
    n_source_sizes = n_source_sizes,
    p              = p
  )
  X  <- cov_list$X
  R  <- cov_list$R
  X_dagger <- X   # No nonlinear transform in this DGP; keep for structural consistency

  # ---- 3. Build design matrices for PS and Outcome (config-dependent) ----
  X_sq <- X^2
  kappa <- cov_list$kappa

  # Z_site: used to fit the merged site / treatment calibration model.
  #   C1, C2: rich quadratic working basis
  #   C3, C4: reduced linear-only basis (drops the squared block)
  Z_site <- if (config %in% c("C1", "C2")) cbind(X, X_sq) else X

  # W_outcome: used to fit the outcome regression.
  #   C1, C3: correctly specified [X - kappa, X^2] basis
  #   C2, C4: misspecified X-only basis (drops the squared block)
  X_centered <- sweep(X, 2L, kappa, "-")
  W_outcome  <- if (config %in% c("C1", "C3")) cbind(X_centered, X_sq) else X

  # ---- 4. Propensity scores and treatment assignment ----
  ps_params <- get_face_ps_parameters(p)
  p_treat   <- calculate_face_propensity(X, ps_params$alpha1, ps_params$alpha2)
  A         <- rbinom(n_total, 1L, p_treat)

  # ---- 5. Site-specific treatment shifts and potential outcomes ----
  #   The treatment shift lives on each outcome's linear-predictor scale: a mean shift
  #   (Δ_T) for the continuous outcome, a log-odds shift (Δ_bin) for the binary
  #   outcome. The binary linear predictor is additionally standardized to a
  #   controlled logit-scale spread (binary_calib) to avoid expit() saturation.
  out_params   <- get_face_outcome_parameters(p)
  is_binary    <- outcome_type == "binary"
  base_ate     <- if (is_binary) FACE_BINARY_ATE_TARGET else FACE_ATE_TARGET
  binary_calib <- if (is_binary) get_face_binary_calibration(p, kappa = kappa) else NULL
  ate_map      <- build_face_ate_map(K, ate_deviation, n_deviated_sites,
                                     base_ate = base_ate)

  # Generate both potential outcomes first, then select the observed outcome.
  # This enforces consistency exactly: Y_i = A_i Y_i(1) + (1-A_i)Y_i(0).
  # The conditional observed-data distribution is unchanged relative to drawing
  # Y directly at the realized treatment, while avoiding a third outcome draw.
  A_ones  <- rep(1L, n_total)
  A_zeros <- rep(0L, n_total)
  Y_1 <- generate_face_outcomes(X, A_ones,  R, out_params$beta_linear,
                                      out_params$beta_squared, ate_map, kappa,
                                      outcome_type = outcome_type,
                                      binary_calib = binary_calib,
                                      effect_mod_strength = effect_mod_strength)
  Y_0 <- generate_face_outcomes(X, A_zeros, R, out_params$beta_linear,
                                      out_params$beta_squared, ate_map, kappa,
                                      outcome_type = outcome_type,
                                      binary_calib = binary_calib,
                                      effect_mod_strength = effect_mod_strength)
  Y <- ifelse(A == 1L, Y_1, Y_0)

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
      # Same deterministic logit transform as the data (binary_calib), so the
      # realized-sample estimand matches how Y was generated.
      g_target <- face_binary_logit(eta_target, binary_calib)
      mu0_true <- mean(logistic(g_target))
      mu1_true <- mean(logistic(g_target + base_ate))
    } else {
      mu0_true <- mean(eta_target)
      mu1_true <- mu0_true + base_ate
    }
    mu1_realized <- mean(Y_1[target_idx])
    mu0_realized <- mean(Y_0[target_idx])

  } else {
    # Superpopulation truth: fixed across simulations (Monte Carlo integration)
    superpop_truth <- calculate_face_truth(
      p = p,
      kappa = kappa,
      outcome_type = outcome_type,
      binary_calib = binary_calib
    )
    mu1_true     <- superpop_truth$mu1_superpop
    mu0_true     <- superpop_truth$mu0_superpop
    mu1_realized <- mean(Y_1[target_idx])   # sample-specific (for reference logging)
    mu0_realized <- mean(Y_0[target_idx])
  }

  # ---- 7. Assemble output list (matches RoCE DGP structure) ----
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
    # Historical `ate_map` name retained for compatibility. For binary outcomes
    # these are conditional log-odds shifts, not marginal risk-difference ATEs.
    ate_map      = ate_map,
    treatment_shift_map = ate_map,
    config       = config,
    Z_site       = Z_site,
    W_outcome    = W_outcome,
    estimand_type = estimand_type,
    outcome_type  = outcome_type,
    dgp_type      = "face"
  )
}
