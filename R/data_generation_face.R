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
calculate_face_propensity <- function(X, alpha1, alpha2, X_dagger = NULL,
                                      misspecification_strength = 0) {
  eta <- .face_mixed_predictor(X, X_dagger, alpha1, alpha2, kappa = 0,
                               strength = misspecification_strength)
  p_treat <- logistic(eta)
  pmax(pmin(p_treat, POSITIVITY_UPPER), POSITIVITY_LOWER)
}

# Kang--Schafer-style transforms of the leading signal coordinates. The true
# C2--C4 mechanisms mix these with the untransformed coordinates so that the
# truth lies outside the quadratic working basis phi(X) (main.tex,
# sec:simulations). With fewer than four covariates only the transforms whose
# inputs exist are used.
.face_transformed_coordinates <- function(X) {
  transforms <- list(
    function(X) exp(X[, 1] / 2),
    function(X) X[, 2] / (1 + exp(X[, 1])) + 10,
    function(X) (X[, 1] * X[, 3] / 25 + 0.6)^3,
    function(X) (X[, 2] + X[, 4] + 20)^2
  )
  n_signal <- min(FACE_SIGNAL_COORDINATES, ncol(X))
  do.call(cbind, lapply(transforms[seq_len(n_signal)], function(transform) transform(X)))
}

.face_standardization <- function(transformed) {
  list(mean = colMeans(transformed), sd = apply(transformed, 2L, stats::sd))
}

# X with its signal coordinates replaced by the transformed coordinates,
# standardized on the target reference population to location kappa and unit
# scale so that the same coefficient vectors apply to X and X_dagger.
.face_x_dagger <- function(X, standardization, kappa) {
  transformed <- .face_transformed_coordinates(X)
  for (j in seq_len(ncol(transformed))) {
    X[, j] <- kappa + (transformed[, j] - standardization$mean[j]) / standardization$sd[j]
  }
  X
}

# Linear-plus-quadratic predictor (X - kappa)' linear + (X^2)' squared, mixed
# between X and X_dagger: eta_omega = (1 - omega) eta(X) + omega eta(X_dagger).
.face_mixed_predictor <- function(X, X_dagger, linear, squared, kappa, strength) {
  predictor <- function(M) {
    as.numeric(sweep(M, 2L, kappa, "-") %*% linear + (M^2) %*% squared)
  }
  if (strength == 0) return(predictor(X))
  if (is.null(X_dagger)) {
    stop(".face_mixed_predictor: X_dagger is required when the misspecification strength is positive.",
         call. = FALSE)
  }
  (1 - strength) * predictor(X) + strength * predictor(X_dagger)
}

# Which true mechanisms a configuration misspecifies.
.face_misspecification_strengths <- function(config, strength) {
  if (length(strength) != 1L || !is.finite(strength) || strength < 0 || strength > 1) {
    stop("misspecification_strength must be a single number in [0, 1].", call. = FALSE)
  }
  list(
    outcome = if (config %in% c("C2", "C4")) strength else 0,
    propensity = if (config %in% c("C3", "C4")) strength else 0
  )
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

#' Reference population of the FACE paper DGP
#'
#' One deterministic draw of the target covariate law (\code{n_ref} units,
#' \code{ref_seed}) supplies the standardization of the transformed
#' coordinates, the standardized-logit calibration of the binary outcome, and
#' the superpopulation truth for the requested configuration. Both the data
#' generator and the truth use this same population, so truth and data share
#' one transform by construction.
#'
#' @param p Number of covariates.
#' @param kappa Location parameter of the covariate law.
#' @param config Configuration string (see \code{generate_face_data()}).
#' @param misspecification_strength Mixing weight omega of the transformed
#'   coordinates in the misspecified mechanisms.
#' @param outcome_type "continuous" or "binary".
#' @param n_ref,ref_seed Size and seed of the reference draw.
#' @return List with \code{standardization}, \code{calibration} (binary only:
#'   \code{eta_mean}, \code{eta_sd}, \code{signal_sd}), \code{mu1_superpop},
#'   \code{mu0_superpop}, \code{ate_superpop}, and the arguments.
#' @keywords internal
.face_reference_population <- function(p, kappa, config, misspecification_strength,
                                       outcome_type, n_ref = 100000L,
                                       ref_seed = 99999L) {
  strengths <- .face_misspecification_strengths(config, misspecification_strength)
  out_params <- get_face_outcome_parameters(p)
  with_seed(ref_seed, {
    X_ref <- generate_face_covariates(n_target       = n_ref,
                                      n_source_sizes = integer(0),
                                      p              = p,
                                      kappa          = kappa)$X
    standardization <- .face_standardization(.face_transformed_coordinates(X_ref))
    X_dagger_ref <- .face_x_dagger(X_ref, standardization, kappa)
    eta_ref <- .face_mixed_predictor(
      X_ref, X_dagger_ref, out_params$beta_linear, out_params$beta_squared,
      kappa, strengths$outcome
    )
    if (outcome_type == "binary") {
      calibration <- list(
        eta_mean = mean(eta_ref),
        eta_sd = max(stats::sd(eta_ref), EPSILON_DEFAULT),
        signal_sd = FACE_BINARY_SIGNAL_SD
      )
      g_ref <- face_binary_logit(eta_ref, calibration)
      mu0_superpop <- mean(logistic(g_ref))
      mu1_superpop <- mean(logistic(g_ref + FACE_BINARY_ATE_TARGET))
    } else {
      calibration <- NULL
      mu0_superpop <- mean(eta_ref)
      mu1_superpop <- mu0_superpop + FACE_ATE_TARGET
    }
    list(
      standardization = standardization,
      calibration = calibration,
      mu1_superpop = mu1_superpop,
      mu0_superpop = mu0_superpop,
      ate_superpop = mu1_superpop - mu0_superpop,
      n_ref = as.integer(n_ref),
      ref_seed = as.integer(ref_seed),
      p = as.integer(p),
      kappa = as.numeric(kappa),
      config = config,
      misspecification_strength = as.numeric(misspecification_strength)
    )
  })
}

#' Binary-outcome calibration of the FACE paper DGP
#'
#' The standardized-logit calibration of \code{.face_reference_population()}
#' together with the reference-population outcome means. The configuration is
#' required because C2 and C4 change the outcome mechanism and hence the
#' calibration.
#'
#' @inheritParams .face_reference_population
#' @param signal_sd Target logit-scale standard deviation of the signal.
#' @keywords internal
get_face_binary_calibration <- function(p, config, kappa = FACE_KAPPA,
                                        signal_sd = FACE_BINARY_SIGNAL_SD,
                                        n_ref = 100000L, ref_seed = 99999L,
                                        misspecification_strength = FACE_MISSPECIFICATION_STRENGTH) {
  if (!identical(as.numeric(signal_sd), as.numeric(FACE_BINARY_SIGNAL_SD))) {
    stop("get_face_binary_calibration: signal_sd is fixed at FACE_BINARY_SIGNAL_SD.",
         call. = FALSE)
  }
  reference <- .face_reference_population(
    p, kappa, config, misspecification_strength, "binary", n_ref, ref_seed
  )
  c(reference$calibration,
    reference[c("mu0_superpop", "mu1_superpop", "n_ref", "ref_seed", "p", "kappa")])
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
#' @param X_dagger Transformed covariates from \code{.face_x_dagger()};
#'   required when \code{misspecification_strength > 0}.
#' @param misspecification_strength Mixing weight of the transformed
#'   coordinates in the true outcome predictor (0 recovers the quadratic
#'   mechanism).
#' @return Numeric vector of outcomes (length n).
generate_face_outcomes <- function(X, A, R, beta_lin, beta_sq, ate_map,
                                         kappa        = FACE_KAPPA,
                                         noise_sd     = FACE_NOISE_SD,
                                         outcome_type = "continuous",
                                         binary_calib = NULL,
                                         effect_mod_strength = 0,
                                         em_direction = NULL,
                                         X_dagger = NULL,
                                         misspecification_strength = 0) {
  eta        <- .face_mixed_predictor(X, X_dagger, beta_lin, beta_sq, kappa,
                                      misspecification_strength)
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
           as.numeric(sweep(X, 2L, kappa, "-") %*% em_direction)
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

#' Superpopulation truth of the FACE paper DGP (target-population specific)
#'
#' Target-population potential-outcome means for a configuration, computed on
#' the deterministic reference population of \code{.face_reference_population()}.
#'
#' @inheritParams .face_reference_population
#' @return Named list with \code{mu1_superpop}, \code{mu0_superpop},
#'   \code{ate_superpop}, \code{n_ref}, and \code{ref_seed}.
calculate_face_truth <- function(p, config, kappa = FACE_KAPPA,
                                 outcome_type = "continuous",
                                 misspecification_strength = FACE_MISSPECIFICATION_STRENGTH,
                                 n_ref = 100000L, ref_seed = 99999L) {
  reference <- .face_reference_population(
    p, kappa, config, misspecification_strength, outcome_type, n_ref, ref_seed
  )
  reference[c("mu1_superpop", "mu0_superpop", "ate_superpop", "n_ref", "ref_seed")]
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
#' \strong{Config semantics for the FACE paper DGP:} both nuisance models
#' always use the working basis \eqn{\phi(X) = [X - \kappa, X^2]}
#' (\code{Z_site} and \code{W_outcome} are identical). Misspecification lives
#' in the true mechanism: the leading signal coordinates are replaced by
#' standardized Kang--Schafer-style transforms \eqn{X^\dagger}
#' (\code{.face_x_dagger()}) and the true predictor becomes
#' \eqn{(1 - \omega)\eta(X) + \omega\eta(X^\dagger)} with
#' \eqn{\omega = } \code{misspecification_strength}.
#' \itemize{
#'   \item C1: outcome and treatment mechanisms use \eqn{\eta(X)} (both
#'         working models correctly specified).
#'   \item C2: outcome mechanism misspecified (\eqn{\eta_\omega}), treatment
#'         mechanism correct.
#'   \item C3: treatment mechanism misspecified, outcome mechanism correct.
#'   \item C4: both mechanisms misspecified.
#' }
#' The treatment propensity is quadratic before the DGP's finite-sample
#' clipping, but the skew-normal source-to-target density ratio is not exactly
#' log-quadratic, so even C1 is a rich working-model setting rather than an
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
#' @param misspecification_strength Mixing weight \eqn{\omega} of the
#'   transformed coordinates in the misspecified mechanisms of C2--C4
#'   (default \code{FACE_MISSPECIFICATION_STRENGTH}); ignored by C1. The
#'   misspecified configurations are validated for binary outcomes only
#'   (HISTORY #0002/#0006); the continuous outcome under C2/C4 inherits the
#'   heavy-tailed transformed predictor without calibration.
#' @return Named list with the same fields as the RoCE DGP output;
#'   \code{X_dagger} holds the transformed covariates and
#'   \code{misspecification_strength} the strength actually applied (0 under
#'   C1).
generate_face_data <- function(n_total = NULL, K = 3, p = 4, config = "C1",
                                     estimand_type    = "superpopulation",
                                     outcome_type     = "continuous",
                                     ate_deviation    = 0.0,
                                     n_deviated_sites = 0L,
                                     effect_mod_strength = 0,
                                     n_target         = NULL,
                                     n_source_sizes   = NULL,
                                     misspecification_strength = FACE_MISSPECIFICATION_STRENGTH) {
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
  kappa <- cov_list$kappa

  # ---- 3. Reference population, transformed covariates, common working basis ----
  # The reference draw fixes the standardization of the transformed coordinates,
  # the binary calibration and the superpopulation truth; it does not advance
  # the caller's RNG stream.
  strengths <- .face_misspecification_strengths(config, misspecification_strength)
  reference <- .face_reference_population(
    p, kappa, config, misspecification_strength, outcome_type
  )
  X_dagger <- .face_x_dagger(X, reference$standardization, kappa)
  basis <- cbind(sweep(X, 2L, kappa, "-"), X^2)

  # ---- 4. Propensity scores and treatment assignment ----
  ps_params <- get_face_ps_parameters(p)
  p_treat   <- calculate_face_propensity(X, ps_params$alpha1, ps_params$alpha2,
                                         X_dagger, strengths$propensity)
  A         <- rbinom(n_total, 1L, p_treat)

  # ---- 5. Site-specific treatment shifts and potential outcomes ----
  #   The treatment shift lives on each outcome's linear-predictor scale: a mean shift
  #   (Δ_T) for the continuous outcome, a log-odds shift (Δ_bin) for the binary
  #   outcome. The binary linear predictor is additionally standardized to a
  #   controlled logit-scale spread (calibration) to avoid expit() saturation.
  out_params   <- get_face_outcome_parameters(p)
  is_binary    <- outcome_type == "binary"
  base_ate     <- if (is_binary) FACE_BINARY_ATE_TARGET else FACE_ATE_TARGET
  binary_calib <- reference$calibration
  ate_map      <- build_face_ate_map(K, ate_deviation, n_deviated_sites,
                                     base_ate = base_ate)

  # Generate both potential outcomes first, then select the observed outcome.
  # This enforces consistency exactly: Y_i = A_i Y_i(1) + (1-A_i)Y_i(0).
  # The conditional observed-data distribution is unchanged relative to drawing
  # Y directly at the realized treatment, while avoiding a third outcome draw.
  draw_outcomes <- function(A_fixed) {
    generate_face_outcomes(X, A_fixed, R, out_params$beta_linear,
                           out_params$beta_squared, ate_map, kappa,
                           outcome_type = outcome_type,
                           binary_calib = binary_calib,
                           effect_mod_strength = effect_mod_strength,
                           X_dagger = X_dagger,
                           misspecification_strength = strengths$outcome)
  }
  Y_1 <- draw_outcomes(rep(1L, n_total))
  Y_0 <- draw_outcomes(rep(0L, n_total))
  Y <- ifelse(A == 1L, Y_1, Y_0)

  # ---- 6. True potential-outcome means ----
  target_idx <- R == "t"

  if (estimand_type == "sample") {
    # Sample-specific truth: computed from the realized target observations
    eta_target <- .face_mixed_predictor(
      X[target_idx, , drop = FALSE], X_dagger[target_idx, , drop = FALSE],
      out_params$beta_linear, out_params$beta_squared, kappa, strengths$outcome
    )
    if (is_binary) {
      # Same deterministic logit transform as the data (binary_calib), so the
      # realized-sample estimand matches how Y was generated.
      g_target <- face_binary_logit(eta_target, binary_calib)
      mu0_true <- mean(logistic(g_target))
      mu1_true <- mean(logistic(g_target + base_ate))
    } else {
      mu0_true <- mean(eta_target)
      mu1_true <- mu0_true + base_ate
    }
  } else {
    # Superpopulation truth: fixed across simulations (Monte Carlo integration)
    mu1_true <- reference$mu1_superpop
    mu0_true <- reference$mu0_superpop
  }
  mu1_realized <- mean(Y_1[target_idx])   # sample-specific (for reference logging)
  mu0_realized <- mean(Y_0[target_idx])

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
    misspecification_strength = if (any(unlist(strengths) > 0)) misspecification_strength else 0,
    Z_site       = basis,
    W_outcome    = basis,
    estimand_type = estimand_type,
    outcome_type  = outcome_type,
    dgp_type      = "face"
  )
}
