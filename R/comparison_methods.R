# comparison_methods.R - Multi-site comparison estimators
#
# =============================================================================
# TARGET ESTIMAND
# =============================================================================
# Arm-specific entry points estimate a target potential-outcome mean,
# μᵃ_t = E_t[Y(a)]. run_all_comparisons_tate() pairs the two arm-specific
# influence blocks before estimating the target average treatment effect,
# τ_t = μ¹_t - μ⁰_t, so the cross-arm covariance is retained.
# =============================================================================
#
# This file contains multi-site estimators that aggregate across sites:
#   1. Sample-size weighted (SS)
#   2. Inverse-variance weighted (IVW) with DerSimonian-Laird heterogeneity
#   3. Tilted AIPW (exponential tilting with unpenalized MLE nuisances)
#   4. Federated DR-AIPW (density-ratio corrected, federated)
#   5. Pooled DR-AIPW (density-ratio corrected, centralized)
#   6. run_all_comparisons (orchestrator)
#
# Helper functions are in estimators_helpers.R.
# Target-only estimators are in estimators_target.R.
# Oracle estimator is in estimators_oracle.R.

.name_influence_blocks <- function(blocks, block_names) {
  if (length(blocks) != length(block_names)) {
    stop("influence block names and blocks must have the same length.",
         call. = FALSE)
  }
  stats::setNames(blocks, block_names)
}

.comparison_influence_variance <- function(blocks) {
  if (length(blocks) == 0L) {
    stop("at least one influence block is required.", call. = FALSE)
  }
  sum(vapply(blocks, function(block) {
    influence <- as.numeric(block$influence)
    weight <- as.numeric(block$weight)
    if (length(weight) != 1L || !is.finite(weight) ||
        length(influence) == 0L || any(!is.finite(influence))) {
      stop("comparison influence blocks must contain finite influence vectors and scalar weights.",
           call. = FALSE)
    }
    centered <- influence - mean(influence)
    weight^2 * mean(centered^2) / length(centered)
  }, numeric(1L)))
}

.density_ratio_clipping_diagnostics <- function(weights_by_site) {
  empty <- list(
    dr_weight_n = NA_integer_,
    dr_weight_n_clipped = NA_integer_,
    dr_weight_fraction_clipped = NA_real_,
    dr_weight_max_site_fraction_clipped = NA_real_,
    dr_weight_min_before_clipping = NA_real_,
    dr_weight_max_before_clipping = NA_real_
  )
  if (is.null(weights_by_site) || length(weights_by_site) == 0L) {
    return(empty)
  }
  diagnostics <- lapply(
    weights_by_site,
    function(weights) attr(weights, "clipping_diagnostics", exact = TRUE)
  )
  if (any(vapply(diagnostics, is.null, logical(1L)))) {
    return(empty)
  }
  required <- c(
    "n", "n_clipped", "fraction_clipped", "preclip_min", "preclip_max"
  )
  invalid <- vapply(diagnostics, function(diagnostic) {
    if (!is.list(diagnostic) ||
        length(setdiff(required, names(diagnostic))) > 0L) {
      return(TRUE)
    }
    values <- stats::setNames(suppressWarnings(as.numeric(unlist(
      diagnostic[required], use.names = FALSE
    ))), required)
    if (length(values) != length(required) || any(!is.finite(values))) {
      return(TRUE)
    }
    n <- values[["n"]]
    n_clipped <- values[["n_clipped"]]
    fraction <- values[["fraction_clipped"]]
    n < 1 || n != floor(n) || n_clipped < 0 ||
      n_clipped != floor(n_clipped) || n_clipped > n ||
      fraction < 0 || fraction > 1 ||
      abs(fraction - n_clipped / n) > 1e-12 ||
      values[["preclip_min"]] <= 0 ||
      values[["preclip_max"]] < values[["preclip_min"]]
  }, logical(1L))
  if (any(invalid)) {
    stop(
      "density-ratio clipping metadata are incomplete or inconsistent.",
      call. = FALSE
    )
  }
  n <- sum(vapply(diagnostics, `[[`, numeric(1L), "n"))
  n_clipped <- sum(vapply(diagnostics, `[[`, numeric(1L), "n_clipped"))
  list(
    dr_weight_n = as.integer(n),
    dr_weight_n_clipped = as.integer(n_clipped),
    dr_weight_fraction_clipped = n_clipped / n,
    dr_weight_max_site_fraction_clipped = max(vapply(
      diagnostics, `[[`, numeric(1L), "fraction_clipped"
    )),
    dr_weight_min_before_clipping = min(vapply(
      diagnostics, `[[`, numeric(1L), "preclip_min"
    )),
    dr_weight_max_before_clipping = max(vapply(
      diagnostics, `[[`, numeric(1L), "preclip_max"
    ))
  )
}

.combine_comparison_arms <- function(
    mu1_result, mu0_result, method,
    variance_method = c("bootstrap", "analytic"),
    n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  variance_method <- match.arg(variance_method)
  blocks1 <- mu1_result$influence_blocks
  blocks0 <- mu0_result$influence_blocks
  if (is.null(blocks1) || is.null(blocks0)) {
    stop(sprintf(
      "%s does not expose both arm-specific influence blocks.", method
    ), call. = FALSE)
  }
  if (is.null(names(blocks1)) || is.null(names(blocks0)) ||
      !setequal(names(blocks1), names(blocks0))) {
    stop(sprintf(
      "%s arm-specific influence blocks have inconsistent site names.", method
    ), call. = FALSE)
  }

  blocks0 <- blocks0[names(blocks1)]
  tate_blocks <- Map(function(block1, block0) {
    influence1 <- as.numeric(block1$influence)
    influence0 <- as.numeric(block0$influence)
    if (length(influence1) != length(influence0)) {
      stop(sprintf(
        "%s arm-specific influence blocks have inconsistent lengths.", method
      ), call. = FALSE)
    }
    list(
      influence = as.numeric(block1$weight) * influence1 -
        as.numeric(block0$weight) * influence0,
      weight = 1
    )
  }, blocks1, blocks0)
  names(tate_blocks) <- names(blocks1)

  mu1_influence_variance <- .comparison_influence_variance(blocks1)
  mu0_influence_variance <- .comparison_influence_variance(blocks0)
  analytic_variance <- max(
    .comparison_influence_variance(tate_blocks), VARIANCE_MIN
  )
  cross_arm_covariance <- (
    mu1_influence_variance + mu0_influence_variance - analytic_variance
  ) / 2
  cross_arm_correlation <- cross_arm_covariance / sqrt(
    max(mu1_influence_variance, VARIANCE_MIN) *
      max(mu0_influence_variance, VARIANCE_MIN)
  )
  variance_result <- .resolve_comparison_variance(
    analytic_variance = analytic_variance,
    blocks = tate_blocks,
    variance_method = variance_method,
    n_bootstrap = n_bootstrap
  )
  estimate <- as.numeric(mu1_result$estimate - mu0_result$estimate)
  dr_diagnostics1 <- mu1_result$components$dr_weight_diagnostics
  dr_diagnostics0 <- mu0_result$components$dr_weight_diagnostics
  if (!is.null(dr_diagnostics1) && !is.null(dr_diagnostics0) &&
      !identical(dr_diagnostics1, dr_diagnostics0)) {
    stop(sprintf(
      "%s treatment arms carry inconsistent density-ratio diagnostics.",
      method
    ), call. = FALSE)
  }
  dr_diagnostics <- dr_diagnostics1 %||% dr_diagnostics0 %||% NULL

  list(
    estimate = estimate,
    variance = variance_result$variance,
    se = variance_result$se,
    method = method,
    n = max(mu1_result$n %||% 0L, mu0_result$n %||% 0L),
    influence_blocks = tate_blocks,
    components = list(
      estimand = "TATE",
      variance_method = variance_result$variance_method,
      variance_analytic = variance_result$variance_analytic,
      se_analytic = variance_result$se_analytic,
      se_bootstrap = variance_result$se_bootstrap,
      n_bootstrap = variance_result$n_bootstrap,
      mu1_estimate = mu1_result$estimate,
      mu0_estimate = mu0_result$estimate,
      mu1_reported_se = mu1_result$se,
      mu0_reported_se = mu0_result$se,
      mu1_influence_se = sqrt(mu1_influence_variance),
      mu0_influence_se = sqrt(mu0_influence_variance),
      cross_arm_covariance = cross_arm_covariance,
      cross_arm_correlation = cross_arm_correlation,
      dr_weight_diagnostics = dr_diagnostics
    )
  )
}


#' Run deterministic independent fits across named sites
#'
#' @param sites Character site identifiers.
#' @param fit_function Function accepting one site identifier.
#' @param n_cores Number of independent fits to run in parallel.
#' @return Named list of site-level results.
#' @keywords internal
.parallel_site_fits <- function(sites, fit_function, n_cores = 1L) {
  sites <- as.character(sites)
  if (length(sites) == 0L) {
    return(stats::setNames(list(), character(0)))
  }
  actual_cores <- min(setup_parallel(n_cores), length(sites))
  # Give every site an isolated seed so fold construction is identical under
  # sequential and forked execution. Restore the parent to the state just after
  # seed allocation; fork bookkeeping must not perturb downstream bootstraps.
  site_seeds <- stats::setNames(
    sample.int(.Machine$integer.max, length(sites)),
    sites
  )
  rng_state_after_seed_allocation <- get(
    ".Random.seed", envir = globalenv(), inherits = FALSE
  )
  on.exit(
    assign(
      ".Random.seed", rng_state_after_seed_allocation,
      envir = globalenv()
    ),
    add = TRUE
  )
  site_fits <- parallel_lapply(sites, function(site) {
    with_seed(site_seeds[[site]], fit_function(site))
  }, n_cores = actual_cores)
  stats::setNames(site_fits, sites)
}

#' Compute site-level AIPW fits once for naive comparison methods
#'
#' Internal helper used by sample-size and inverse-variance baselines to avoid
#' duplicated nuisance-model fitting work.
#'
#' @param data_split Split data by site.
#' @param family GLM family.
#' @param use_rcal Logical.
#' @param use_crossfit Logical.
#' @param n_folds Optional number of folds.
#' @param A_val Treatment value.
#' @param n_cores Number of site-level fits to run in parallel.
#' @return Named list of per-site \code{fit_site_aipw} outputs.
#' @keywords internal
.fit_site_aipw_all_sites <- function(data_split, family = "binomial",
                                     use_rcal = FALSE, use_crossfit = TRUE,
                                     n_folds = NULL, A_val = 1L,
                                     n_cores = 1L) {
  sites <- names(data_split)
  .parallel_site_fits(sites, function(site) {
    fit_site_aipw(
      data_split[[site]],
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }, n_cores = n_cores)
}


#' Fit source-to-target density-ratio weights once for DR baselines
#'
#' Federated-DR and Pooled-DR use the same density-ratio weights, and those
#' weights do not depend on the treatment arm. This internal helper avoids
#' repeating the same cross-validation up to four times per source in a TATE
#' comparison.
#'
#' @param data_split Split data by site.
#' @param dr_lambda Optional fixed density-ratio penalty.
#' @param dr_lambda_rule Density-ratio CV rule.
#' @param n_cores Number of source fits to run in parallel.
#' @param dr_weights_by_site Optional precomputed named list.
#' @param caller Calling function name used in validation errors.
#' @return Named list of normalized density-ratio weight vectors.
.resolve_dr_weights_by_site <- function(
    data_split, dr_lambda = NULL, dr_lambda_rule = c("min", "1se"),
    n_cores = 1L, dr_weights_by_site = NULL,
    caller = ".resolve_dr_weights_by_site") {
  dr_lambda_rule <- .match_nuisance_lambda_rule(
    dr_lambda_rule, caller, arg = "dr_lambda_rule"
  )
  source_sites <- setdiff(names(data_split), "t")

  if (is.null(dr_weights_by_site)) {
    if (length(source_sites) == 0L) {
      return(stats::setNames(list(), character(0)))
    }
    target_covariates <- as.matrix(data_split[["t"]]$Z_site)
    actual_cores <- min(setup_parallel(n_cores), length(source_sites))
    # Density-ratio CV may use random fold assignment. Give every source a
    # stable isolated seed and restore the caller's RNG exactly, so caching and
    # source-level parallelism cannot perturb subsequent outcome-model CV.
    had_rng_state <- exists(
      ".Random.seed", envir = globalenv(), inherits = FALSE
    )
    if (had_rng_state) {
      rng_state <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    }
    on.exit({
      if (had_rng_state) {
        assign(".Random.seed", rng_state, envir = globalenv())
      } else if (exists(
        ".Random.seed", envir = globalenv(), inherits = FALSE
      )) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
    source_seeds <- stats::setNames(
      104729L + seq_along(source_sites),
      source_sites
    )
    fitted_weights <- parallel_lapply(source_sites, function(site) {
      with_seed(source_seeds[[site]], {
        calculate_dr_weights(
          Z_source = as.matrix(data_split[[site]]$Z_site),
          Z_target = target_covariates,
          lambda = dr_lambda,
          lambda_rule = dr_lambda_rule
        )
      })
    }, n_cores = actual_cores)
    return(stats::setNames(fitted_weights, source_sites))
  }

  if (!is.list(dr_weights_by_site) ||
      is.null(names(dr_weights_by_site)) ||
      anyDuplicated(names(dr_weights_by_site))) {
    stop(
      caller,
      ": dr_weights_by_site must be a uniquely named list.",
      call. = FALSE
    )
  }
  missing_sites <- setdiff(source_sites, names(dr_weights_by_site))
  extra_sites <- setdiff(names(dr_weights_by_site), source_sites)
  if (length(missing_sites) > 0L || length(extra_sites) > 0L) {
    stop(
      sprintf(
        "%s: dr_weights_by_site names must match source sites; missing={%s}, extra={%s}.",
        caller,
        paste(missing_sites, collapse = ","),
        paste(extra_sites, collapse = ",")
      ),
      call. = FALSE
    )
  }
  for (site in source_sites) {
    weights <- dr_weights_by_site[[site]]
    expected_length <- data_split[[site]]$n
    if (!is.numeric(weights) || length(weights) != expected_length ||
        any(!is.finite(weights)) || any(weights <= 0)) {
      stop(
        sprintf(
          "%s: density-ratio weights for source '%s' must contain %d positive finite values.",
          caller, site, expected_length
        ),
        call. = FALSE
      )
    }
  }
  dr_weights_by_site[source_sites]
}


#' Fit reusable site-level components for the DR comparison estimators
#'
#' Federated-DR and Pooled-DR share the same source-level weighted AIPW fit and
#' density-ratio influence correction. Computing those components once avoids
#' duplicate high-dimensional nuisance CV and matrix solves.
#'
#' @param data_split Split data by site.
#' @param dr_weights_by_site Validated source density-ratio weights.
#' @param family Outcome GLM family.
#' @param A_val Treatment arm.
#' @param required_sites Sites needed by the caller.
#' @param n_cores Number of independent site fits to run in parallel.
#' @param dr_site_components Optional precomputed named list.
#' @param caller Calling function name for validation errors.
#' @return Named list of reusable target/source DR components.
#' @keywords internal
.resolve_dr_site_components <- function(
    data_split, dr_weights_by_site, family = "binomial", A_val = 1L,
    required_sites = names(data_split), n_cores = 1L,
    dr_site_components = NULL, caller = ".resolve_dr_site_components") {
  required_sites <- unique(as.character(required_sites))
  unknown_sites <- setdiff(required_sites, names(data_split))
  if (length(unknown_sites) > 0L) {
    stop(
      sprintf("%s: unknown required site(s): %s.", caller,
              paste(unknown_sites, collapse = ", ")),
      call. = FALSE
    )
  }

  if (is.null(dr_site_components)) {
    target_data <- data_split[["t"]]
    Z_target <- as.matrix(target_data$Z_site)
    mean_phi_target <- c(1, colMeans(Z_target))
    Z_target_centered <- sweep(
      cbind(1, Z_target), 2, mean_phi_target, "-"
    )

    return(.parallel_site_fits(required_sites, function(site) {
      site_data <- data_split[[site]]
      is_target <- identical(site, "t")
      dr_weights <- if (is_target) NULL else dr_weights_by_site[[site]]
      base_result <- calculate_weighted_site_aipw(
        y = site_data$Y,
        a = site_data$A,
        X = as.matrix(site_data$W_outcome),
        weights = dr_weights,
        family = family,
        A_val = A_val
      )

      if (is_target) {
        return(list(
          estimate = base_result$estimate,
          variance = base_result$variance,
          n = site_data$n,
          varphi_ot = base_result$varphi_ot,
          base_result = base_result
        ))
      }

      Z_source <- as.matrix(site_data$Z_site)
      Z_centered <- sweep(cbind(1, Z_source), 2, mean_phi_target, "-")
      n_source <- site_data$n
      d_alpha <- max(mean(dr_weights), DIVISION_FLOOR)
      phi_centered <- base_result$phi - base_result$estimate
      A_s <- -1 / d_alpha * colMeans(
        as.numeric(dr_weights) * Z_centered * as.numeric(phi_centered)
      )
      # calculate_dr_weights() uses A_dummy = 1 for every source observation.
      # Its score and Jacobian must therefore use the full source sample. The
      # AIPW pseudo-outcome already carries the requested-arm indicator in its
      # residual term.
      density_score_weights <- as.numeric(dr_weights)
      M_alpha <- t(Z_centered) %*%
        (Z_centered * density_score_weights) /
        max(1, n_source)
      adjustment <- solve_with_ridge(M_alpha) %*% A_s
      source_correction <- as.numeric(
        density_score_weights * (Z_centered %*% adjustment)
      )
      target_correction <- -as.numeric(
        Z_target_centered %*% adjustment
      )
      source_influence <-
        as.numeric(base_result$varphi_ot) + source_correction
      source_influence <- source_influence - mean(source_influence)

      list(
        estimate = base_result$estimate,
        variance = mean(source_influence^2) / max(1, n_source),
        n = n_source,
        varphi_ot = source_influence,
        target_if_component = target_correction,
        source_if_correction = source_correction,
        lambda = as.numeric(
          attr(dr_weights, "lambda_used") %||% NA_real_
        ),
        base_result = base_result
      )
    }, n_cores = n_cores))
  }

  if (!is.list(dr_site_components) ||
      is.null(names(dr_site_components)) ||
      anyDuplicated(names(dr_site_components))) {
    stop(
      caller, ": dr_site_components must be a uniquely named list.",
      call. = FALSE
    )
  }
  missing_sites <- setdiff(required_sites, names(dr_site_components))
  if (length(missing_sites) > 0L) {
    stop(
      sprintf("%s: dr_site_components is missing site(s): %s.", caller,
              paste(missing_sites, collapse = ", ")),
      call. = FALSE
    )
  }
  required_fields <- c("estimate", "variance", "n", "varphi_ot")
  for (site in required_sites) {
    component <- dr_site_components[[site]]
    missing_fields <- required_fields[
      vapply(required_fields, function(field) is.null(component[[field]]), logical(1L))
    ]
    if (length(missing_fields) > 0L) {
      stop(
        sprintf("%s: component for site '%s' is missing field(s): %s.",
                caller, site, paste(missing_fields, collapse = ", ")),
        call. = FALSE
      )
    }
  }
  dr_site_components[required_sites]
}


#' Sample-size adjusted estimator (SS) with correct variance estimation
#' 
#' This estimator computes a sample-size weighted average of site-specific
#' AIPW estimates. By default (\code{variance_method = "bootstrap"}) the standard
#' error is the multiplier (wild) bootstrap, i.e. the fixed-effects sampling variance
#' of the weighted average; the DerSimonian-Laird random-effects (\eqn{\tau^2}) variance
#' is available via \code{variance_method = "analytic"} and is retained in
#' \code{components} for reference.
#' 
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   for variance-valid inference (use_rcal is ignored in this mode).
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param A_val Treatment arm, either 0 or 1.
#' @param site_fits Optional precomputed per-site outputs from
#'   \code{fit_site_aipw}. When provided, avoids refitting nuisances.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws when
#'   \code{variance_method = "bootstrap"}.
#' @return estimate with variance
#' @export
estimate_sample_size_weighted <- function(data_split, family = "binomial",
                                          use_rcal = FALSE, use_crossfit = TRUE,
                                          n_folds = NULL, A_val = 1L,
                                          site_fits = NULL,
                                          variance_method = c("bootstrap", "analytic"),
                                          n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  sites <- names(data_split)
  site_estimates <- list()
  total_n <- 0

  if (is.null(site_fits)) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }
  
  for (site in sites) {
    res <- site_fits[[site]]
    site_estimates[[site]] <- list(
      estimate = res$estimate,
      variance = res$variance,
      n = res$n,
      V_ot = res$V_ot
    )
    total_n <- total_n + res$n
  }
  
  # ==========================================================================
  # VARIANCE ESTIMATION WITH HETEROGENEITY CORRECTION
  # ==========================================================================
  # Similar to IVW, we use a random-effects model to account for between-site
  # heterogeneity. When sites have different covariate distributions (covariate
  # shift), they may estimate different population parameters, introducing
  # extra variance that the naive pooled-IF formula doesn't capture.
  #
  # Random-effects model:
  #   μ̂_j ~ N(μ, Var_j + τ²)
  # where τ² is estimated via DerSimonian-Laird using Cochran's Q.
  #
  # The random-effects variance for sample-size weighted estimator is:
  #   Var_RE = Σ_j w_j² × (Var_j + τ²)
  # ==========================================================================
  
  K <- length(site_estimates)
  
  # Compute sample-size weighted estimate
  weighted_sum <- 0
  for (site in names(site_estimates)) {
    w <- site_estimates[[site]]$n / total_n
    weighted_sum <- weighted_sum + w * site_estimates[[site]]$estimate
  }
  
  # ==========================================================================
  # Estimate between-site heterogeneity using DerSimonian-Laird
  # ==========================================================================
  site_est_vec <- sapply(site_estimates, function(x) x$estimate)
  site_variances <- sapply(site_estimates, function(x) x$V_ot / x$n)
  het <- calculate_dl_heterogeneity(site_est_vec, site_variances)
  tau_sq <- het$tau_sq
  
  # ==========================================================================
  # Compute variance with heterogeneity correction
  # For sample-size weights: Var_RE = Σ_j (n_j/N)² × (Var_j + τ²)
  # ==========================================================================
  
  var_fe <- 0   # Fixed-effects variance (naive)
  var_re <- 0   # Random-effects variance (with τ²)
  
  for (site in names(site_estimates)) {
    n_j <- site_estimates[[site]]$n
    w_j <- n_j / total_n
    var_j <- site_variances[site]
    
    var_fe <- var_fe + w_j^2 * var_j
    var_re <- var_re + w_j^2 * (var_j + tau_sq)
  }
  
  # Variance: bootstrap (default) over per-site influence blocks with sample-size
  # weights, or the analytic random-effects (tau^2) variance. The bootstrap reproduces
  # the fixed-effects sampling variance and avoids the tau^2 inflation.
  boot_blocks <- .name_influence_blocks(
    lapply(names(site_estimates), function(site) {
    list(influence = site_fits[[site]]$varphi_ot,
         weight = site_estimates[[site]]$n / total_n)
    }),
    names(site_estimates)
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = var_re, blocks = boot_blocks,
    variance_method = variance_method, n_bootstrap = n_bootstrap
  )

  return(list(
    estimate = weighted_sum,
    variance = var_res$variance,
    se = var_res$se,
    method = "sample_size",
    n = total_n,
    influence_blocks = boot_blocks,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_bootstrap = var_res$n_bootstrap,
      var_fixed_effects = as.numeric(var_fe),
      var_random_effects = as.numeric(var_re),
      tau_squared = as.numeric(tau_sq),
      Q_statistic = het$Q,
      df = K - 1,
      I_squared = het$I_squared
    )
  ))
}

#' Inverse-variance weighted estimator (IVW) with heterogeneity correction
#' 
#' This estimator computes an inverse-variance weighted average of site-specific
#' AIPW estimates. By default (\code{variance_method = "bootstrap"}) the standard
#' error is the multiplier (wild) bootstrap, i.e. the fixed-effects sampling variance
#' \eqn{1/\sum_j \mathrm{precision}_j}; the random-effects (\eqn{\tau^2}) variance is
#' available via \code{variance_method = "analytic"} and retained in \code{components}.
#'
#' The analytic random-effects path is motivated as follows. When sites have different
#' covariate distributions, they may estimate slightly different population parameters,
#' and the fixed-effects IVW variance (\eqn{1/\sum_j \mathrm{precision}_j}) would then
#' understate uncertainty ABOUT THAT COMMON PARAMETER; Cochran's Q estimates the
#' between-site component (\eqn{\tau^2}). For the fixed target estimand, however, that
#' between-site spread is transport bias rather than sampling variance, so the bootstrap
#' (fixed-effects) SE is the honest default and the \eqn{\tau^2} inflation is opt-in.
#' 
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances
#'   for variance-valid inference (use_rcal is ignored in this mode).
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param A_val Treatment arm, either 0 or 1.
#' @param site_fits Optional precomputed per-site outputs from
#'   \code{fit_site_aipw}. When provided, avoids refitting nuisances.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws when
#'   \code{variance_method = "bootstrap"}.
#' @return estimate with variance
#' @export
estimate_inverse_variance_weighted <- function(data_split, family = "binomial",
                                               use_rcal = FALSE, use_crossfit = TRUE,
                                               n_folds = NULL, A_val = 1L,
                                               site_fits = NULL,
                                               variance_method = c("bootstrap", "analytic"),
                                               n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  sites <- names(data_split)
  site_results <- list()
  total_n <- 0

  if (is.null(site_fits)) {
    site_fits <- .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val
    )
  }
  
  for (site in sites) {
    res <- site_fits[[site]]
    safe_variance <- max(res$variance, VARIANCE_MIN)
    safe_theta <- max(res$V_ot, VARIANCE_MIN)
    site_results[[site]] <- list(
      estimate = res$estimate,
      variance = safe_variance,
      precision = 1 / safe_variance,
      theta = safe_theta,
      n = res$n
    )
    total_n <- total_n + res$n
  }
  
  # ==========================================================================
  # VARIANCE ESTIMATION WITH HETEROGENEITY CORRECTION
  # ==========================================================================
  # Standard IVW (fixed-effects) assumes all sites estimate the same parameter.
  # When there's covariate shift, sites may estimate different parameters,
  # introducing between-site heterogeneity that inflates variance.
  #
  # We use a random-effects approach:
  #   μ̂_j ~ N(μ, Var_j + τ²)
  # where τ² is the between-site variance estimated via DerSimonian-Laird.
  #
  # Random-effects variance:
  #   Var_RE = 1 / Σ_j (1 / (Var_j + τ²))
  # ==========================================================================
  
  K <- length(site_results)
  
  # Calculate fixed-effects IVW weights
  total_precision <- sum(sapply(site_results, function(x) x$precision))
  
  # Fixed-effects weighted estimate
  weighted_sum <- 0
  for (site in names(site_results)) {
    w <- site_results[[site]]$precision / total_precision
    weighted_sum <- weighted_sum + w * site_results[[site]]$estimate
  }
  
  # Fixed-effects variance (base variance)
  var_fe <- 1 / total_precision
  
  # ==========================================================================
  # Estimate between-site heterogeneity using DerSimonian-Laird
  # ==========================================================================
  site_est_vec <- sapply(site_results, function(x) x$estimate)
  site_var_vec <- sapply(site_results, function(x) x$variance)
  het <- calculate_dl_heterogeneity(site_est_vec, site_var_vec)
  tau_sq <- het$tau_sq
  
  # ==========================================================================
  # Random-effects variance
  # Var_RE = 1 / Σ_j (1 / (Var_j + τ²))
  # ==========================================================================
  if (tau_sq > 0) {
    # Random-effects weights
    re_precision <- 0
    for (site in names(site_results)) {
      var_j <- site_results[[site]]$variance
      re_precision <- re_precision + 1 / (var_j + tau_sq)
    }
    var_re <- 1 / re_precision
    
    # Re-compute weighted estimate with random-effects weights
    weighted_sum_re <- 0
    for (site in names(site_results)) {
      var_j <- site_results[[site]]$variance
      w_re <- (1 / (var_j + tau_sq)) / re_precision
      weighted_sum_re <- weighted_sum_re + w_re * site_results[[site]]$estimate
    }
    
    # Use random-effects estimate and variance
    final_estimate <- weighted_sum_re
    final_variance <- var_re
  } else {
    # No heterogeneity detected, use fixed-effects
    final_estimate <- weighted_sum
    final_variance <- var_fe
  }
  
  # ==========================================================================
  # Additional pooled IF variance as a check (like sample_size method)
  # This provides an alternative variance estimate based on pooled IFs
  # ==========================================================================
  
  # Compute weights for pooled variance (use RE weights if tau_sq > 0)
  weights <- numeric(K)
  for (idx in seq_along(site_results)) {
    site <- names(site_results)[idx]
    if (tau_sq > 0) {
      var_j <- site_results[[site]]$variance
      re_precision_total <- sum(sapply(site_results, function(x) 1/(x$variance + tau_sq)))
      weights[idx] <- (1 / (var_j + tau_sq)) / re_precision_total
    } else {
      weights[idx] <- site_results[[site]]$precision / total_precision
    }
  }
  
  # Pooled IF variance with between-site component
  within_site_var <- 0
  between_site_var <- 0
  
  for (idx in seq_along(site_results)) {
    site <- names(site_results)[idx]
    w_j <- weights[idx]
    n_j <- site_results[[site]]$n
    theta_j <- site_results[[site]]$theta
    mu_j <- site_results[[site]]$estimate
    
    # Within-site: w_j² × Var_j = w_j² × θ_j / n_j
    within_site_var <- within_site_var + w_j^2 * theta_j / n_j
    
    # Between-site: w_j² × (μ̂_j - μ̂)²
    between_site_var <- between_site_var + w_j^2 * (mu_j - final_estimate)^2
  }
  
  pooled_variance <- within_site_var + between_site_var
  
  # Use the larger of RE variance and pooled variance for robustness (analytic path)
  final_variance <- max(final_variance, pooled_variance)

  # Variance: bootstrap (default) over per-site influence blocks with fixed-effects
  # precision weights, or the analytic random-effects variance above.
  boot_blocks <- .name_influence_blocks(
    lapply(seq_along(site_results), function(idx) {
    site <- names(site_results)[idx]
    list(influence = site_fits[[site]]$varphi_ot,
         weight = weights[idx])
    }),
    names(site_results)
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = final_variance, blocks = boot_blocks,
    variance_method = variance_method, n_bootstrap = n_bootstrap
  )

  return(list(
    estimate = final_estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "inverse_variance",
    n = total_n,
    influence_blocks = boot_blocks,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_bootstrap = var_res$n_bootstrap,
      var_fixed_effects = var_fe,
      var_random_effects = if (tau_sq > 0) var_re else var_fe,
      var_pooled = pooled_variance,
      tau_squared = tau_sq,
      Q_statistic = het$Q,
      df = K - 1,
      I_squared = het$I_squared,
      within_site_var = within_site_var,
      between_site_var = between_site_var
    )
  ))
}

#' Tilted AIPW estimator with unpenalized MLE nuisances
#' @param data_split split data by site
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#'        Supports both "binomial" (logit) and "gaussian" (identity).
#' @param A_val Treatment arm, either 0 or 1.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws when
#'   \code{variance_method = "bootstrap"}.
#' @return estimate with variance
#' @export
estimate_tilted_aipw <- function(data_split, family = "binomial", A_val = 1L,
                                 variance_method = c("bootstrap", "analytic"),
                                 n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT) {
  variance_method <- match.arg(variance_method)
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)
  glm_spec <- resolve_glm_family(family)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")

  # Target-only estimate using unpenalized MLE for strict variance derivation.
  # validate_algorithm_inputs() already guarantees at least one A == A_val unit
  # in every site, so the inline emptiness checks live in the validator now.
  target_x <- as.matrix(target_data$W_outcome)
  target_y <- as.numeric(target_data$Y)
  target_a <- as.numeric(target_data$A)
  target_n <- target_data$n

  # Fit nuisance models and compute AIPW estimate
  target_est <- tryCatch({
    ps_fit_t <- fit_logit_mle(target_x, target_a)
    or_fit_t <- fit_glm_mle(target_x[target_a == A_val, , drop = FALSE],
                            target_y[target_a == A_val],
                            family = glm_spec$family)
    eta_t <- as.numeric(cbind(1, target_x) %*% or_fit_t$coefficients)
    m_hat_t <- switch(glm_spec$link,
                      "logit" = logistic(eta_t),
                      "identity" = eta_t,
                      eta_t)
    m_hat_t <- clip_outcome_pred(m_hat_t, glm_spec$family)
    target_aipw <- calculate_aipw_influence(target_y, target_a, target_x, m_hat_t,
                                            ps_fit_t$fitted, A_val = A_val,
                                            family = glm_spec$family)
    list(
      estimate = target_aipw$estimate,
      variance = target_aipw$variance,
      varphi_ot = target_aipw$influence
    )
  }, error = function(e) {
    stop(sprintf("estimate_tilted_aipw: target-site AIPW nuisance fitting failed: %s",
                 conditionMessage(e)), call. = FALSE)
  })

  # Calculate source site estimates with unpenalized MLE nuisances
  source_estimates <- list()
  Z_target <- as.matrix(target_data$Z_site)
  mean_phi_target <- c(1, colMeans(Z_target))
  # Pre-compute target-side matrices (constant across source sites)
  Z_target_int <- cbind(1, Z_target)
  Z_target_centered <- sweep(Z_target_int, 2, mean_phi_target, "-")
  
  for (site in source_sites) {
    source_data <- data_split[[site]]
    
    y_source <- source_data$Y
    tr_source <- source_data$A
    x_source <- as.matrix(source_data$W_outcome)
    n_source <- source_data$n

    # validate_algorithm_inputs() ensures every source site has >= 1 unit with
    # A == A_val, so treated_idx is guaranteed non-empty here.
    treated_idx <- which(tr_source == A_val)

    # Fit nuisance models and compute weighted AIPW estimate
    source_estimates[[site]] <- tryCatch({
      X_matrix <- as.matrix(x_source)
      Z_source <- as.matrix(source_data$Z_site)

      # Step 1: Fit outcome model (unpenalized MLE)
      or_fit <- fit_glm_mle(X_matrix[treated_idx, , drop = FALSE],
                as.numeric(y_source[treated_idx]),
                family = glm_spec$family)
      eta_source <- as.numeric(cbind(1, X_matrix) %*% or_fit$coefficients)
      m1_pred <- switch(glm_spec$link,
            "logit" = logistic(eta_source),
            "identity" = eta_source,
            eta_source)
      m1_pred <- clip_outcome_pred(m1_pred, glm_spec$family)

      # Step 2: Fit propensity model (unpenalized MLE)
      ps_fit <- fit_logit_mle(X_matrix, as.numeric(tr_source))
      prop_scores <- clip_propensity(ps_fit$fitted)

      # Step 3: Fit density ratio model via exponential tilting (unpenalized)
      # γ_{s,A_val} is the arm-specific density ratio
      alpha <- fit_initial_density_ratio(Z_source, tr_source, mean_phi_target, lambda = 0.0,
                                         A_val = A_val)

      # Density ratio weights
      Z_int <- cbind(1, Z_source)
      eta_i <- as.numeric(Z_int %*% alpha)
      w_i <- .normalize_log_weights(-eta_i, "estimate_tilted_aipw")

      # Step 4: AIPW estimate with nuisance-adjusted influence
      aipw_res <- calculate_aipw_influence(as.numeric(y_source), as.numeric(tr_source),
                                         X_matrix, m1_pred, prop_scores, w_i,
                                         A_val = A_val, family = glm_spec$family)

      weighted_estimate <- aipw_res$estimate
      phi_i <- aipw_res$phi

      # Density ratio adjustment for influence function
      d_alpha <- aipw_res$d
      phi_centered <- phi_i - weighted_estimate
      
      # Center source Z by mean_phi_target for consistent score function
      # This ensures source and target influence components are on same scale
      Z_int_centered <- sweep(Z_int, 2, mean_phi_target, "-")
      
      A_s <- -1 / d_alpha * colMeans(as.numeric(w_i) * Z_int_centered * as.numeric(phi_centered))
      # Use A-weighted score to match the density ratio estimating equation
      aw_i <- as.numeric(tr_source == A_val) * as.numeric(w_i)
      M_alpha <- t(Z_int_centered) %*% (Z_int_centered * aw_i) / n_source
      adj_alpha <- solve_with_ridge(M_alpha) %*% A_s
      infl_alpha_source <- as.numeric(aw_i * (Z_int_centered %*% adj_alpha))

      influence <- aipw_res$influence + infl_alpha_source
      weighted_variance <- mean(influence^2) / n_source

      # Target-side influence component for cross-site covariance
      target_if_component <- -as.numeric(Z_target_centered %*% adj_alpha)

      # Return the result (will be assigned to source_estimates[[site]])
      list(
        estimate = weighted_estimate,
        variance = weighted_variance,
        varphi_ot = influence,
        weights = w_i,
        n = n_source,
        A_s = A_s,
        alpha = alpha,
        d_alpha = d_alpha,
        target_if_component = target_if_component
      )
      
    }, error = function(e) {
      stop(sprintf("estimate_tilted_aipw: source site '%s' nuisance fitting failed: %s",
                   site, conditionMessage(e)), call. = FALSE)
    })
  }
  
  # Aggregate using sample size weights
  total_n <- target_data$n + sum(sapply(source_estimates, function(x) x$n))
  
  # Calculate sample size weights
  w_target <- target_data$n / total_n
  w_sources <- sapply(source_estimates, function(x) x$n / total_n)
  
  # Aggregated estimate
  weighted_sum <- w_target * target_est$estimate
  
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    weighted_sum <- weighted_sum + w * source_estimates[[site]]$estimate
  }
  
  # ============================================================
  # Calculate variance for tilted AIPW using influence functions
  # ============================================================
  # The aggregated estimator is a weighted sum of target-only and source
  # estimators. We combine target-side IF components (shared target data)
  # and source-side IF components (site-specific data) to capture
  # cross-site covariance induced by shared target information.
  # ============================================================
  
  target_if <- target_est$varphi_ot
  if (length(target_if) != target_n) {
    target_if <- rep(0, target_n)
  }
  
  aggregated_target_if <- w_target * target_if
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    target_if_component <- source_estimates[[site]]$target_if_component
    if (is.null(target_if_component) || length(target_if_component) != target_n) {
      target_if_component <- rep(0, target_n)
    }
    aggregated_target_if <- aggregated_target_if + w * target_if_component
  }
  
  aggregated_target_if <- aggregated_target_if - mean(aggregated_target_if)
  var_target_component <- mean(aggregated_target_if^2) / max(1, target_n)
  
  var_source_component <- 0
  for (site in names(source_estimates)) {
    w <- source_estimates[[site]]$n / total_n
    n_source <- source_estimates[[site]]$n
    source_if <- source_estimates[[site]]$varphi_ot
    if (!is.null(source_if) && length(source_if) > 0 && n_source > 0) {
      source_if <- source_if - mean(source_if)
      var_source_component <- var_source_component +
        (w^2) * mean(source_if^2) / n_source
    }
  }
  
  var_total <- max(var_target_component + var_source_component, VARIANCE_MIN)

  # Variance: bootstrap (default) over the shared-target influence block plus per-source
  # influence blocks, or the analytic influence-function variance above.
  boot_blocks <- .name_influence_blocks(
    c(
      list(list(influence = aggregated_target_if, weight = 1)),
      lapply(names(source_estimates), function(site) {
      list(influence = source_estimates[[site]]$varphi_ot,
           weight = source_estimates[[site]]$n / total_n)
      })
    ),
    c("t", names(source_estimates))
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = var_total, blocks = boot_blocks,
    variance_method = variance_method, n_bootstrap = n_bootstrap
  )

  return(list(
    estimate = weighted_sum,
    variance = var_res$variance,
    se = var_res$se,
    method = "tilted_aipw",
    n = total_n,
    influence_blocks = boot_blocks,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_bootstrap = var_res$n_bootstrap,
      target_weight = w_target,
      source_weights = w_sources,
      target_variance_component = var_target_component,
      source_variance_component = var_source_component
    )
  ))
}


#' Federated DR-AIPW estimator
#' 
#' Each source site uses density ratio weighting to estimate E_t[Y(1)].
#' All sites estimate the SAME estimand, then aggregate with IVW.
#' This is the theoretically correct federated baseline.
#' 
#' @param data_split Split data by site
#' @param dr_lambda Regularization for density ratio. If NULL, selected via CV.
#' @param dr_lambda_rule CV selection rule when \code{dr_lambda = NULL}:
#'   \code{"min"} (default) selects \code{lambda.min}; \code{"1se"} selects
#'   \code{lambda.1se}.
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment arm, either 0 or 1.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws when
#'   \code{variance_method = "bootstrap"}.
#' @param dr_weights_by_site Optional named list of precomputed source-to-target
#'   density-ratio weights.
#' @param dr_site_components Optional reusable site-level weighted-AIPW and
#'   density-ratio correction components.
#' @param n_cores Number of independent source density-ratio fits to run in
#'   parallel when \code{dr_weights_by_site} is not supplied.
#' @return List with estimate, variance, se
#' @export
estimate_federated_dr <- function(data_split, dr_lambda = NULL,
                                  dr_lambda_rule = c("min", "1se"),
                                  A_val = 1L, family = "binomial",
                                  variance_method = c("bootstrap", "analytic"),
                                  n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
                                  dr_weights_by_site = NULL,
                                  dr_site_components = NULL,
                                  n_cores = 1L) {
  variance_method <- match.arg(variance_method)
  dr_lambda_rule <- .match_nuisance_lambda_rule(
    dr_lambda_rule, "estimate_federated_dr", arg = "dr_lambda_rule"
  )
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  n_target <- target_data$n
  dr_weights_by_site <- .resolve_dr_weights_by_site(
    data_split = data_split,
    dr_lambda = dr_lambda,
    dr_lambda_rule = dr_lambda_rule,
    n_cores = n_cores,
    dr_weights_by_site = dr_weights_by_site,
    caller = "estimate_federated_dr"
  )
  dr_weight_diagnostics <-
    .density_ratio_clipping_diagnostics(dr_weights_by_site)
  
  site_results <- .resolve_dr_site_components(
    data_split = data_split,
    dr_weights_by_site = dr_weights_by_site,
    family = family,
    A_val = A_val,
    required_sites = c("t", source_sites),
    n_cores = n_cores,
    dr_site_components = dr_site_components,
    caller = "estimate_federated_dr"
  )
  
  # Filter valid sites and aggregate with IVW
  valid_sites <- names(site_results)[sapply(site_results, function(x) 
    !is.na(x$estimate) && is.finite(x$variance) && x$variance > 0)]
  
  if (length(valid_sites) == 0) {
    return(list(estimate = NA, variance = Inf, se = Inf, method = "federated_dr"))
  }
  
  precisions <- sapply(site_results[valid_sites], function(x) 1/x$variance)
  total_precision <- sum(precisions)
  ivw_weights <- precisions / total_precision
  
  estimate <- sum(ivw_weights * sapply(site_results[valid_sites], function(x) x$estimate))
  
  # Shared-target + source-side variance decomposition (same style as main methods)
  weight_by_site <- setNames(as.numeric(ivw_weights), valid_sites)
  aggregated_target_if <- rep(0, n_target)
  
  if ("t" %in% valid_sites) {
    aggregated_target_if <- aggregated_target_if + weight_by_site[["t"]] * site_results[["t"]]$varphi_ot
  }
  for (site in setdiff(valid_sites, "t")) {
    target_if_component <- site_results[[site]]$target_if_component
    if (!is.null(target_if_component) && length(target_if_component) == n_target) {
      aggregated_target_if <- aggregated_target_if + weight_by_site[[site]] * target_if_component
    }
  }
  aggregated_target_if <- aggregated_target_if - mean(aggregated_target_if)
  var_target_component <- mean(aggregated_target_if^2) / max(1, n_target)

  var_source_component <- 0
  for (site in setdiff(valid_sites, "t")) {
    w_site <- weight_by_site[[site]]
    n_site <- site_results[[site]]$n
    varphi_site <- site_results[[site]]$varphi_ot
    if (!is.null(varphi_site) && length(varphi_site) > 0) {
      varphi_site <- varphi_site - mean(varphi_site)
      var_source_component <- var_source_component + w_site^2 * mean(varphi_site^2) / max(1, n_site)
    }
  }

  variance <- max(var_target_component + var_source_component, VARIANCE_MIN)

  # Variance: bootstrap (default) over the shared-target influence block plus per-source
  # influence blocks (IVW weights), or the analytic influence-function variance above.
  source_block_sites <- setdiff(valid_sites, "t")
  boot_blocks <- .name_influence_blocks(
    c(
      list(list(influence = aggregated_target_if, weight = 1)),
      lapply(source_block_sites, function(site) {
      list(influence = site_results[[site]]$varphi_ot,
           weight = weight_by_site[[site]])
      })
    ),
    c("t", source_block_sites)
  )
  var_res <- .resolve_comparison_variance(
    analytic_variance = variance, blocks = boot_blocks,
    variance_method = variance_method, n_bootstrap = n_bootstrap
  )

  return(list(
    estimate = estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "federated_dr",
    n = sum(sapply(site_results[valid_sites], function(x) x$n)),
    influence_blocks = boot_blocks,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_bootstrap = var_res$n_bootstrap,
      ivw_weights = weight_by_site,
      var_target_component = var_target_component,
      var_source_component = var_source_component,
      dr_weight_diagnostics = dr_weight_diagnostics
    )
  ))
}

#' Pooled DR-AIPW estimator (centralized version)
#' 
#' Pools all data and uses density ratio weighting for each observation.
#' This is a centralized (non-federated) method for comparison.
#' 
#' @param data_split Split data by site
#' @param dr_lambda Regularization for density ratio. If NULL, selected via CV.
#' @param dr_lambda_rule CV selection rule when \code{dr_lambda = NULL}:
#'   \code{"min"} (default) selects \code{lambda.min}; \code{"1se"} selects
#'   \code{lambda.1se}.
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment arm, either 0 or 1.
#' @param variance_method Standard-error method: \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap over influence-function blocks; \code{"analytic"} uses
#'   the method's analytic variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws when
#'   \code{variance_method = "bootstrap"}.
#' @param dr_weights_by_site Optional named list of precomputed source-to-target
#'   density-ratio weights.
#' @param dr_site_components Optional reusable site-level weighted-AIPW and
#'   density-ratio correction components.
#' @param n_cores Number of independent source density-ratio fits to run in
#'   parallel when \code{dr_weights_by_site} is not supplied.
#' @return List with estimate, variance, se
#' @export
estimate_pooled_dr <- function(data_split, dr_lambda = NULL,
                               dr_lambda_rule = c("min", "1se"),
                               A_val = 1L, family = "binomial",
                               variance_method = c("bootstrap", "analytic"),
                               n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
                               dr_weights_by_site = NULL,
                               dr_site_components = NULL,
                               n_cores = 1L) {
  variance_method <- match.arg(variance_method)
  dr_lambda_rule <- .match_nuisance_lambda_rule(
    dr_lambda_rule, "estimate_pooled_dr", arg = "dr_lambda_rule"
  )
  validate_algorithm_inputs(data_split, family = family, A_val = A_val)

  target_data <- data_split[["t"]]
  source_sites <- setdiff(names(data_split), "t")
  n_target <- target_data$n
  dr_weights_by_site <- .resolve_dr_weights_by_site(
    data_split = data_split,
    dr_lambda = dr_lambda,
    dr_lambda_rule = dr_lambda_rule,
    n_cores = n_cores,
    dr_weights_by_site = dr_weights_by_site,
    caller = "estimate_pooled_dr"
  )
  dr_weight_diagnostics <-
    .density_ratio_clipping_diagnostics(dr_weights_by_site)
  dr_site_components <- .resolve_dr_site_components(
    data_split = data_split,
    dr_weights_by_site = dr_weights_by_site,
    family = family,
    A_val = A_val,
    required_sites = source_sites,
    n_cores = n_cores,
    dr_site_components = dr_site_components,
    caller = "estimate_pooled_dr"
  )
  
  # Pool all data
  all_y <- target_data$Y
  all_a <- target_data$A
  all_X <- as.matrix(target_data$W_outcome)
  all_weights <- rep(1, target_data$n)
  all_sites <- rep("t", target_data$n)
  source_meta <- list()
  row_cursor <- target_data$n
  
  for (site in source_sites) {
    source_data <- data_split[[site]]
    dr_weights <- dr_weights_by_site[[site]]
    source_component <- dr_site_components[[site]]
    n_source <- source_data$n
    
    idx_start <- row_cursor + 1L
    idx_end <- row_cursor + n_source
    row_cursor <- idx_end
    source_meta[[site]] <- list(
      idx = idx_start:idx_end,
      infl_alpha_source = source_component$source_if_correction,
      target_if_component = source_component$target_if_component,
      lambda = as.numeric(
        source_component$lambda %||%
          attr(dr_weights, "lambda_used") %||% dr_lambda %||% NA_real_
      )
    )
    
    all_y <- c(all_y, source_data$Y)
    all_a <- c(all_a, source_data$A)
    all_X <- rbind(all_X, as.matrix(source_data$W_outcome))
    all_weights <- c(all_weights, dr_weights)
    all_sites <- c(all_sites, rep(site, n_source))
  }
  
  result <- calculate_weighted_site_aipw(y = all_y, a = all_a, X = all_X, weights = all_weights, family = family, A_val = A_val)

  # Stacked IF: base weighted AIPW IF + DR-weight estimation correction
  N_all <- length(all_y)
  total_correction <- rep(0, N_all)
  for (site in names(source_meta)) {
    meta <- source_meta[[site]]
    total_correction[meta$idx] <- total_correction[meta$idx] + meta$infl_alpha_source
    total_correction[seq_len(n_target)] <- total_correction[seq_len(n_target)] + meta$target_if_component
  }
  varphi_total <- as.numeric(result$varphi_ot) + total_correction

  # Site-stratified variance (same style as main cross-fitting estimators):
  #   Var = (1/N_all^2) * sum_g sum_{i in g} (varphi_i - mean_g(varphi))^2
  wss <- 0
  for (site in unique(all_sites)) {
    idx <- which(all_sites == site)
    if (length(idx) == 0) next
    varphi_g <- varphi_total[idx]
    varphi_g <- varphi_g - mean(varphi_g)
    wss <- wss + sum(varphi_g^2)
  }
  final_variance <- max(wss / (N_all^2), VARIANCE_MIN)

  # Variance: bootstrap (default) over per-site blocks of the stacked influence function
  # (one block per site, weight n_site/N_all so the wild bootstrap reproduces the
  # site-stratified within-group variance above), or that analytic variance.
  block_sites <- unique(all_sites)
  boot_blocks <- .name_influence_blocks(lapply(block_sites, function(site) {
    idx <- which(all_sites == site)
    list(influence = varphi_total[idx], weight = length(idx) / N_all)
  }), block_sites)
  var_res <- .resolve_comparison_variance(
    analytic_variance = final_variance, blocks = boot_blocks,
    variance_method = variance_method, n_bootstrap = n_bootstrap
  )

  source_lambdas <- setNames(
    sapply(source_sites, function(site) {
      meta <- source_meta[[site]]
      if (is.null(meta)) NA_real_ else as.numeric(meta$lambda)
    }),
    source_sites
  )
  
  return(list(
    estimate = result$estimate,
    variance = var_res$variance,
    se = var_res$se,
    method = "pooled_dr",
    n = N_all,
    varphi_ot = varphi_total,
    influence_blocks = boot_blocks,
    components = list(
      variance_method = var_res$variance_method,
      variance_analytic = var_res$variance_analytic,
      se_analytic = var_res$se_analytic,
      se_bootstrap = var_res$se_bootstrap,
      n_bootstrap = var_res$n_bootstrap,
      n_target = n_target,
      n_source_total = N_all - n_target,
      mean_weight = mean(all_weights),
      source_lambdas = source_lambdas,
      variance_wss = wss,
      dr_weight_diagnostics = dr_weight_diagnostics
    )
  ))
}

#' Run all comparison methods
#' 
#' Returns both:
#' - Naive methods (sample_size, inverse_variance): may have bias due to covariate shift
#' - DR-corrected methods (federated_dr, pooled_dr): theoretically correct
#' 
#' @param data_split split data by site
#' @param use_rcal Logical. If TRUE, use RCAL. If FALSE (default), use glmnet.
#' @param use_crossfit Logical. If TRUE (default), use cross-fitted nuisances.
#' @param n_folds Number of cross-fitting folds (default uses data-driven value).
#' @param family GLM family ("binomial", "gaussian", etc.). Default "binomial".
#' @param A_val Treatment arm, either 0 or 1.
#' @param variance_method Standard-error method for the comparison baselines
#'   (\code{sample_size}, \code{inverse_variance}, \code{federated_dr},
#'   \code{pooled_dr}, \code{tilted_aipw}): \code{"bootstrap"} (default) uses the
#'   multiplier (wild) bootstrap; \code{"analytic"} uses each method's analytic
#'   variance. \code{target_only} always uses its analytic influence-function variance.
#' @param n_bootstrap Number of multiplier-bootstrap draws used when
#'   \code{variance_method = "bootstrap"}. The production default is 5,000.
#' @param include_tilted Whether the optional unpenalized tilted-AIPW baseline
#'   may be included. It is normally disabled in high-dimensional studies.
#' @return Named list of per-method results. The five comparison baselines additionally
#'   carry \code{components$variance_method}, \code{components$variance_analytic},
#'   \code{components$se_analytic}, and \code{components$se_bootstrap}; the
#'   \code{target_only} benchmark does not. Consumers should read \code{$estimate} and
#'   \code{$se} uniformly and treat those four \code{components} fields as baseline-only.
#' @param methods Optional character vector selecting methods to compute. The
#'   default \code{NULL} computes all available methods (subject to
#'   \code{include_tilted}). Restricting this vector avoids unnecessary nuisance
#'   fits in large simulation studies.
#' @param n_cores Number of cores used to fit the independent site-specific
#'   AIPW nuisances shared by SS and IVW. Defaults to one.
#' @param dr_weights_by_site Optional named list of source-to-target
#'   density-ratio weights shared by Federated-DR and Pooled-DR.
#' @param dr_site_components Optional precomputed site-level DR components
#'   shared by Federated-DR and Pooled-DR for the requested treatment arm.
#' @export
run_all_comparisons <- function(data_split, use_rcal = FALSE,
                                use_crossfit = TRUE, n_folds = NULL,
                                family = "binomial", A_val = 1L,
                                variance_method = c("bootstrap", "analytic"),
                                n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
                                include_tilted = TRUE,
                                methods = NULL,
                                n_cores = 1L,
                                dr_weights_by_site = NULL,
                                dr_site_components = NULL) {
  variance_method <- match.arg(variance_method)
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, "run_all_comparisons"
  )
  available <- c(
    "target_only", "sample_size", "inverse_variance",
    "federated_dr", "pooled_dr", "tilted_aipw"
  )
  requested <- if (is.null(methods)) available else unique(as.character(methods))
  unknown <- setdiff(requested, available)
  if (length(unknown) > 0L) {
    stop(sprintf(
      "run_all_comparisons: unknown method(s): %s.",
      paste(unknown, collapse = ", ")
    ))
  }
  if (!isTRUE(include_tilted)) {
    requested <- setdiff(requested, "tilted_aipw")
  }

  site_fit_methods <- c("sample_size", "inverse_variance")
  precomputed_site_fits <- if (any(requested %in% site_fit_methods)) {
    .fit_site_aipw_all_sites(
      data_split = data_split,
      family = family,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      A_val = A_val,
      n_cores = n_cores
    )
  } else {
    NULL
  }
  dr_methods <- c("federated_dr", "pooled_dr")
  shared_dr_weights <- if (any(requested %in% dr_methods)) {
    .resolve_dr_weights_by_site(
      data_split = data_split,
      n_cores = n_cores,
      dr_weights_by_site = dr_weights_by_site,
      caller = "run_all_comparisons"
    )
  } else {
    NULL
  }
  shared_dr_components <- if (any(requested %in% dr_methods)) {
    required_dr_sites <- if ("federated_dr" %in% requested) {
      names(data_split)
    } else {
      setdiff(names(data_split), "t")
    }
    .resolve_dr_site_components(
      data_split = data_split,
      dr_weights_by_site = shared_dr_weights,
      family = family,
      A_val = A_val,
      required_sites = required_dr_sites,
      n_cores = n_cores,
      dr_site_components = dr_site_components,
      caller = "run_all_comparisons"
    )
  } else {
    NULL
  }

  results <- list()
  if ("target_only" %in% requested) {
    target_result <- if (!is.null(precomputed_site_fits[["t"]])) {
      precomputed_site_fits[["t"]]
    } else {
      estimate_target_only(
        data_split, family, use_rcal = use_rcal,
        use_crossfit = use_crossfit, n_folds = n_folds, A_val = A_val
      )
    }
    target_result$influence_blocks <- list(
      t = list(influence = target_result$varphi_ot, weight = 1)
    )
    results$target_only <- target_result
  }
  if ("sample_size" %in% requested) {
    results$sample_size <- estimate_sample_size_weighted(
      data_split, family, use_rcal = use_rcal,
      use_crossfit = use_crossfit, n_folds = n_folds,
      A_val = A_val, site_fits = precomputed_site_fits,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap
    )
  }
  if ("inverse_variance" %in% requested) {
    results$inverse_variance <- estimate_inverse_variance_weighted(
      data_split, family, use_rcal = use_rcal,
      use_crossfit = use_crossfit, n_folds = n_folds,
      A_val = A_val, site_fits = precomputed_site_fits,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap
    )
  }
  if ("federated_dr" %in% requested) {
    results$federated_dr <- estimate_federated_dr(
      data_split, dr_lambda = NULL, A_val = A_val, family = family,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap,
      dr_weights_by_site = shared_dr_weights,
      dr_site_components = shared_dr_components,
      n_cores = n_cores
    )
  }
  if ("pooled_dr" %in% requested) {
    results$pooled_dr <- estimate_pooled_dr(
      data_split, dr_lambda = NULL, A_val = A_val, family = family,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap,
      dr_weights_by_site = shared_dr_weights,
      dr_site_components = shared_dr_components,
      n_cores = n_cores
    )
  }

  # Tilted AIPW uses unpenalized MLE nuisances and is degenerate in high
  # dimension (perfect separation / collinearity makes the MLE GLM fail), so it
  # is opt-in. The main manuscript omits it; set include_tilted = TRUE to keep it.
  if ("tilted_aipw" %in% requested) {
    results$tilted_aipw <- estimate_tilted_aipw(
      data_split, family, A_val = A_val, variance_method = variance_method,
      n_bootstrap = n_bootstrap)
  }

  return(results)
}

#' Run comparison methods for the target average treatment effect
#'
#' Fits the two treatment arms and forms each TATE estimate and standard error
#' from paired, site-aligned influence-function blocks. This retains the
#' treatment--control covariance within every site.
#'
#' @inheritParams run_all_comparisons
#' @param mu1_results Optional arm-1 results previously returned by
#'   \code{run_all_comparisons}; supplying them avoids refitting that arm.
#' @return Named list of TATE results for the requested methods.
#' @export
run_all_comparisons_tate <- function(
    data_split, use_rcal = FALSE, use_crossfit = TRUE, n_folds = NULL,
    family = "binomial",
    variance_method = c("bootstrap", "analytic"),
    n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
    include_tilted = TRUE, methods = NULL, mu1_results = NULL,
    n_cores = 1L, dr_weights_by_site = NULL) {
  variance_method <- match.arg(variance_method)
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, "run_all_comparisons_tate"
  )
  available <- c(
    "target_only", "sample_size", "inverse_variance",
    "federated_dr", "pooled_dr", "tilted_aipw"
  )
  requested <- if (is.null(methods)) available else unique(as.character(methods))
  unknown <- setdiff(requested, available)
  if (length(unknown) > 0L) {
    stop(sprintf(
      "run_all_comparisons_tate: unknown method(s): %s.",
      paste(unknown, collapse = ", ")
    ), call. = FALSE)
  }
  if (!isTRUE(include_tilted)) {
    requested <- setdiff(requested, "tilted_aipw")
  }
  dr_methods <- c("federated_dr", "pooled_dr")
  shared_dr_weights <- if (any(requested %in% dr_methods)) {
    .resolve_dr_weights_by_site(
      data_split = data_split,
      n_cores = n_cores,
      dr_weights_by_site = dr_weights_by_site,
      caller = "run_all_comparisons_tate"
    )
  } else {
    NULL
  }

  if (is.null(mu1_results)) {
    mu1_results <- run_all_comparisons(
      data_split = data_split,
      use_rcal = use_rcal,
      use_crossfit = use_crossfit,
      n_folds = n_folds,
      family = family,
      A_val = 1L,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap,
      include_tilted = include_tilted,
      methods = requested,
      n_cores = n_cores,
      dr_weights_by_site = shared_dr_weights
    )
  }
  missing_mu1 <- setdiff(requested, names(mu1_results))
  if (length(missing_mu1) > 0L) {
    stop(sprintf(
      "run_all_comparisons_tate: mu1_results is missing method(s): %s.",
      paste(missing_mu1, collapse = ", ")
    ), call. = FALSE)
  }

  mu0_results <- run_all_comparisons(
    data_split = data_split,
    use_rcal = use_rcal,
    use_crossfit = use_crossfit,
    n_folds = n_folds,
    family = family,
    A_val = 0L,
    variance_method = variance_method,
    n_bootstrap = n_bootstrap,
    include_tilted = include_tilted,
    methods = requested,
    n_cores = n_cores,
    dr_weights_by_site = shared_dr_weights
  )

  stats::setNames(lapply(requested, function(method) {
    .combine_comparison_arms(
      mu1_result = mu1_results[[method]],
      mu0_result = mu0_results[[method]],
      method = method,
      variance_method = variance_method,
      n_bootstrap = n_bootstrap
    )
  }), requested)
}
