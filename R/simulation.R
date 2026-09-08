# simulation.R - Simulation experiment functions for RoCE
#
# This file contains reusable simulation functions extracted from main.R.
# The orchestration script (main.R) calls these functions with experiment-
# specific configuration and checkpoint infrastructure.
#
# Contents:
#   1. run_single_simulation - Run one Monte Carlo replicate
#   2. summarize_results - Aggregate simulation results into summary statistics
#   3. run_simulation_study - Full study orchestration with checkpointing

.is_tate_method <- function(method) {
  grepl("_ate(?:_armwise|_hard_threshold|_quadratic_bias)?$", method, perl = TRUE)
}

.bind_sim_result_list <- function(rows) {
  if (length(rows) == 0) {
    return(data.frame())
  }
  keep <- vapply(rows, function(x) is.data.frame(x) && nrow(x) > 0,
                 logical(1))
  if (!any(keep)) {
    return(data.frame())
  }
  rows <- rows[keep]
  all_columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing_columns <- setdiff(all_columns, names(row))
    for (column in missing_columns) {
      row[[column]] <- NA
    }
    row[all_columns]
  })
  do.call(rbind, rows)
}

.face_heterogeneity_type <- function(
    ate_deviation, n_deviated_sites, effect_mod_strength = 0,
    deviation_mechanism = c("treated_arm", "both_arms")) {
  deviation_mechanism <- match.arg(deviation_mechanism)
  if (!is.numeric(ate_deviation) || length(ate_deviation) != 1L ||
      !is.finite(ate_deviation) || ate_deviation < 0 ||
      !is.numeric(n_deviated_sites) || length(n_deviated_sites) != 1L ||
      !is.finite(n_deviated_sites) || n_deviated_sites < 0 ||
      n_deviated_sites != floor(n_deviated_sites) ||
      !is.numeric(effect_mod_strength) || length(effect_mod_strength) != 1L ||
      !is.finite(effect_mod_strength)) {
    stop(
      paste0(
        ".face_heterogeneity_type: ate_deviation must be non-negative and ",
        "n_deviated_sites must be a non-negative integer, and ",
        "effect_mod_strength must be finite."
      ),
      call. = FALSE
    )
  }
  has_deviation <- ate_deviation > 0 && n_deviated_sites > 0
  has_effect_modification <- effect_mod_strength != 0
  shared_shift <- identical(deviation_mechanism, "both_arms")
  if (!has_deviation && !has_effect_modification) {
    return("none")
  }
  if (has_deviation && has_effect_modification) {
    return(if (shared_shift) "shared_shift_and_effect_modification"
           else "deviation_and_effect_modification")
  }
  if (has_effect_modification) {
    return("source_effect_modification")
  }
  if (n_deviated_sites == 1) {
    return(if (shared_shift) "one_shared_shift_source" else "one_deviated_source")
  }
  if (shared_shift) "multiple_shared_shift_sources" else "multiple_deviated_sources"
}

.abort_with_context <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

.warn_with_context <- function(fmt, ...) {
  warning(sprintf(fmt, ...), call. = FALSE)
}

.require_method_result <- function(results, method_name, required = c("estimate", "se")) {
  method_res <- results[[method_name]]
  if (is.null(method_res)) {
    .abort_with_context(
      "Method '%s' did not return a result. Check earlier warnings from its fitting routine.",
      method_name
    )
  }
  missing <- required[vapply(required, function(nm) is.null(method_res[[nm]]), logical(1))]
  if (length(missing) > 0L) {
    .abort_with_context(
      "Method '%s' returned an incomplete result; missing field(s): %s.",
      method_name, paste(missing, collapse = ", ")
    )
  }
  method_res
}

.comparison_variance_diagnostics <- function(method_result) {
  components <- method_result$components
  if (is.null(components) || is.null(components$variance_method)) {
    return(list(
      comparison_variance_method = NA_character_,
      comparison_se_analytic = NA_real_,
      comparison_se_bootstrap = NA_real_,
      comparison_n_bootstrap = NA_integer_,
      dr_weight_n = NA_integer_,
      dr_weight_n_clipped = NA_integer_,
      dr_weight_fraction_clipped = NA_real_,
      dr_weight_max_site_fraction_clipped = NA_real_,
      dr_weight_min_before_clipping = NA_real_,
      dr_weight_max_before_clipping = NA_real_
    ))
  }
  dr <- components$dr_weight_diagnostics %||% list()
  list(
    comparison_variance_method = as.character(components$variance_method),
    comparison_se_analytic = as.numeric(components$se_analytic),
    comparison_se_bootstrap = as.numeric(components$se_bootstrap),
    comparison_n_bootstrap = as.integer(components$n_bootstrap),
    dr_weight_n = as.integer(dr$dr_weight_n %||% NA_integer_),
    dr_weight_n_clipped = as.integer(
      dr$dr_weight_n_clipped %||% NA_integer_
    ),
    dr_weight_fraction_clipped = as.numeric(
      dr$dr_weight_fraction_clipped %||% NA_real_
    ),
    dr_weight_max_site_fraction_clipped = as.numeric(
      dr$dr_weight_max_site_fraction_clipped %||% NA_real_
    ),
    dr_weight_min_before_clipping = as.numeric(
      dr$dr_weight_min_before_clipping %||% NA_real_
    ),
    dr_weight_max_before_clipping = as.numeric(
      dr$dr_weight_max_before_clipping %||% NA_real_
    )
  )
}

.binary_cell_diagnostics <- function(data_split, family) {
  empty <- list(
    min_site_arm_outcome_cell_n = NA_integer_,
    min_target_arm_outcome_cell_n = NA_integer_,
    n_site_arm_outcome_cells_below_8 = NA_integer_
  )
  if (!identical(family, "binomial")) {
    return(empty)
  }
  if (!is.list(data_split) || is.null(data_split[["t"]]) ||
      any(vapply(data_split, function(site) {
        is.null(site$A) || is.null(site$Y) || length(site$A) != length(site$Y) ||
          anyNA(site$A) || anyNA(site$Y) ||
          any(!site$A %in% c(0, 1)) || any(!site$Y %in% c(0, 1))
      }, logical(1L)))) {
    stop(
      ".binary_cell_diagnostics: binomial sites require aligned binary A and Y.",
      call. = FALSE
    )
  }
  counts <- lapply(data_split, function(site) {
    as.integer(table(
      factor(site$A, levels = c(0, 1)),
      factor(site$Y, levels = c(0, 1))
    ))
  })
  all_counts <- unlist(counts, use.names = FALSE)
  list(
    min_site_arm_outcome_cell_n = as.integer(min(all_counts)),
    min_target_arm_outcome_cell_n = as.integer(min(counts[["t"]])),
    n_site_arm_outcome_cells_below_8 = as.integer(sum(all_counts < 8L))
  )
}

.make_simulation_result_row <- function(
    sim_id, method, estimate, se, truth,
    n_total, K, p, config, heterogeneity_type, estimand_type) {
  estimate <- as.numeric(estimate)
  se <- as.numeric(se)
  truth <- as.numeric(truth)
  if (length(estimate) != 1L || length(se) != 1L || length(truth) != 1L ||
      any(!is.finite(c(estimate, se, truth))) || se < 0) {
    .abort_with_context(
      "Cannot create result row for '%s': estimate, SE, and truth must be finite scalars with SE >= 0.",
      method
    )
  }

  ci_lower <- estimate - Z_ALPHA_05 * se
  ci_upper <- estimate + Z_ALPHA_05 * se
  data.frame(
    sim_id = sim_id,
    method = method,
    estimate = estimate,
    se = se,
    bias = estimate - truth,
    coverage = truth >= ci_lower && truth <= ci_upper,
    ci_width = ci_upper - ci_lower,
    n_total = n_total,
    K = K,
    p = p,
    config = config,
    heterogeneity_type = heterogeneity_type,
    estimand_type = estimand_type,
    stringsAsFactors = FALSE
  )
}

.validate_weight_bootstrap_replicates <- function(
    n_weight_bootstrap, caller = "run_single_simulation") {
  if (length(n_weight_bootstrap) != 1L ||
      !is.numeric(n_weight_bootstrap) || is.na(n_weight_bootstrap) ||
      !is.finite(n_weight_bootstrap) || n_weight_bootstrap < 0 ||
      n_weight_bootstrap != floor(n_weight_bootstrap) ||
      n_weight_bootstrap > .Machine$integer.max) {
    stop(sprintf(
      "%s: n_weight_bootstrap must be one non-negative integer.", caller
    ), call. = FALSE)
  }
  n_weight_bootstrap <- as.integer(n_weight_bootstrap)
  if (identical(n_weight_bootstrap, 1L)) {
    stop(sprintf(
      "%s: n_weight_bootstrap must be 0 (disabled) or at least 2.", caller
    ), call. = FALSE)
  }
  n_weight_bootstrap
}

# Derive an RNG stream without consuming or changing the caller's RNG state.
# Deliberately omit rho: a same-simulation rho group then uses common
# multiplier draws, preserving the paired design used by the FACE experiment.
.tate_weight_bootstrap_seed <- function(
    sim_id, screening_rule = "soft_penalty") {
  if (length(sim_id) != 1L || !is.numeric(sim_id) || is.na(sim_id) ||
      !is.finite(sim_id) || sim_id < 1 || sim_id != floor(sim_id)) {
    stop("sim_id must be one positive integer for weight-bootstrap seeding.",
         call. = FALSE)
  }
  screening_rule <- match.arg(
    screening_rule, c("soft_penalty", "hard_threshold", "quadratic_bias")
  )
  stream_offset <- switch(
    screening_rule,
    soft_penalty = 13007,
    hard_threshold = 17011,
    quadratic_bias = 19013
  )
  modulus <- as.double(.Machine$integer.max) - 1
  as.integer((as.double(sim_id) * 104729 + stream_offset) %% modulus + 1)
}

.run_tate_weight_bootstrap <- function(
    tate_result, n_weight_bootstrap, sim_id,
    screening_rule = "soft_penalty") {
  B <- .validate_weight_bootstrap_replicates(
    n_weight_bootstrap, ".run_tate_weight_bootstrap"
  )
  if (B == 0L) return(list())

  bootstrap_seed <- .tate_weight_bootstrap_seed(sim_id, screening_rule)
  fit <- estimate_tate_weight_bootstrap(
    tate_result = tate_result,
    B = B,
    seed = bootstrap_seed,
    relearn_weights = TRUE,
    screening_rule = screening_rule,
    save_draws = FALSE
  )
  required <- c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "n_bootstrap", "bootstrap_multiplier", "bootstrap_seed",
    "bootstrap_failures"
  )
  missing <- setdiff(required, names(fit))
  if (length(missing) > 0L) {
    stop(
      "estimate_tate_weight_bootstrap() omitted required fields: ",
      paste(missing, collapse = ", "), call. = FALSE
    )
  }
  numeric_scalars <- c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "n_bootstrap", "bootstrap_seed", "bootstrap_failures"
  )
  invalid_numeric <- vapply(numeric_scalars, function(name) {
    value <- fit[[name]]
    length(value) != 1L || !is.numeric(value) || !is.finite(value)
  }, logical(1L))
  if (any(invalid_numeric)) {
    stop(
      "estimate_tate_weight_bootstrap() returned invalid scalar fields: ",
      paste(numeric_scalars[invalid_numeric], collapse = ", "),
      call. = FALSE
    )
  }
  if (fit$variance_weight_relearn_bootstrap <= 0 ||
      fit$se_weight_relearn_bootstrap <= 0 ||
      fit$se_fixed_weight_bootstrap <= 0 ||
      fit$weight_uncertainty_ratio <= 0 ||
      fit$n_bootstrap != B || fit$bootstrap_seed != bootstrap_seed ||
      fit$bootstrap_failures < 0 || fit$bootstrap_failures > B ||
      fit$bootstrap_failures != floor(fit$bootstrap_failures) ||
      length(fit$bootstrap_multiplier) != 1L ||
      is.na(fit$bootstrap_multiplier) ||
      !as.character(fit$bootstrap_multiplier) %in%
        c("exponential", "site_stratified_exponential")) {
    stop("estimate_tate_weight_bootstrap() returned inconsistent metadata.",
         call. = FALSE)
  }

  list(
    variance_weight_relearn_bootstrap =
      as.numeric(fit$variance_weight_relearn_bootstrap),
    se_weight_relearn_bootstrap = as.numeric(fit$se_weight_relearn_bootstrap),
    se_fixed_weight_bootstrap = as.numeric(fit$se_fixed_weight_bootstrap),
    weight_uncertainty_ratio = as.numeric(fit$weight_uncertainty_ratio),
    weight_relearn_n_bootstrap = as.integer(fit$n_bootstrap),
    weight_bootstrap_multiplier = as.character(fit$bootstrap_multiplier),
    weight_bootstrap_seed = as.integer(fit$bootstrap_seed),
    weight_bootstrap_failures = as.integer(fit$bootstrap_failures),
    weight_bootstrap_relearn_weights = TRUE,
    weight_bootstrap_screening_rule = screening_rule
  )
}

.summarize_tate_crossfit_timing <- function(mu1_result, mu0_result) {
  arm_results <- list(mu1 = mu1_result, mu0 = mu0_result)
  if (any(vapply(
    arm_results,
    function(result) is.null(result$timing) || is.null(result$timing$folds),
    logical(1L)
  ))) {
    return(numeric(0))
  }

  fold_timing <- do.call(rbind, lapply(arm_results, function(result) {
    as.data.frame(result$timing$folds)
  }))
  additive_columns <- setdiff(names(fold_timing), "fold")
  fold_sums <- vapply(
    fold_timing[additive_columns],
    function(values) sum(as.numeric(values)),
    numeric(1L)
  )
  names(fold_sums) <- paste0("face_", names(fold_sums))

  c(
    face_mu1_total_seconds = as.numeric(mu1_result$timing$total_seconds),
    face_mu0_total_seconds = as.numeric(mu0_result$timing$total_seconds),
    face_sum_arm_seconds =
      as.numeric(mu1_result$timing$total_seconds) +
      as.numeric(mu0_result$timing$total_seconds),
    fold_sums
  )
}

.summarize_tate_nuisance_diagnostics <- function(mu1_result, mu0_result) {
  arm_results <- list(mu1 = mu1_result, mu0 = mu0_result)
  if (any(vapply(
    arm_results,
    function(result) is.null(result$nuisance_fit_diagnostics),
    logical(1L)
  ))) {
    return(numeric(0))
  }

  diagnostics <- do.call(rbind, lapply(
    arm_results,
    function(result) as.data.frame(result$nuisance_fit_diagnostics)
  ))
  count_columns <- c(
    "initial_outcome_degenerate",
    "target_only_outcome_degenerate",
    "initial_dr_nonconverged",
    "calibrated_dr_nonconverged",
    "calibrated_outcome_nonconverged",
    "initial_dr_cv_invalid_fold_fits",
    "initial_dr_cv_invalid_lambdas",
    "initial_dr_cv_path_tail_skipped_fold_fits",
    "calibrated_dr_cv_invalid_fold_fits",
    "calibrated_dr_cv_invalid_lambdas",
    "calibrated_dr_cv_path_tail_skipped_fold_fits",
    "calibrated_outcome_cv_invalid_fold_fits",
    "calibrated_outcome_cv_invalid_lambdas",
    "calibrated_outcome_cv_path_tail_skipped_fold_fits",
    "initial_dr_line_search_failures",
    "calibrated_dr_line_search_failures",
    "calibrated_outcome_line_search_failures",
    "initial_dr_support_floor_applied",
    "calibrated_dr_support_floor_applied"
  )
  iteration_columns <- c(
    "max_initial_dr_iterations",
    "max_calibrated_dr_iterations",
    "max_calibrated_outcome_iterations"
  )
  continuous_max_columns <- c(
    "max_initial_dr_update_ratio",
    "max_initial_dr_abs_coefficient",
    "max_calibrated_dr_update_ratio",
    "max_calibrated_dr_abs_coefficient",
    "max_calibrated_outcome_update_ratio",
    "max_calibrated_outcome_abs_coefficient",
    "max_initial_dr_support_floor",
    "max_calibrated_dr_support_floor"
  )
  missing_columns <- setdiff(
    c(count_columns, iteration_columns, continuous_max_columns),
    names(diagnostics)
  )
  if (length(missing_columns) > 0L) {
    stop(sprintf(
      ".summarize_tate_nuisance_diagnostics: missing diagnostic column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ), call. = FALSE)
  }

  counts <- vapply(
    diagnostics[count_columns],
    function(values) sum(as.numeric(values), na.rm = TRUE),
    numeric(1L)
  )
  arm_counts <- unlist(lapply(names(arm_results), function(arm_name) {
    arm_diagnostics <- as.data.frame(
      arm_results[[arm_name]]$nuisance_fit_diagnostics
    )
    values <- vapply(
      arm_diagnostics[count_columns],
      function(column_values) sum(as.numeric(column_values), na.rm = TRUE),
      numeric(1L)
    )
    names(values) <- paste0("face_", arm_name, "_", names(values))
    values
  }), use.names = TRUE)
  maxima <- vapply(
    diagnostics[iteration_columns],
    function(values) .max_iterations_or_na(values),
    integer(1L)
  )
  continuous_maxima <- vapply(
    continuous_max_columns,
    function(column) {
      if (column %in% names(diagnostics)) {
        .max_numeric_or_na(diagnostics[[column]])
      } else {
        NA_real_
      }
    },
    numeric(1L)
  )
  names(counts) <- paste0("face_", names(counts))
  names(maxima) <- paste0("face_", names(maxima))
  names(continuous_maxima) <- paste0("face_", names(continuous_maxima))
  c(counts, arm_counts, maxima, continuous_maxima)
}

.count_direct_tate_degenerate_fits <- function(direct_tate_results) {
  if (length(direct_tate_results) == 0L) {
    return(NA_integer_)
  }
  diagnostic_columns <- c(
    "initial_outcome_degenerate", "target_only_outcome_degenerate"
  )
  counts <- vapply(direct_tate_results, function(result) {
    diagnostics <- result$nuisance_fit_diagnostics
    if (is.null(diagnostics)) {
      stop(
        "TATE result is missing nuisance_fit_diagnostics.",
        call. = FALSE
      )
    }
    missing_columns <- setdiff(diagnostic_columns, names(diagnostics))
    if (length(missing_columns) > 0L) {
      stop(
        "TATE nuisance diagnostics are missing: ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    values <- unlist(diagnostics[diagnostic_columns], use.names = FALSE)
    values <- suppressWarnings(as.numeric(values))
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
      stop(
        "TATE degenerate-fit diagnostics must be finite and nonnegative.",
        call. = FALSE
      )
    }
    sum(values)
  }, numeric(1L))
  as.integer(sum(counts))
}

.named_source_diagnostics <- function(result, prefix = "") {
  required <- c(
    "weights", "source_estimates", "target_only", "fold_weights",
    "fold_wald_statistics", "fold_penalty_coefficients"
  )
  missing <- required[vapply(required, function(field) {
    is.null(result[[field]])
  }, logical(1L))]
  if (length(missing) > 0L || is.null(result$target_only$estimate)) {
    stop(
      ".named_source_diagnostics: missing aggregation field(s): ",
      paste(c(missing, if (is.null(result$target_only$estimate)) {
        "target_only$estimate"
      }), collapse = ", "),
      call. = FALSE
    )
  }
  K <- length(result$weights)
  source_names <- names(result$weights)
  if (K == 0L) return(numeric(0L))
  if (is.null(source_names) || anyNA(source_names) ||
      any(!nzchar(source_names)) || anyDuplicated(source_names)) {
    stop(
      ".named_source_diagnostics: nonempty weights require unique source names.",
      call. = FALSE
    )
  }
  if (!identical(names(result$source_estimates), source_names)) {
    stop(
      ".named_source_diagnostics: source estimate names must match weight names.",
      call. = FALSE
    )
  }
  scalar_values <- c(
    result$weights, result$source_estimates,
    target_estimate = result$target_only$estimate
  )
  if (!is.numeric(scalar_values) || any(!is.finite(scalar_values))) {
    stop(
      ".named_source_diagnostics: estimates and weights must be finite numeric values.",
      call. = FALSE
    )
  }
  source_keys <- gsub("[^A-Za-z0-9]+", "_", source_names)
  if (any(!nzchar(source_keys)) || anyDuplicated(source_keys)) {
    stop(
      ".named_source_diagnostics: source names produce ambiguous diagnostic keys.",
      call. = FALSE
    )
  }
  fold_weights <- as.matrix(result$fold_weights)
  fold_wald <- as.matrix(result$fold_wald_statistics)
  fold_penalty <- as.matrix(result$fold_penalty_coefficients)
  expected_fold_dimension <- dim(fold_weights)
  if (length(expected_fold_dimension) != 2L ||
      expected_fold_dimension[[1L]] < 1L ||
      expected_fold_dimension[[2L]] != K ||
      !identical(dim(fold_wald), expected_fold_dimension) ||
      !identical(dim(fold_penalty), expected_fold_dimension) ||
      any(!is.finite(fold_weights)) || any(!is.finite(fold_wald)) ||
      any(!is.finite(fold_penalty))) {
    stop(".named_source_diagnostics: source diagnostic dimensions disagree.",
         call. = FALSE)
  }
  matrix_names <- list(
    weights = colnames(fold_weights),
    wald = colnames(fold_wald),
    penalty = colnames(fold_penalty)
  )
  if (any(vapply(matrix_names, function(values) {
    !is.null(values) && !identical(values, source_names)
  }, logical(1L)))) {
    stop(
      ".named_source_diagnostics: matrix source names must match weight names.",
      call. = FALSE
    )
  }

  phase1 <- result$intermediates
  phase1b <- result$intermediates$phase1b
  required_phase1 <- c("source_estimates_matrix", "target_estimates")
  required_phase1b <- c(
    "avg_source_est", "avg_target_est", "V_ot", "V_t", "V_s", "C_ot",
    "n_t", "n_s"
  )
  if (is.null(phase1) || is.null(phase1b) ||
      any(!required_phase1 %in% names(phase1)) ||
      any(!required_phase1b %in% names(phase1b))) {
    stop(
      ".named_source_diagnostics: phase1 and phase1b discrepancies are required.",
      call. = FALSE
    )
  }
  evaluation_source <- as.matrix(phase1$source_estimates_matrix)
  evaluation_target <- as.numeric(phase1$target_estimates)
  screening_source <- as.matrix(phase1b$avg_source_est)
  screening_target <- as.numeric(phase1b$avg_target_est)
  if (!identical(dim(evaluation_source), expected_fold_dimension) ||
      length(evaluation_target) != expected_fold_dimension[[1L]] ||
      !identical(dim(screening_source), expected_fold_dimension) ||
      length(screening_target) != expected_fold_dimension[[1L]] ||
      any(!is.finite(evaluation_source)) || any(!is.finite(evaluation_target)) ||
      any(!is.finite(screening_source)) || any(!is.finite(screening_target))) {
    stop(
      ".named_source_diagnostics: discrepancy dimensions disagree with folds.",
      call. = FALSE
    )
  }
  if (!identical(colnames(evaluation_source), source_names) ||
      !identical(colnames(screening_source), source_names)) {
    stop(
      ".named_source_diagnostics: discrepancy source names must match weights.",
      call. = FALSE
    )
  }
  evaluation_discrepancy <- sweep(
    evaluation_source, 1L, evaluation_target, "-"
  )
  screening_discrepancy <- sweep(
    screening_source, 1L, screening_target, "-"
  )
  screening_V_ot <- as.numeric(phase1b$V_ot)
  screening_V_t <- as.matrix(phase1b$V_t)
  screening_V_s <- as.matrix(phase1b$V_s)
  screening_C_ot <- as.matrix(phase1b$C_ot)
  screening_n_t <- as.numeric(phase1b$n_t)
  screening_n_s <- as.matrix(phase1b$n_s)
  if (length(screening_V_ot) != expected_fold_dimension[[1L]] ||
      !identical(dim(screening_V_t), expected_fold_dimension) ||
      !identical(dim(screening_V_s), expected_fold_dimension) ||
      !identical(dim(screening_C_ot), expected_fold_dimension) ||
      length(screening_n_t) != expected_fold_dimension[[1L]] ||
      !identical(dim(screening_n_s), expected_fold_dimension) ||
      any(!is.finite(screening_V_ot)) || any(!is.finite(screening_V_t)) ||
      any(!is.finite(screening_V_s)) || any(!is.finite(screening_C_ot)) ||
      any(!is.finite(screening_n_t)) || any(screening_n_t <= 0) ||
      any(!is.finite(screening_n_s)) || any(screening_n_s <= 0)) {
    stop(
      ".named_source_diagnostics: Wald component dimensions disagree with folds.",
      call. = FALSE
    )
  }
  screening_discrepancy_variance <-
    sweep(screening_V_t - 2 * screening_C_ot, 1L, screening_n_t, "/") +
    screening_V_s / screening_n_s +
    screening_V_ot / screening_n_t
  screening_discrepancy_se <- sqrt(pmax(
    screening_discrepancy_variance, VARIANCE_MIN
  ))
  reconstructed_wald <- abs(screening_discrepancy) /
    screening_discrepancy_se
  wald_identity_error <- abs(reconstructed_wald - fold_wald)
  if (any(wald_identity_error > 1e-8 * (1 + abs(fold_wald)))) {
    stop(
      ".named_source_diagnostics: saved Wald statistics do not match phase1b components.",
      call. = FALSE
    )
  }

  fold_included <- NULL
  if (!is.null(result$fold_source_included)) {
    fold_included <- as.matrix(result$fold_source_included)
    if (!identical(dim(fold_included), expected_fold_dimension) ||
        anyNA(fold_included) ||
        (!is.null(colnames(fold_included)) &&
          !identical(colnames(fold_included), source_names)) ||
        any(!fold_included %in% c(TRUE, FALSE))) {
      stop(
        ".named_source_diagnostics: inclusion and weight dimensions disagree.",
        call. = FALSE
      )
    }
  }
  diagnostics <- numeric(0L)
  for (j in seq_len(K)) {
    key <- paste0(prefix, "source_", source_keys[[j]], "_")
    source_values <- c(
      target_estimate = unname(result$target_only$estimate),
      source_estimate = unname(result$source_estimates[[j]]),
      weight = unname(result$weights[[j]]),
      fold_weight_mean = mean(fold_weights[, j]),
      fold_weight_sd = stats::sd(fold_weights[, j]),
      fold_weight_min = min(fold_weights[, j]),
      fold_weight_max = max(fold_weights[, j]),
      wald_mean = mean(fold_wald[, j]),
      wald_max = max(fold_wald[, j]),
      penalty_activation_fraction = mean(fold_penalty[, j] > 0)
    )
    if (!is.null(fold_included)) {
      source_values <- c(
        source_values,
        inclusion_fraction = mean(fold_included[, j])
      )
    }
    source_values <- c(
      source_values,
      screening_discrepancy_mean = mean(screening_discrepancy[, j]),
      screening_discrepancy_mean_abs =
        mean(abs(screening_discrepancy[, j])),
      screening_discrepancy_max_abs =
        max(abs(screening_discrepancy[, j])),
      evaluation_discrepancy_mean = mean(evaluation_discrepancy[, j]),
      evaluation_discrepancy_mean_abs =
        mean(abs(evaluation_discrepancy[, j])),
      evaluation_discrepancy_max_abs =
        max(abs(evaluation_discrepancy[, j])),
      screening_discrepancy_se_mean =
        mean(screening_discrepancy_se[, j]),
      wald_identity_error_max = max(wald_identity_error[, j])
    )
    fold_values <- c(fold_weights[, j], fold_wald[, j], fold_penalty[, j])
    names(fold_values) <- c(
      paste0("fold_", seq_len(nrow(fold_weights)), "_weight"),
      paste0("fold_", seq_len(nrow(fold_wald)), "_wald"),
      paste0("fold_", seq_len(nrow(fold_penalty)), "_penalty")
    )
    fold_screening_values <- screening_discrepancy[, j]
    names(fold_screening_values) <- paste0(
      "fold_", seq_len(nrow(screening_discrepancy)),
      "_screening_discrepancy"
    )
    fold_evaluation_values <- evaluation_discrepancy[, j]
    names(fold_evaluation_values) <- paste0(
      "fold_", seq_len(nrow(evaluation_discrepancy)),
      "_evaluation_discrepancy"
    )
    fold_values <- c(
      fold_values, fold_screening_values, fold_evaluation_values
    )
    fold_screening_se_values <- screening_discrepancy_se[, j]
    names(fold_screening_se_values) <- paste0(
      "fold_", seq_len(nrow(screening_discrepancy_se)),
      "_screening_discrepancy_se"
    )
    fold_values <- c(fold_values, fold_screening_se_values)
    if (!is.null(fold_included)) {
      fold_inclusion_values <- as.numeric(fold_included[, j])
      names(fold_inclusion_values) <- paste0(
        "fold_", seq_len(nrow(fold_included)), "_included"
      )
      fold_values <- c(fold_values, fold_inclusion_values)
    }
    names(source_values) <- paste0(key, names(source_values))
    names(fold_values) <- paste0(key, names(fold_values))
    diagnostics <- c(diagnostics, source_values, fold_values)
  }
  diagnostics
}

.summarize_tate_aggregation_diagnostics <- function(result) {
  required <- c(
    "weights", "fold_wald_statistics", "fold_penalty_coefficients",
    "fold_weight_optimizer_iterations", "fold_weight_psd_ridge",
    "clip_diagnostics", "se_fixed_weights", "weight_layer"
  )
  missing <- required[vapply(required, function(field) {
    is.null(result[[field]])
  }, logical(1L))]
  if (length(missing) > 0L) {
    stop(sprintf(
      ".summarize_tate_aggregation_diagnostics: missing field(s): %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  clip_total <- result$clip_diagnostics$total
  c(
    target_anchor_weight = 1 - sum(result$weights),
    mean_abs_source_weight = mean(abs(result$weights)),
    max_abs_source_weight = max(abs(result$weights)),
    max_wald_statistic = max(result$fold_wald_statistics),
    mean_wald_statistic = mean(result$fold_wald_statistics),
    penalized_source_fold_fraction =
      mean(result$fold_penalty_coefficients > 0),
    max_weight_optimizer_iterations =
      max(result$fold_weight_optimizer_iterations),
    max_weight_psd_ridge = max(result$fold_weight_psd_ridge),
    weight_psd_ridge_fold_fraction =
      mean(result$fold_weight_psd_ridge > 0),
    se_fixed_weights = result$se_fixed_weights,
    weight_layer_indirect_variance = result$weight_layer$indirect_variance,
    weight_layer_cross_term = result$weight_layer$cross_term,
    weight_layer_kink_cells = result$weight_layer$kink_cells,
    inference_logit_truncated = clip_total$logit_truncated,
    inference_logit_truncation_fraction =
      clip_total$logit_truncation_fraction,
    inference_max_abs_logit = clip_total$max_abs_logit,
    inference_safety_clip_count =
      clip_total$weight_min_clipped + clip_total$weight_max_clipped +
      clip_total$ratio_min_clipped + clip_total$ratio_max_clipped,
    .named_source_diagnostics(result)
  )
}

#' Run a single Monte Carlo simulation replicate
#'
#' Generates data, runs the requested estimators, and returns a data frame
#' of per-method results (estimate, SE, bias, coverage, CI width).  Core method
#' failures are fail-fast: they raise an error with the current simulation and
#' method context instead of being converted to partial result rows.
#'
#' @param sim_id Integer simulation replicate ID (also used as RNG seed)
#' @param n_total Integer total sample size across all sites
#' @param K Integer number of source sites
#' @param p Integer number of covariates
#' @param config configuration ("C1", "C2", "C3", or "C4"). For the
#'   FACE negative-transfer DGP, both nuisance models always use the quadratic
#'   working basis; C2 misspecifies the true outcome mechanism, C3 the true
#'   treatment mechanism, and C4 both, through the transformed covariates of
#'   \code{generate_face_data()}. The skew-normal source-to-target density
#'   ratio is not exactly log-quadratic, so even C1 is a rich working-model
#'   setting rather than literal joint parametric correctness. For the
#'   separate RoCE DGP, the analogous C1--C4 bases are exact or misspecified as
#'   documented by \code{generate_simulation_data()}.
#' @param methods vector of methods to run
#' @param verbose Logical. Print progress messages via \code{log_info}.
#' @param n_cores_internal number of cores for internal parallelization (source sites).
#'        NULL or 1 for sequential, -1 for all cores minus 1.
#'        Note: When running simulations in parallel, set to 1 to avoid nested parallelism.
#' @param nlambda_init Integer. Number of lambda candidates for initial outcome
#'        model CV (glmnet). Lower values (e.g. 20) speed up fitting; default 100.
#' @param nuisance_lambda_rule Nuisance-model cross-validation selection rule:
#'   \code{"min"} selects the minimum validation loss and \code{"1se"}
#'   selects the largest penalty within one standard error of that minimum.
#' @param estimand_type Type of estimand for true value calculation:
#'        - "superpopulation" (default): fixed superpopulation parameter (same across all simulations)
#'        - "sample": sample-specific true value that varies across simulations
#' @param site_allocation Method for site allocation:
#'        - "model" (default): multinomial logistic model based on covariates (realistic)
#'        - "uniform": uniform random allocation (simpler, for variance validation)
#'        - "balanced": equal sample sizes per site (simplest, for baseline)
#' @param transform_type Type of covariate transformation:
#'        - "mild" (default): gentle transformations that preserve IF orthogonality
#'        - "strong": aggressive nonlinear transformations
#'        - "none": no transformation
#' @param outcome_type "binary" (logistic link) or "continuous" (identity link)
#' @param heterogeneity_type Level of outcome heterogeneity across sites for
#'   \code{dgp_type = "roce"}:
#'        - "none": homogeneous outcome models
#'        - "mild": approximately 40\% change in 2 non-zero coefficients
#'        - "strong": approximately 80\% change in 2 non-zero coefficients
#'        - "partial": only first half of source sites have mild heterogeneity.
#'   For the FACE DGP, result metadata is derived instead from
#'   \code{ate_deviation}, \code{n_deviated_sites}, and
#'   \code{effect_mod_strength}.
#' @param shift_strength Numeric multiplier for covariate shift intensity.
#' @param n_folds Integer. Number of cross-fitting folds (>= 3, default 10)
#' @param use_lambda_cache Logical. If TRUE, enable lambda caching
#'   within cross-fitting nuisance-model loops.
#' @param aggregation_lambda Positive truncated-Wald penalty multiplier.
#'   Its reciprocal is the source penalty-activation cutoff. The locked
#'   manuscript implementation uses \code{1}, corresponding to cutoff
#'   \code{1}; this value was selected without inspecting coverage in a
#'   disjoint pilot experiment.
#' @param n_bootstrap Number of multiplier-bootstrap draws used by comparison
#'   methods. The production default is 5,000.
#' @param M_tau Positive fitting-stage truncation radius.
#' @param M_tau_inference Positive or infinite inference-stage truncation
#'   radius used in source corrections and influence-function variances.
#' @param estimate_ate Logical. If TRUE, also estimate target-only ATE
#' @param return_fitted_tate Logical. If \code{TRUE}, attach the generated
#'   site-split data and fitted TATE objects in the
#'   \code{roce_simulation_artifacts} attribute. This opt-in path supports
#'   post-fit cutoff and inference-radius sensitivity analyses without changing
#'   the primary result rows or refitting high-dimensional nuisances.
#' @param rho_reuse_reference Optional internal reuse specification for a
#'   same-seed FACE-DGP rho group. It must contain the rho-zero
#'   \code{data_split}, its one-round \code{fitted_tate}, and the source names
#'   whose treated outcomes change. Strict equality checks reject changes to
#'   the target, covariates, treatment, control outcomes, or non-refitted
#'   sources. The default \code{NULL} performs a fully independent fit.
#' @param parallel_treatment_arms Logical. If \code{TRUE} and
#'   \code{estimate_ate = TRUE}, fit the treated and control RoCE nuisance
#'   pipelines concurrently. The caller should allocate approximately twice
#'   \code{n_cores_internal}; defaults to \code{FALSE}.
#' @param n_target Optional integer target-site sample size for explicit
#'   per-site allocation (FACE DGP only; paired with \code{n_source_sizes}).
#' @param n_source_sizes Optional integer vector of per-site source sample
#'   sizes (FACE DGP only); when supplied, \code{K} and the total are derived
#'   from it and \code{n_total} is ignored.
#' @param dgp_type Data-generating process, \code{"face"} or \code{"roce"}.
#' @param ate_deviation Non-negative additive source treatment-shift deviation
#'   under the FACE DGP: a mean shift for Gaussian outcomes and a log-odds
#'   shift for binary outcomes. The historical argument name is retained for
#'   compatibility; it is not generally equal to the induced marginal TATE
#'   difference for binary outcomes.
#' @param n_deviated_sites Number of leading deviated sources under the FACE DGP.
#' @param deviation_mechanism \code{"treated_arm"} (default) or
#'   \code{"both_arms"}: how the deviated sources deviate under the FACE DGP
#'   (see \code{\link{generate_face_data}}). Under \code{"both_arms"} the
#'   positive-rho reuse refits both arms of the changed sources.
#' @param effect_mod_strength Source-only treatment-effect-modification strength
#'   under the FACE DGP; zero recovers the standard DGP.
#' @param include_hard_threshold_diagnostic Logical. Also evaluate a diagnostic
#'   TATE aggregation that sets sources above the foldwise Wald cutoff to zero
#'   before variance minimization. Defaults to \code{FALSE}; this diagnostic
#'   reuses the fitted nuisances and does not replace the primary soft-penalty
#'   estimator.
#' @param include_quadratic_bias_rule Logical. Also evaluate the pre-specified
#'   smooth quadratic-bias sensitivity rule (\code{screening_rule =
#'   "quadratic_bias"}) as the \code{<method>_ate_quadratic_bias} row. Defaults
#'   to \code{TRUE}; it reuses the fitted nuisances and does not replace the
#'   primary soft-penalty estimator.
#' @param n_weight_bootstrap Number of multiplier-bootstrap draws that relearn
#'   the common source weights for the primary soft-penalty direct-TATE row.
#'   Use 0 (the default) to disable this diagnostic without changing existing
#'   estimates, standard errors, or RNG streams; positive values must be at
#'   least 2. The analytic standard error remains in \code{se}.
#' @return data frame with results
#' @export
run_single_simulation <- function(sim_id, n_total = 1000, K = 3, p = 4,
                                 config = "C1",
                                 methods = c("two_round_crossfit", "one_round_crossfit",
                                           "target_only", "sample_size", "inverse_variance",
                                           "federated_dr", "pooled_dr", "tilted_aipw",
                                           "oracle_dr"),
                                 verbose = TRUE, n_cores_internal = NULL,
                                 nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                                 estimand_type = "superpopulation",
                                 site_allocation = "model",
                                 transform_type = "mild",
                                 outcome_type = "binary",
                                 heterogeneity_type = "none",
                                 shift_strength = ROCE_SHIFT_STRENGTH_DEFAULT,
                                 n_folds = N_FOLDS_DEFAULT,
                                 use_lambda_cache = TRUE,
                                 aggregation_lambda = AGG_WALD_LAMBDA,
                                 n_bootstrap = BOOTSTRAP_REPLICATES_DEFAULT,
                                 M_tau = M_TAU_DEFAULT,
                                 M_tau_inference = M_TAU_INFERENCE_DEFAULT,
                                 estimate_ate = FALSE,
                                 return_fitted_tate = FALSE,
                                 rho_reuse_reference = NULL,
                                 parallel_treatment_arms = FALSE,
                                 # FACE paper DGP parameters
                                 dgp_type = "face",
                                 ate_deviation    = 0.0,
                                 n_deviated_sites = 0L,
                                 deviation_mechanism = c("treated_arm", "both_arms"),
                                 effect_mod_strength = 0,
                                 # Explicit per-site sample sizes (FACE DGP only)
                                 n_target         = NULL,
                                 n_source_sizes   = NULL,
                                 nuisance_lambda_rule = c("min", "1se"),
                                 include_hard_threshold_diagnostic = FALSE,
                                 include_quadratic_bias_rule = TRUE,
                                 n_weight_bootstrap = 0L) {

  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "run_single_simulation",
    arg = "nuisance_lambda_rule"
  )
  deviation_mechanism <- match.arg(deviation_mechanism)
  if (length(include_hard_threshold_diagnostic) != 1L ||
      is.na(include_hard_threshold_diagnostic) ||
      !is.logical(include_hard_threshold_diagnostic)) {
    stop(
      "include_hard_threshold_diagnostic must be TRUE or FALSE.",
      call. = FALSE
    )
  }
  if (length(include_quadratic_bias_rule) != 1L ||
      is.na(include_quadratic_bias_rule) ||
      !is.logical(include_quadratic_bias_rule)) {
    stop("include_quadratic_bias_rule must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(aggregation_lambda) != 1L ||
      !is.finite(aggregation_lambda) || aggregation_lambda <= 0) {
    stop("aggregation_lambda must be one finite positive number.",
         call. = FALSE)
  }
  n_bootstrap <- .validate_bootstrap_replicates(
    n_bootstrap, "run_single_simulation"
  )
  n_weight_bootstrap <- .validate_weight_bootstrap_replicates(
    n_weight_bootstrap, "run_single_simulation"
  )
  validate_truncation_parameters(M_tau, M_tau_inference)
  if (length(parallel_treatment_arms) != 1L ||
      is.na(parallel_treatment_arms) ||
      !is.logical(parallel_treatment_arms)) {
    stop("parallel_treatment_arms must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(return_fitted_tate) != 1L || is.na(return_fitted_tate) ||
      !is.logical(return_fitted_tate)) {
    stop("return_fitted_tate must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(return_fitted_tate) && !isTRUE(estimate_ate)) {
    stop(
      "return_fitted_tate = TRUE requires estimate_ate = TRUE.",
      call. = FALSE
    )
  }
  if (n_weight_bootstrap > 0L &&
      (!isTRUE(estimate_ate) || !"one_round_crossfit" %in% methods)) {
    stop(
      paste0(
        "n_weight_bootstrap > 0 requires estimate_ate = TRUE and the ",
        "primary one_round_crossfit method."
      ),
      call. = FALSE
    )
  }
  if (!is.null(rho_reuse_reference)) {
    required_reuse_fields <- c(
      "data_split", "fitted_tate", "changed_sources"
    )
    missing_reuse_fields <- setdiff(
      required_reuse_fields, names(rho_reuse_reference)
    )
    if (!is.list(rho_reuse_reference) ||
        length(missing_reuse_fields) > 0L) {
      stop(
        "rho_reuse_reference must contain data_split, fitted_tate, and changed_sources.",
        call. = FALSE
      )
    }
    if (!isTRUE(estimate_ate) || !identical(dgp_type, "face") ||
        ate_deviation <= 0 || as.integer(n_deviated_sites) < 1L ||
        !"one_round_crossfit" %in% methods) {
      stop(
        paste0(
          "rho_reuse_reference is restricted to positive-rho FACE-DGP ",
          "TATE runs that request one_round_crossfit."
        ),
        call. = FALSE
      )
    }
  }

  # Resolve explicit per-site sample sizes up front so logging and result
  # bookkeeping use a concrete (n_total, K) consistent with the allocation.
  if (dgp_type == "face" && (!is.null(n_source_sizes) || !is.null(n_target))) {
    site_sizes <- resolve_face_site_sizes(n_total, n_target, n_source_sizes, K)
    n_total    <- site_sizes$n_total
    K          <- site_sizes$K
  }
  result_heterogeneity_type <- if (identical(dgp_type, "face")) {
    .face_heterogeneity_type(
      ate_deviation, n_deviated_sites, effect_mod_strength, deviation_mechanism
    )
  } else {
    heterogeneity_type
  }

  # Track timing for each stage
  sim_start_time <- Sys.time()
  stage_times <- list()

  # Reset graceful-degradation counters for this simulation (see fit_glmnet_cv:
  # single-class outcome folds fall back to the constant mean and are counted).
  .reset_fit_diagnostics()

  log_info(verbose, "\n[%s] Simulation %d (config=%s, dgp=%s, heterogeneity=%s)\n",
           format(Sys.time(), "%H:%M:%S"), sim_id, config, dgp_type,
           result_heterogeneity_type)

  # Generate data
  data_gen_start <- Sys.time()
  set.seed(sim_id)
  data <- generate_simulation_data(n_total, K, p, config,
                                   estimand_type    = estimand_type,
                                   site_allocation  = site_allocation,
                                   transform_type   = transform_type,
                                   outcome_type     = outcome_type,
                                   heterogeneity_type = heterogeneity_type,
                                   shift_strength   = shift_strength,
                                   dgp_type         = dgp_type,
                                   ate_deviation    = ate_deviation,
                                   n_deviated_sites = as.integer(n_deviated_sites),
                                   deviation_mechanism = deviation_mechanism,
                                   effect_mod_strength = effect_mod_strength,
                                   n_target         = n_target,
                                   n_source_sizes   = n_source_sizes,
                                   warn_ignored     = FALSE)
  data_split <- split_data_by_site(data)
  stage_times$data_gen <- as.numeric(difftime(Sys.time(), data_gen_start, units = "secs"))

  log_info(verbose, "    Data generated in %.2fs (n_t=%d, n_s=%s)\n",
           stage_times$data_gen, data_split$t$n,
           paste(sapply(setdiff(names(data_split), "t"), function(s) data_split[[s]]$n), collapse=","))

  # Get true potential outcome from generated data
  # The interpretation depends on estimand_type:
  #   - "sample": sample-specific E_n[Y(1)] (varies across simulations)
  #   - "superpopulation": fixed E[Y(1)] (same for all simulations)
  true_potential_outcome <- data$mu1_true

  # Find target sample indices (for logging)
  target_idx <- which(data$R == "t")

  if (estimand_type == "superpopulation") {
    log_info(verbose, "  True potential outcome (superpop): %.4f, Sample-specific: %.4f\n",
             true_potential_outcome, data$mu1_realized)
  } else {
    log_info(verbose, "  True potential outcome (sample): %.4f\n", true_potential_outcome)
  }
  log_info(verbose, "  Target sample size: %d, proportion: %.3f\n",
           length(target_idx), length(target_idx) / n_total)

  # Initialize results as a list (avoid O(n²) rbind-in-loop)
  results_list <- list()
  crossfit_mu1_results <- list()
  crossfit_mu0_results <- list()
  direct_tate_results <- list()

  # Precompute fold partition once and reuse across one/two-round and ATE reruns
  precomputed_folds <- build_crossfit_folds(data_split, n_folds)
  target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
  target_only_fit_cache <- new.env(hash = TRUE, parent = emptyenv())

  # Map the realized DGP outcome type to the GLM family for cross-fitting.  Use
  # the field returned by the generator so future DGP-specific defaults cannot
  # make the fitted family diverge from the generated outcome.
  family <- switch(data$outcome_type,
    "binary"     = "binomial",
    "continuous" = "gaussian",
    tolower(data$outcome_type)
  )
  if (!family %in% VALID_GLM_FAMILIES) {
    stop(sprintf("Unsupported GLM family '%s'. Supported families: %s",
                 family, paste(VALID_GLM_FAMILIES, collapse = ", ")))
  }
  cell_diagnostics <- .binary_cell_diagnostics(data_split, family)

  # Run one-round cross-fitting algorithm (if requested)
  if ("one_round_crossfit" %in% methods) {
    one_round_start <- Sys.time()
    log_info(verbose, "    Running one-round cross-fitting (nlambda=%d)...\n", nlambda_init)
    # Route both scheduling topologies through the same TATE wrapper.
    # Besides assembling the joint contrast, run_tate_crossfit gives each arm
    # a deterministic work-unit RNG stream.  Bypassing it in sequential mode
    # would make a resource choice change the fitted estimator.
    if (isTRUE(estimate_ate)) {
      one_round_tate_res <- if (is.null(rho_reuse_reference)) {
        run_tate_crossfit(
          data_split = data_split,
          n_folds = n_folds,
          communication_mode = "one_round",
          lambda_selection = aggregation_lambda,
          verbose = FALSE,
          n_cores = n_cores_internal,
          nlambda_init = nlambda_init,
          family = family,
          M_tau = M_tau,
          M_tau_inference = M_tau_inference,
          use_lambda_cache = use_lambda_cache,
          precomputed_folds = precomputed_folds,
          target_only_ps_cache = target_only_ps_cache,
          target_only_fit_cache = target_only_fit_cache,
          nuisance_lambda_rule = nuisance_lambda_rule,
          parallel_arms = parallel_treatment_arms
        )
      } else {
        reference_tate <- rho_reuse_reference$fitted_tate
        reference_nlambda <- as.integer(
          reference_tate$arm_results$mu1$fold_results[[1L]]$
            source_results[[1L]]$nuisance_nlambda
        )
        if (!identical(reference_nlambda, as.integer(nlambda_init)) ||
            !identical(reference_tate$nuisance_lambda_rule,
                       nuisance_lambda_rule) ||
            !isTRUE(all.equal(reference_tate$M_tau, M_tau)) ||
            !isTRUE(all.equal(reference_tate$M_tau_inference,
                             M_tau_inference))) {
          stop(
            paste0(
              "rho_reuse_reference nuisance grid, CV rule, and truncation ",
              "parameters must match the current run."
            ),
            call. = FALSE
          )
        }
        # The shared-shift mechanism changes both arms of the deviated
        # sources, so their control-arm fits are refitted as well.
        .reuse_one_round_tate_across_rho(
          reference_data_split = rho_reuse_reference$data_split,
          data_split = data_split,
          fitted_tate = reference_tate,
          changed_sources = rho_reuse_reference$changed_sources,
          lambda_selection = aggregation_lambda,
          refit_control_arm = identical(deviation_mechanism, "both_arms"),
          verbose = FALSE,
          n_cores = n_cores_internal
        )
      }
      one_round_cf_res <- one_round_tate_res$arm_results$mu1
      crossfit_mu0_results[["one_round_crossfit"]] <-
        one_round_tate_res$arm_results$mu0
      direct_tate_results[["one_round_crossfit"]] <- one_round_tate_res
    } else {
      one_round_cf_res <- run_crossfit(
        data_split, n_folds = n_folds,
        communication_mode = "one_round",
        lambda_selection = aggregation_lambda,
        verbose = FALSE, n_cores = n_cores_internal,
        nlambda_init = nlambda_init,
        family = family,
        M_tau = M_tau,
        M_tau_inference = M_tau_inference,
        use_lambda_cache = use_lambda_cache,
        precomputed_folds = precomputed_folds,
        target_only_ps_cache = target_only_ps_cache,
        target_only_fit_cache = target_only_fit_cache,
        nuisance_lambda_rule = nuisance_lambda_rule
      )
    }
    stage_times$one_round <- as.numeric(difftime(Sys.time(), one_round_start, units = "secs"))

    log_info(verbose, "    [OK] One-round done in %.2fs: est=%.4f, bias=%.4f, SE=%.4f\n",
             stage_times$one_round,
             one_round_cf_res$estimate, one_round_cf_res$estimate - true_potential_outcome,
             one_round_cf_res$se)
    crossfit_mu1_results[["one_round_crossfit"]] <- one_round_cf_res

    results_list[[length(results_list) + 1]] <- data.frame(
      sim_id = sim_id,
      method = "one_round_crossfit",
      estimate = one_round_cf_res$estimate,
      se = one_round_cf_res$se,
      bias = one_round_cf_res$estimate - true_potential_outcome,
      coverage = (true_potential_outcome >= one_round_cf_res$ci_lower) &
                (true_potential_outcome <= one_round_cf_res$ci_upper),
      ci_width = one_round_cf_res$ci_upper - one_round_cf_res$ci_lower,
      n_total = n_total,
      K = K,
      p = p,
      config = config,
      heterogeneity_type = heterogeneity_type,
      estimand_type = estimand_type,
      stringsAsFactors = FALSE
    )
  }

  # Run two-round cross-fitting algorithm
  # Note: n_folds >= 3 is required for proper two-level cross-fitting calibration
  if ("two_round_crossfit" %in% methods) {
    two_round_start <- Sys.time()
    log_info(verbose, "    Running two-round cross-fitting (nlambda=%d)...\n", nlambda_init)
    if (isTRUE(estimate_ate)) {
      two_round_tate_res <- run_tate_crossfit(
        data_split = data_split,
        n_folds = n_folds,
        communication_mode = "two_round",
        lambda_selection = aggregation_lambda,
        verbose = FALSE,
        n_cores = n_cores_internal,
        nlambda_init = nlambda_init,
        family = family,
        M_tau = M_tau,
        M_tau_inference = M_tau_inference,
        use_lambda_cache = use_lambda_cache,
        precomputed_folds = precomputed_folds,
        target_only_ps_cache = target_only_ps_cache,
        target_only_fit_cache = target_only_fit_cache,
        nuisance_lambda_rule = nuisance_lambda_rule,
        parallel_arms = parallel_treatment_arms
      )
      two_round_cf_res <- two_round_tate_res$arm_results$mu1
      crossfit_mu0_results[["two_round_crossfit"]] <-
        two_round_tate_res$arm_results$mu0
      direct_tate_results[["two_round_crossfit"]] <- two_round_tate_res
    } else {
      two_round_cf_res <- run_crossfit(
        data_split, n_folds = n_folds,
        communication_mode = "two_round",
        lambda_selection = aggregation_lambda,
        verbose = FALSE, n_cores = n_cores_internal,
        nlambda_init = nlambda_init,
        family = family,
        M_tau = M_tau,
        M_tau_inference = M_tau_inference,
        use_lambda_cache = use_lambda_cache,
        precomputed_folds = precomputed_folds,
        target_only_ps_cache = target_only_ps_cache,
        target_only_fit_cache = target_only_fit_cache,
        nuisance_lambda_rule = nuisance_lambda_rule
      )
    }
    stage_times$two_round <- as.numeric(difftime(Sys.time(), two_round_start, units = "secs"))

    log_info(verbose, "    [OK] Two-round done in %.2fs: est=%.4f, bias=%.4f, SE=%.4f\n",
             stage_times$two_round,
             two_round_cf_res$estimate, two_round_cf_res$estimate - true_potential_outcome,
             two_round_cf_res$se)
    crossfit_mu1_results[["two_round_crossfit"]] <- two_round_cf_res

    results_list[[length(results_list) + 1]] <- data.frame(
      sim_id = sim_id,
      method = "two_round_crossfit",
      estimate = two_round_cf_res$estimate,
      se = two_round_cf_res$se,
      bias = two_round_cf_res$estimate - true_potential_outcome,
      coverage = (true_potential_outcome >= two_round_cf_res$ci_lower) &
                (true_potential_outcome <= two_round_cf_res$ci_upper),
      ci_width = two_round_cf_res$ci_upper - two_round_cf_res$ci_lower,
      n_total = n_total,
      K = K,
      p = p,
      config = config,
      heterogeneity_type = heterogeneity_type,
      estimand_type = estimand_type,
      stringsAsFactors = FALSE
    )
  }

  # Run comparison methods
  comp_start <- Sys.time()
  comparison_method_names <- c(
    "target_only", "sample_size", "inverse_variance",
    "federated_dr", "pooled_dr", "tilted_aipw"
  )
  requested_comparisons <- intersect(methods, comparison_method_names)
  comparison_dr_methods <- c("federated_dr", "pooled_dr")
  comparison_dr_weights <- if (
      any(requested_comparisons %in% comparison_dr_methods)
  ) {
    .resolve_dr_weights_by_site(
      data_split = data_split,
      n_cores = n_cores_internal,
      caller = "run_single_simulation"
    )
  } else {
    NULL
  }
  log_info(
    verbose, "    Running comparison methods: %s\n",
    if (length(requested_comparisons) > 0L) {
      paste(requested_comparisons, collapse = ", ")
    } else {
      "none"
    }
  )
  comparison_results <- if (length(requested_comparisons) > 0L) {
    run_all_comparisons(
      data_split,
      family = family,
      n_folds = n_folds,
      include_tilted = "tilted_aipw" %in% requested_comparisons,
      methods = requested_comparisons,
      n_cores = n_cores_internal,
      n_bootstrap = n_bootstrap,
      dr_weights_by_site = comparison_dr_weights
    )
  } else {
    list()
  }
  stage_times$comparison <- as.numeric(difftime(Sys.time(), comp_start, units = "secs"))

  for (method_name in comparison_method_names) {
    if (method_name %in% methods) {
      method_res <- .require_method_result(comparison_results, method_name)

      # Calculate confidence interval
      ci_lower <- method_res$estimate - Z_ALPHA_05 * method_res$se
      ci_upper <- method_res$estimate + Z_ALPHA_05 * method_res$se

      comparison_row <- data.frame(
        sim_id = sim_id,
        method = method_name,
        estimate = method_res$estimate,
        se = method_res$se,
        bias = method_res$estimate - true_potential_outcome,
        coverage = (true_potential_outcome >= ci_lower) & (true_potential_outcome <= ci_upper),
        ci_width = ci_upper - ci_lower,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type,
        stringsAsFactors = FALSE
      )
      comparison_diagnostics <- .comparison_variance_diagnostics(method_res)
      for (diagnostic_name in names(comparison_diagnostics)) {
        comparison_row[[diagnostic_name]] <-
          comparison_diagnostics[[diagnostic_name]]
      }
      results_list[[length(results_list) + 1]] <- comparison_row
    }
  }

  # Run oracle DR estimator (uses known true parameters)
  # Not available for the FACE paper DGP (gamma_params and alpha1_true are NULL).
  if ("oracle_dr" %in% methods) {
    oracle_start <- Sys.time()
    if (is.null(data$gamma_params) || is.null(data$alpha1_true)) {
      log_info(verbose, "    [--] Oracle DR skipped: true parameters not available for dgp_type='%s'\n",
               data$dgp_type %||% dgp_type)
    } else {
      log_info(verbose, "    Running oracle DR estimator...\n")
      target_propensity_true <- NULL
      if (!is.null(data$p_treat_true)) {
        target_propensity_true <- data$p_treat_true[target_idx]
      }
      oracle_res <- estimate_oracle_dr(data_split, data$alpha1_true, data$gamma_params,
                                       outcome_type = data$outcome_type,
                                       target_propensity_true = target_propensity_true)
      stage_times$oracle <- as.numeric(difftime(Sys.time(), oracle_start, units = "secs"))

      ci_lower <- oracle_res$estimate - Z_ALPHA_05 * oracle_res$se
      ci_upper <- oracle_res$estimate + Z_ALPHA_05 * oracle_res$se

      log_info(verbose, "    [OK] Oracle DR done in %.2fs: est=%.4f, bias=%.4f\n",
               stage_times$oracle, oracle_res$estimate,
               oracle_res$estimate - true_potential_outcome)

      results_list[[length(results_list) + 1]] <- data.frame(
        sim_id = sim_id,
        method = "oracle_dr",
        estimate = oracle_res$estimate,
        se = oracle_res$se,
        bias = oracle_res$estimate - true_potential_outcome,
        coverage = (true_potential_outcome >= ci_lower) & (true_potential_outcome <= ci_upper),
        ci_width = ci_upper - ci_lower,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type,
        stringsAsFactors = FALSE
      )
    }  # end else (gamma_params guard)
  }

  # TATE estimation. The `_ate` suffix is retained in method labels for
  # compatibility with existing simulation CSVs and plotting scripts.
  if (isTRUE(estimate_ate)) {
    log_info(verbose, "    Assembling TATE and comparison contrasts...\n")
    tate_start <- Sys.time()
    tate_truth <- data$mu1_true - data$mu0_true
    target_only_tate <- NULL

    for (crossfit_method in c("two_round_crossfit", "one_round_crossfit")) {
      if (crossfit_method %in% methods) {
        communication_mode <- if (crossfit_method == "two_round_crossfit") {
          "two_round"
        } else {
          "one_round"
        }
        mu0_res <- crossfit_mu0_results[[crossfit_method]]
        if (is.null(mu0_res)) {
          mu0_res <- run_crossfit(
            data_split, n_folds = n_folds,
            communication_mode = communication_mode,
            lambda_selection = aggregation_lambda,
            verbose = FALSE, n_cores = n_cores_internal,
            nlambda_init = nlambda_init, family = family, A_val = 0L,
            M_tau = M_tau, M_tau_inference = M_tau_inference,
            use_lambda_cache = use_lambda_cache,
            precomputed_folds = precomputed_folds,
            target_only_ps_cache = target_only_ps_cache,
            target_only_fit_cache = target_only_fit_cache,
            nuisance_lambda_rule = nuisance_lambda_rule
          )
        }

        mu1_cf <- crossfit_mu1_results[[crossfit_method]]
        if (is.null(mu1_cf) || is.null(mu1_cf$all_phi_agg) ||
            is.null(mu0_res$all_phi_agg)) {
          .abort_with_context(
            "arm-specific pseudo-values for %s are unavailable; cannot construct TATE.",
            crossfit_method
          )
        }

        group_sizes <- vapply(data_split, function(site) site$n, integer(1L))
        armwise_phi_tau <- mu1_cf$all_phi_agg - mu0_res$all_phi_agg
        armwise_tate_estimate <- mu1_cf$estimate - mu0_res$estimate
        armwise_tate_variance <- .multisite_pseudovalue_variance(
          armwise_phi_tau, group_sizes
        )
        armwise_tate_se <- sqrt(armwise_tate_variance)

        # Preserve the previous arm-wise construction as an explicitly labeled
        # diagnostic. It now uses the correct within-site-centered joint IF
        # variance rather than globally squaring uncentered pseudo-values.
        armwise_tate_row <- .make_simulation_result_row(
          sim_id = sim_id,
          method = paste0(crossfit_method, "_ate_armwise"),
          estimate = armwise_tate_estimate,
          se = armwise_tate_se,
          truth = tate_truth,
          n_total = n_total,
          K = K,
          p = p,
          config = config,
          heterogeneity_type = heterogeneity_type,
          estimand_type = estimand_type
        )
        armwise_source_diagnostics <- c(
          .named_source_diagnostics(mu1_cf, prefix = "mu1_"),
          .named_source_diagnostics(mu0_res, prefix = "mu0_")
        )
        for (diagnostic_name in names(armwise_source_diagnostics)) {
          armwise_tate_row[[diagnostic_name]] <-
            armwise_source_diagnostics[[diagnostic_name]]
        }
        results_list[[length(results_list) + 1L]] <- armwise_tate_row

        direct_tate_res <- direct_tate_results[[crossfit_method]]
        if (is.null(direct_tate_res)) {
          direct_tate_res <- calculate_tate_crossfit_aggregation(
            data_split = data_split,
            mu1_result = mu1_cf,
            mu0_result = mu0_res,
            lambda_selection = aggregation_lambda,
            verbose = FALSE
          )
        }
        direct_tate_res$communication_mode <- communication_mode
        direct_tate_res$parallel_arms <- parallel_treatment_arms
        direct_tate_res$nuisance_lambda_rule <-
          mu1_cf$nuisance_lambda_rule %||% "min"
        direct_tate_res$family <- family
        direct_tate_res$M_tau <- M_tau
        direct_tate_res$M_tau_inference <- M_tau_inference
        direct_tate_res$aggregation_lambda_selection <- aggregation_lambda
        if (is.null(direct_tate_res$nuisance_fit_diagnostics)) {
          direct_tate_res$nuisance_fit_diagnostics <- rbind(
            transform(mu1_cf$nuisance_fit_diagnostics, A_val = 1L),
            transform(mu0_res$nuisance_fit_diagnostics, A_val = 0L)
          )
          rownames(direct_tate_res$nuisance_fit_diagnostics) <- NULL
        }
        # Retain the assembled TATE object for opt-in post-fit
        # sensitivity reuse even when the two arms were fitted sequentially.
        direct_tate_results[[crossfit_method]] <- direct_tate_res
        if (is.null(target_only_tate)) {
          target_only_tate <- direct_tate_res$target_only
        }

        direct_tate_row <- .make_simulation_result_row(
          sim_id = sim_id,
          method = paste0(crossfit_method, "_ate"),
          estimate = direct_tate_res$estimate,
          se = direct_tate_res$se,
          truth = tate_truth,
          n_total = n_total,
          K = K,
          p = p,
          config = config,
          heterogeneity_type = heterogeneity_type,
          estimand_type = estimand_type
        )
        aggregation_diagnostics <-
          .summarize_tate_aggregation_diagnostics(direct_tate_res)
        for (diagnostic_name in names(aggregation_diagnostics)) {
          direct_tate_row[[diagnostic_name]] <-
            aggregation_diagnostics[[diagnostic_name]]
        }
        crossfit_diagnostics <- c(
          .summarize_tate_crossfit_timing(mu1_cf, mu0_res),
          .summarize_tate_nuisance_diagnostics(mu1_cf, mu0_res)
        )
        for (diagnostic_name in names(crossfit_diagnostics)) {
          direct_tate_row[[diagnostic_name]] <-
            crossfit_diagnostics[[diagnostic_name]]
        }
        if (n_weight_bootstrap > 0L &&
            identical(crossfit_method, "one_round_crossfit")) {
          weight_bootstrap <- .run_tate_weight_bootstrap(
            tate_result = direct_tate_res,
            n_weight_bootstrap = n_weight_bootstrap,
            sim_id = sim_id,
            screening_rule = "soft_penalty"
          )
          for (diagnostic_name in names(weight_bootstrap)) {
            direct_tate_row[[diagnostic_name]] <-
              weight_bootstrap[[diagnostic_name]]
          }
        }
        results_list[[length(results_list) + 1L]] <- direct_tate_row

        # The smooth quadratic-bias rule is the pre-specified sensitivity
        # estimator; it reuses the same arm-specific nuisance fits. The gate
        # lets a production run recover from a failure in this path without
        # touching the primary row.
        if (isTRUE(include_quadratic_bias_rule)) {
          quadratic_tate_res <- calculate_tate_crossfit_aggregation(
            data_split = data_split,
            mu1_result = mu1_cf,
            mu0_result = mu0_res,
            lambda_selection = aggregation_lambda,
            verbose = FALSE,
            screening_rule = "quadratic_bias"
          )
          quadratic_tate_row <- .make_simulation_result_row(
            sim_id = sim_id,
            method = paste0(crossfit_method, "_ate_quadratic_bias"),
            estimate = quadratic_tate_res$estimate,
            se = quadratic_tate_res$se,
            truth = tate_truth,
            n_total = n_total,
            K = K,
            p = p,
            config = config,
            heterogeneity_type = heterogeneity_type,
            estimand_type = estimand_type
          )
          quadratic_diagnostics <- c(
            .summarize_tate_aggregation_diagnostics(quadratic_tate_res),
            crossfit_diagnostics
          )
          for (diagnostic_name in names(quadratic_diagnostics)) {
            quadratic_tate_row[[diagnostic_name]] <-
              quadratic_diagnostics[[diagnostic_name]]
          }
          results_list[[length(results_list) + 1L]] <- quadratic_tate_row
        }

        if (isTRUE(include_hard_threshold_diagnostic)) {
          hard_tate_res <- calculate_tate_crossfit_aggregation(
            data_split = data_split,
            mu1_result = mu1_cf,
            mu0_result = mu0_res,
            lambda_selection = aggregation_lambda,
            verbose = FALSE,
            screening_rule = "hard_threshold"
          )
          hard_tate_row <- .make_simulation_result_row(
            sim_id = sim_id,
            method = paste0(crossfit_method, "_ate_hard_threshold"),
            estimate = hard_tate_res$estimate,
            se = hard_tate_res$se,
            truth = tate_truth,
            n_total = n_total,
            K = K,
            p = p,
            config = config,
            heterogeneity_type = heterogeneity_type,
            estimand_type = estimand_type
          )
          hard_diagnostics <-
            .summarize_tate_aggregation_diagnostics(hard_tate_res)
          for (diagnostic_name in names(hard_diagnostics)) {
            hard_tate_row[[diagnostic_name]] <-
              hard_diagnostics[[diagnostic_name]]
          }
          # The hard-screen diagnostic reuses the exact same arm-specific
          # nuisance fits as the primary common-weight TATE estimator.  Copy
          # their diagnostics so fail-closed replicate QC can audit this row
          # independently instead of mistaking shared fits for missing fits.
          for (diagnostic_name in names(crossfit_diagnostics)) {
            hard_tate_row[[diagnostic_name]] <-
              crossfit_diagnostics[[diagnostic_name]]
          }
          results_list[[length(results_list) + 1L]] <- hard_tate_row
        }
        log_info(
          verbose,
          "    [OK] %s TATE: est=%.4f, true=%.4f, bias=%.4f; arm-wise diagnostic=%.4f\n",
          crossfit_method, direct_tate_res$estimate, tate_truth,
          direct_tate_res$estimate - tate_truth, armwise_tate_estimate
        )
      }
    }

    # The direct RoCE optimization uses target-only TATE as its eta = 0
    # reference. Reusing that exact cross-fitted result keeps the baseline and
    # aggregation objective on the same folds and nuisance fits.
    if ("target_only" %in% methods) {
      if (is.null(target_only_tate)) {
        .abort_with_context(
          paste0(
            "target-only TATE requires at least one requested RoCE cross-fit ",
            "method so the eta=0 reference uses the same folds."
          )
        )
      }
      results_list[[length(results_list) + 1L]] <- .make_simulation_result_row(
        sim_id = sim_id,
        method = "target_only_ate",
        estimate = target_only_tate$estimate,
        se = target_only_tate$se,
        truth = tate_truth,
        n_total = n_total,
        K = K,
        p = p,
        config = config,
        heterogeneity_type = heterogeneity_type,
        estimand_type = estimand_type
      )
      log_info(
        verbose,
        "    [OK] Target-only TATE: est=%.4f, true=%.4f, bias=%.4f\n",
        target_only_tate$estimate, tate_truth,
        target_only_tate$estimate - tate_truth
      )
    }

    tate_comparison_methods <- setdiff(
      requested_comparisons, "target_only"
    )
    if (length(tate_comparison_methods) > 0L) {
      comparison_tate_results <- run_all_comparisons_tate(
        data_split = data_split,
        family = family,
        n_folds = n_folds,
        include_tilted = "tilted_aipw" %in% tate_comparison_methods,
        methods = tate_comparison_methods,
        mu1_results = comparison_results,
        n_cores = n_cores_internal,
        n_bootstrap = n_bootstrap,
        dr_weights_by_site = comparison_dr_weights
      )
      for (method_name in tate_comparison_methods) {
        method_result <- .require_method_result(
          comparison_tate_results, method_name
        )
        comparison_tate_row <- .make_simulation_result_row(
            sim_id = sim_id,
            method = paste0(method_name, "_ate"),
            estimate = method_result$estimate,
            se = method_result$se,
            truth = tate_truth,
            n_total = n_total,
            K = K,
            p = p,
            config = config,
            heterogeneity_type = heterogeneity_type,
            estimand_type = estimand_type
          )
        comparison_diagnostics <-
          .comparison_variance_diagnostics(method_result)
        for (diagnostic_name in names(comparison_diagnostics)) {
          comparison_tate_row[[diagnostic_name]] <-
            comparison_diagnostics[[diagnostic_name]]
        }
        results_list[[length(results_list) + 1L]] <- comparison_tate_row
      }
    }

    stage_times$tate <- as.numeric(difftime(Sys.time(), tate_start, units = "secs"))
  }

  # Print total simulation time summary
  sim_total_time <- as.numeric(difftime(Sys.time(), sim_start_time, units = "secs"))
  log_info(verbose, "    [OK] Comparison methods done in %.2fs\n", stage_times$comparison)
  log_info(verbose, "    [TIME] Sim %d total: %.2fs (data=%.1f%%, 1rnd=%.1f%%, 2rnd=%.1f%%, comp=%.1f%%)\n",
           sim_id, sim_total_time,
           100 * (stage_times$data_gen %||% 0) / sim_total_time,
           100 * (stage_times$one_round %||% 0) / sim_total_time,
           100 * (stage_times$two_round %||% 0) / sim_total_time,
           100 * (stage_times$comparison %||% 0) / sim_total_time)

  # Combine results from list into a data.frame (single rbind at the end)
  results <- .bind_sim_result_list(results_list)

  # When TATE was fitted, derive the fallback count from returned
  # arm/fold diagnostics. Process-global counters lose child-worker updates
  # under fork parallelism and therefore cannot be used for topology-invariant
  # production QC. Keep the process-local fallback only for legacy runs that
  # did not request TATE.
  direct_degenerate_count <-
    .count_direct_tate_degenerate_fits(direct_tate_results)
  n_or_degenerate <- if (is.na(direct_degenerate_count)) {
    .get_fit_diagnostics()$or_degenerate_folds
  } else {
    direct_degenerate_count
  }
  if (nrow(results) > 0) {
    # Keep legacy potential-outcome rows for backwards compatibility, but mark
    # their estimand explicitly so diagnostics never rank them against TATE
    # rows.  All result constructors use bias = estimate - truth.
    results$estimand_scope <- ifelse(
      .is_tate_method(results$method),
      "tate",
      "treated_mean"
    )
    results$truth <- results$estimate - results$bias
    results$ci_lower <- results$estimate - Z_ALPHA_05 * results$se
    results$ci_upper <- results$estimate + Z_ALPHA_05 * results$se
    results$dgp_type <- dgp_type
    # DGP provenance: the misspecification strength actually applied (0 under
    # C1 and for the roce DGP), so cells from different DGP definitions cannot
    # be merged silently under the same config label.
    results$misspecification_strength <- data$misspecification_strength %||% 0
    results$outcome_family <- family
    results$heterogeneity_type <- result_heterogeneity_type
    results$or_degenerate_folds <- n_or_degenerate
    results$aggregation_lambda <- aggregation_lambda
    results$aggregation_cutoff <- 1 / aggregation_lambda
    # Persist the nuisance-CV grid size. Runtime experiments with different
    # grids must remain distinguishable in downstream aggregation and audits.
    results$nlambda_init <- as.integer(nlambda_init)
    results$nuisance_lambda_rule <- nuisance_lambda_rule
    results$n_bootstrap <- n_bootstrap
    results$n_weight_bootstrap <- n_weight_bootstrap
    results$M_tau <- M_tau
    results$M_tau_inference <- M_tau_inference
    results$hard_threshold_diagnostic_requested <-
      include_hard_threshold_diagnostic
    results$quadratic_bias_rule_requested <- include_quadratic_bias_rule
    results$deviation_mechanism <- deviation_mechanism
    weight_bootstrap_numeric <- c(
      "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
      "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
      "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
      "weight_bootstrap_failures"
    )
    for (column in setdiff(weight_bootstrap_numeric, names(results))) {
      results[[column]] <- NA_real_
    }
    if (!"weight_bootstrap_relearn_weights" %in% names(results)) {
      results$weight_bootstrap_relearn_weights <- NA
    }
    weight_bootstrap_character <- c(
      "weight_bootstrap_multiplier", "weight_bootstrap_screening_rule"
    )
    for (column in setdiff(weight_bootstrap_character, names(results))) {
      results[[column]] <- NA_character_
    }
    for (diagnostic_name in names(cell_diagnostics)) {
      results[[diagnostic_name]] <- cell_diagnostics[[diagnostic_name]]
    }
    results$rho_reuse_enabled <- !is.null(rho_reuse_reference)
    results$rho_reuse_changed_sources <- if (is.null(rho_reuse_reference)) {
      ""
    } else {
      paste(
        unique(as.character(rho_reuse_reference$changed_sources)),
        collapse = ";"
      )
    }
    results$data_generation_seconds <- stage_times$data_gen %||% NA_real_
    # These are wall-clock stages.  With concurrent treatment arms they cover
    # both arm-specific RoCE nuisance pipelines, so labeling them "mu1"
    # would be misleading.
    results$one_round_face_wall_seconds <-
      stage_times$one_round %||% NA_real_
    results$two_round_face_wall_seconds <-
      stage_times$two_round %||% NA_real_
    results$comparison_mu1_seconds <- stage_times$comparison %||% NA_real_
    results$tate_contrast_stage_seconds <- stage_times$tate %||% NA_real_
    results$simulation_elapsed_seconds <- sim_total_time
  }
  if (n_or_degenerate > 0L) {
    .warn_with_context(
      "sim %d (config=%s, p=%d): %d outcome-model fold(s) had a single-class response; used the constant degenerate fallback.",
      sim_id, config, p, n_or_degenerate)
  }

  if (isTRUE(return_fitted_tate)) {
    if (length(direct_tate_results) == 0L) {
      stop(
        paste0(
          "return_fitted_tate = TRUE requires at least one requested ",
          "cross-fit method with a TATE result."
        ),
        call. = FALSE
      )
    }
    attr(results, "roce_simulation_artifacts") <- list(
      data_split = data_split,
      tate_truth = data$mu1_true - data$mu0_true,
      direct_tate_results = direct_tate_results
    )
  }

  return(results)
}

.run_face_positive_rho_update <- function(
    rho, simulation_args, reuse_reference, keep_artifact = FALSE) {
  current <- do.call(
    run_single_simulation,
    c(simulation_args, list(
      ate_deviation = rho,
      n_deviated_sites = length(reuse_reference$changed_sources),
      return_fitted_tate = keep_artifact,
      rho_reuse_reference = reuse_reference
    ))
  )
  current_artifacts <- attr(current, "roce_simulation_artifacts")
  attr(current, "roce_simulation_artifacts") <- NULL
  current$rho <- rho
  list(result = current, artifacts = current_artifacts)
}

.run_face_rho_group <- function(
    simulation_args,
    rho_values = c(0, 0.5, 1, 1.5, 2, 2.5),
    changed_sources = "s1",
    artifact_rhos = numeric(0),
    positive_rho_workers = 1L) {
  caller <- ".run_face_rho_group"
  if (!is.list(simulation_args) || is.null(names(simulation_args)) ||
      any(!nzchar(names(simulation_args))) || anyDuplicated(names(simulation_args))) {
    stop(caller, ": simulation_args must be a uniquely named list.",
         call. = FALSE)
  }
  forbidden <- intersect(
    names(simulation_args),
    c(
      "ate_deviation", "n_deviated_sites", "return_fitted_tate",
      "rho_reuse_reference"
    )
  )
  if (length(forbidden) > 0L) {
    stop(
      caller, ": simulation_args must not override grouped field(s): ",
      paste(forbidden, collapse = ", "), ".", call. = FALSE
    )
  }
  rho_values <- as.numeric(rho_values)
  if (length(rho_values) < 2L || any(!is.finite(rho_values)) ||
      any(rho_values < 0) || anyDuplicated(rho_values) ||
      !identical(rho_values[[1L]], 0)) {
    stop(
      caller,
      ": rho_values must be unique, non-negative, and start with zero.",
      call. = FALSE
    )
  }
  methods <- as.character(simulation_args$methods %||% character(0))
  if (!identical(simulation_args$dgp_type %||% "face", "face") ||
      !isTRUE(simulation_args$estimate_ate) ||
      !"one_round_crossfit" %in% methods ||
      "two_round_crossfit" %in% methods) {
    stop(
      caller,
      ": grouped rho reuse requires FACE TATE with one_round_crossfit only.",
      call. = FALSE
    )
  }
  artifact_rhos <- as.numeric(artifact_rhos)
  if (any(!artifact_rhos %in% rho_values) || anyDuplicated(artifact_rhos)) {
    stop(caller, ": artifact_rhos must be a unique subset of rho_values.",
         call. = FALSE)
  }
  if (length(positive_rho_workers) != 1L ||
      !is.numeric(positive_rho_workers) || is.na(positive_rho_workers) ||
      !is.finite(positive_rho_workers) || positive_rho_workers < 1 ||
      abs(positive_rho_workers - round(positive_rho_workers)) >
        sqrt(.Machine$double.eps)) {
    stop(caller, ": positive_rho_workers must be a positive integer.",
         call. = FALSE)
  }
  positive_rho_workers <- min(
    as.integer(positive_rho_workers), length(rho_values) - 1L
  )

  results <- vector("list", length(rho_values))
  # Format each rho separately: vector formatting pads integer values (0.0),
  # whereas later scalar lookups use 0. Mixing the two leaves NULL slots and
  # appends new entries when all six fitted artifacts are requested.
  rho_keys <- vapply(rho_values, format, character(1L),
                     scientific = FALSE, trim = TRUE)
  if (anyDuplicated(rho_keys)) {
    stop(caller, ": rho values must have distinct formatted keys.",
         call. = FALSE)
  }
  names(results) <- rho_keys
  artifacts <- stats::setNames(
    vector("list", length(artifact_rhos)),
    rho_keys[match(artifact_rhos, rho_values)]
  )

  reference <- do.call(
    run_single_simulation,
    c(simulation_args, list(
      ate_deviation = 0,
      n_deviated_sites = 0L,
      return_fitted_tate = TRUE
    ))
  )
  reference_artifacts <- attr(reference, "roce_simulation_artifacts")
  if (is.null(reference_artifacts) ||
      is.null(reference_artifacts$direct_tate_results$one_round_crossfit)) {
    stop(caller, ": rho-zero reference did not return reusable artifacts.",
         call. = FALSE)
  }
  reference$rho <- 0
  if (0 %in% artifact_rhos) {
    artifacts[[rho_keys[[1L]]]] <-
      reference_artifacts
  }
  attr(reference, "roce_simulation_artifacts") <- NULL
  results[[1L]] <- reference

  reuse_reference <- list(
    data_split = reference_artifacts$data_split,
    fitted_tate =
      reference_artifacts$direct_tate_results$one_round_crossfit,
    changed_sources = unique(as.character(changed_sources))
  )
  # The rho-zero reference may initialize OpenMP worker pools while fitting
  # multiple sources. Forking afterwards and then entering an OpenMP nuisance-
  # CV kernel can deadlock on Linux. Each positive-rho update therefore uses
  # one source worker per fit. Independent rho updates may run concurrently in
  # fresh PSOCK processes, which retain OpenMP parallelism within each CV path
  # without inheriting the reference process's OpenMP state.
  positive_rho_args <- simulation_args
  positive_rho_args$n_cores_internal <- 1L
  positive_indices <- seq_along(rho_values)[-1L]
  positive_tasks <- lapply(positive_indices, function(index) {
    list(
      rho = rho_values[[index]],
      keep_artifact = rho_values[[index]] %in% artifact_rhos
    )
  })
  names(positive_tasks) <- names(results)[positive_indices]
  if (positive_rho_workers == 1L) {
    positive_updates <- lapply(positive_tasks, function(task) {
      .run_face_positive_rho_update(
        rho = task$rho,
        simulation_args = positive_rho_args,
        reuse_reference = reuse_reference,
        keep_artifact = task$keep_artifact
      )
    })
    positive_backend <- "serial"
  } else {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop(caller, ": parallel is required for positive-rho PSOCK workers.",
           call. = FALSE)
    }
    worker_package_library <- dirname(normalizePath(
      find.package("RoCE"), mustWork = TRUE
    ))
    cluster <- parallel::makePSOCKcluster(
      positive_rho_workers, outfile = ""
    )
    on.exit({
      if (!is.null(cluster)) parallel::stopCluster(cluster)
    }, add = TRUE)
    worker_packages <- parallel::clusterCall(
      cluster,
      function(package_library) {
        if (isNamespaceLoaded("RoCE")) {
          unloadNamespace("RoCE")
        }
        .libPaths(unique(c(package_library, .libPaths())))
        suppressPackageStartupMessages(library(
          "RoCE", lib.loc = package_library, character.only = TRUE
        ))
        list(
          package_path = normalizePath(
            find.package("RoCE"), mustWork = TRUE
          ),
          has_update = exists(
            ".run_face_positive_rho_update",
            envir = asNamespace("RoCE"), inherits = FALSE
          )
        )
      },
      worker_package_library
    )
    expected_worker_package <- normalizePath(
      file.path(worker_package_library, "RoCE"), mustWork = TRUE
    )
    worker_package_paths <- vapply(
      worker_packages, `[[`, character(1L), "package_path"
    )
    worker_has_update <- vapply(
      worker_packages, `[[`, logical(1L), "has_update"
    )
    if (any(worker_package_paths != expected_worker_package) ||
        !all(worker_has_update)) {
      stop(
        caller,
        ": positive-rho workers did not load the current RoCE package.",
        call. = FALSE
      )
    }
    parallel::clusterExport(
      cluster,
      c("positive_rho_args", "reuse_reference"),
      envir = environment()
    )
    positive_updates <- parallel::parLapply(
      cluster, positive_tasks, function(task) {
        try(
          get(
            ".run_face_positive_rho_update",
            envir = asNamespace("RoCE"), inherits = FALSE
          )(
            rho = task$rho,
            simulation_args = positive_rho_args,
            reuse_reference = reuse_reference,
            keep_artifact = task$keep_artifact
          ),
          silent = TRUE
        )
      }
    )
    positive_updates <- .validate_parallel_results(
      positive_updates, positive_tasks,
      caller = paste0(caller, " positive-rho PSOCK")
    )
    parallel::stopCluster(cluster)
    cluster <- NULL
    positive_backend <- "psock"
  }
  for (offset in seq_along(positive_indices)) {
    index <- positive_indices[[offset]]
    rho <- rho_values[[index]]
    update <- positive_updates[[offset]]
    if (rho %in% artifact_rhos) {
      artifacts[[rho_keys[[index]]]] <-
        update$artifacts
    }
    results[[index]] <- update$result
  }

  list(
    rho_values = rho_values,
    results = results,
    artifacts = artifacts,
    reuse = list(
      reference_rho = 0,
      changed_sources = reuse_reference$changed_sources,
      independent_reference_fit = TRUE,
      positive_rho_workers = positive_rho_workers,
      positive_rho_backend = positive_backend,
      positive_rho_source_workers_per_fit = 1L
    )
  )
}

#' Summarize simulation results into aggregate statistics
#'
#' Computes mean bias, standard deviation of bias, RMSE, mean SE,
#' coverage probability, and CI width statistics, grouped by method
#' and configuration.
#'
#' @param results Data frame of per-simulation results (output of
#'   \code{run_single_simulation}). Must contain columns: method, config,
#'   heterogeneity_type, n_total, K, bias, se, coverage, ci_width.
#' @return Data frame of summary statistics with one row per
#'   method-config-heterogeneity-n_total-K combination.
#' @export
summarize_results <- function(results) {
  # Handle empty or invalid results
  if (is.null(results) || nrow(results) == 0) {
    .warn_with_context("No results to summarize (empty data frame).")
    return(data.frame())
  }

  # Check required columns exist
  required_cols <- c("method", "config", "heterogeneity_type", "n_total", "K", "bias", "se", "coverage", "ci_width")
  missing_cols <- setdiff(required_cols, names(results))
  if (length(missing_cols) > 0) {
    .warn_with_context("Missing columns in results: %s",
                       paste(missing_cols, collapse = ", "))
    return(data.frame())
  }

  # Determine grouping columns (include estimand_type and dgp_type if present)
  group_cols <- c("method", "config", "heterogeneity_type", "n_total", "K")
  if ("estimand_type" %in% names(results)) {
    group_cols <- c(group_cols, "estimand_type")
  }
  if ("dgp_type" %in% names(results)) {
    group_cols <- c(group_cols, "dgp_type")
  }
  # Keep runtime and inference configurations separate when this general
  # summarizer is called directly on task-level output.
  group_cols <- unique(c(
    group_cols,
    intersect(
      c(
        "nlambda_init", "nuisance_lambda_rule", "n_bootstrap",
        "n_weight_bootstrap",
        "M_tau", "M_tau_inference"
      ),
      names(results)
    )
  ))
  group_formula <- stats::as.formula(
    paste("~ ", paste(group_cols, collapse = " + "))
  )

  # Calculate summary statistics by method and configuration
  # For bias and se, we need mean, sd, and rmse
  bias_stats <- aggregate(
    stats::update(group_formula, bias ~ .),
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      rmse = sqrt(mean(x^2, na.rm = TRUE))  # RMSE for bias
    )
  )

  # For se (IF-based), we only need mean
  se_stats <- aggregate(
    stats::update(group_formula, se ~ .),
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  # For coverage (IF-based), we only need mean (proportion of TRUE values)
  coverage_stats <- aggregate(
    stats::update(group_formula, coverage ~ .),
    data = results,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  # For CI width (IF-based), we need mean and sd
  ci_width_stats <- aggregate(
    stats::update(group_formula, ci_width ~ .),
    data = results,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE)
    )
  )

  n_success_stats <- aggregate(
    stats::update(group_formula, estimate ~ .),
    data = results,
    FUN = function(x) sum(is.finite(x))
  )

  # Build base summary data frame
  summary_df <- data.frame(
    method = bias_stats$method,
    config = bias_stats$config,
    heterogeneity_type = bias_stats$heterogeneity_type,
    n_total = bias_stats$n_total,
    K = bias_stats$K,
    bias_mean = bias_stats$bias[, "mean"],
    bias_sd = bias_stats$bias[, "sd"],
    rmse = bias_stats$bias[, "rmse"],
    se_mean = se_stats$se,
    coverage = coverage_stats$coverage,
    ci_width_mean = ci_width_stats$ci_width[, "mean"],
    ci_width_sd = ci_width_stats$ci_width[, "sd"],
    n_success = n_success_stats$estimate
  )

  # Add estimand_type column if present in results
  if ("estimand_type" %in% names(bias_stats)) {
    summary_df$estimand_type <- bias_stats$estimand_type
  }
  # Add dgp_type column if present in results
  if ("dgp_type" %in% names(bias_stats)) {
    summary_df$dgp_type <- bias_stats$dgp_type
  }
  extra_group_cols <- setdiff(group_cols, names(summary_df))
  for (column in extra_group_cols) {
    summary_df[[column]] <- bias_stats[[column]]
  }

  # Surface the graceful-degradation rate if present: mean number of
  # single-class outcome-model folds per simulation (0 = never; a large value
  # flags a saturated DGP that warrants revisiting -- see fit_glmnet_cv).
  if ("or_degenerate_folds" %in% names(results)) {
    or_degen_stats <- aggregate(
      stats::update(group_formula, or_degenerate_folds ~ .),
      data = results,
      FUN = function(x) mean(x, na.rm = TRUE)
    )
    summary_df$or_degenerate_folds_mean <- or_degen_stats$or_degenerate_folds
  }

  return(summary_df)
}


# =============================================================================
# SIMULATION STUDY ORCHESTRATION
# =============================================================================

#' Run a full simulation study with checkpoint/restart support
#'
#' Iterates over a parameter grid (n, K, p, config), running
#' \code{n_sims} Monte Carlo replicates per setting. Supports
#' SLURM preemption-safe checkpointing and nested parallelism.
#'
#' @param n_sims Integer number of MC replicates per setting.
#' @param n_total_vec Integer vector of sample sizes.
#' @param K_vec Integer vector of source-site counts.
#' @param p_vec Integer vector of covariate dimensions.
#' @param configs Character vector of configurations (e.g., "C1").
#' @param n_cores Integer number of available cores.
#' @param checkpoint_config List returned by \code{init_checkpoint_config()},
#'   or NULL to disable checkpointing. Fields: \code{file},
#'   \code{setting_interval}, \code{sim_interval}, \code{preempt_signal},
#'   \code{saved_signal}.
#' @param nlambda_init Integer lambda grid size for initial outcome CV.
#' @param nuisance_lambda_rule Nuisance-model cross-validation selection rule
#'   passed to every simulation replicate: \code{"min"} or \code{"1se"}.
#' @param nested_parallel Logical. Enable sim + source-site parallelism.
#' @param estimand_type "superpopulation" or "sample".
#' @param site_allocation "model", "uniform", or "balanced".
#' @param transform_type "strong", "mild", or "none".
#' @param outcome_type "binary" or "continuous".
#' @param heterogeneity_type "none", "mild", "strong", or "partial" for the
#'   RoCE DGP. FACE-DGP result metadata is derived from its explicit source
#'   deviation and effect-modification arguments.
#' @param shift_strength Numeric covariate shift multiplier.
#' @param n_folds Integer cross-fitting fold count.
#' @param use_lambda_cache Logical. If TRUE, enable lambda caching
#'   within cross-fitting nuisance-model loops.
#' @param aggregation_lambda Positive truncated-Wald penalty multiplier passed
#'   to every simulation replicate.
#' @param verbose_every Integer. In sequential mode, run detailed simulation
#'   logging every \code{verbose_every} simulations (default 10).
#' @param parallel_strategy Character. Nested parallel allocation strategy:
#'   \code{"outer_priority"} (default), \code{"balanced"}, or
#'   \code{"outer_only"}.
#' @param estimate_ate Logical. Also estimate ATE via A_val=0?
#' @param dgp_type "roce" or "face".
#' @param ate_deviation Additive source treatment-shift deviation on the FACE
#'   outcome linear-predictor scale. The historical argument name is retained
#'   for compatibility; for binary outcomes this is a log-odds shift, not the
#'   induced marginal TATE difference.
#' @param n_deviated_sites Integer deviated sites (FACE DGP only).
#' @param deviation_mechanism \code{"treated_arm"} (default) or
#'   \code{"both_arms"}: how the deviated sources deviate under the FACE DGP
#'   (see \code{\link{generate_face_data}}).
#' @param n_target Optional integer target-site sample size for explicit
#'   per-site allocation (FACE DGP only; paired with \code{n_source_sizes}).
#'   When supplied, the \code{n_total} and \code{K} grid axes collapse to this
#'   single allocation while \code{p} and \code{config} are still swept.
#' @param n_source_sizes Optional integer vector of per-site source sample
#'   sizes (FACE DGP only).
#' @param n_weight_bootstrap Number of opt-in weight-relearning bootstrap draws
#'   per primary soft direct-TATE fit; 0 disables the diagnostic.
#' @return Data frame of simulation results.
#' @export
run_simulation_study <- function(n_sims = 500,
                                n_total_vec = c(500, 1000, 2000),
                                K_vec = c(2, 3, 5),
                                p_vec = c(4, 8),
                                configs = c("C1", "C2", "C3", "C4"),
                                n_cores = 1,
                                checkpoint_config = NULL,
                                nlambda_init = LAMBDA_GRID_SIZE_STANDARD,
                                nested_parallel = FALSE,
                                estimand_type = "superpopulation",
                                site_allocation = "model",
                                transform_type = "mild",
                                outcome_type = "binary",
                                heterogeneity_type = "none",
                                shift_strength = ROCE_SHIFT_STRENGTH_DEFAULT,
                                n_folds = N_FOLDS_DEFAULT,
                                use_lambda_cache = TRUE,
                                aggregation_lambda = AGG_WALD_LAMBDA,
                                verbose_every = 10L,
                                parallel_strategy = c("outer_priority", "balanced", "outer_only"),
                                estimate_ate = FALSE,
                                dgp_type = "face",
                                ate_deviation    = 0.0,
                                n_deviated_sites = 0L,
                                deviation_mechanism = c("treated_arm", "both_arms"),
                                # Explicit per-site sample sizes (FACE DGP only)
                                n_target         = NULL,
                                n_source_sizes   = NULL,
                                nuisance_lambda_rule = c("min", "1se"),
                                n_weight_bootstrap = 0L) {
  deviation_mechanism <- match.arg(deviation_mechanism)

  parallel_strategy <- match.arg(parallel_strategy)
  nuisance_lambda_rule <- .match_nuisance_lambda_rule(
    nuisance_lambda_rule, "run_simulation_study",
    arg = "nuisance_lambda_rule"
  )
  n_weight_bootstrap <- .validate_weight_bootstrap_replicates(
    n_weight_bootstrap, "run_simulation_study"
  )
  if (length(aggregation_lambda) != 1L ||
      !is.finite(aggregation_lambda) || aggregation_lambda <= 0) {
    stop("aggregation_lambda must be one finite positive number.",
         call. = FALSE)
  }

  # Explicit per-site sample sizes (FACE DGP only): collapse the n_total and K
  # grid axes to the single allocation they imply, so the study still sweeps p
  # and config while every cell uses the requested per-site sizes.
  if (!is.null(n_source_sizes) || !is.null(n_target)) {
    if (dgp_type != "face") {
      stop("n_target / n_source_sizes are supported only for dgp_type = 'face'.",
           call. = FALSE)
    }
    site_sizes  <- resolve_face_site_sizes(NULL, n_target, n_source_sizes, K = NULL)
    n_total_vec <- site_sizes$n_total
    K_vec       <- site_sizes$K
    cat(sprintf(
      "  [per-site sizes] n_total -> %d, K -> %d (sources %s, target %d)\n",
      site_sizes$n_total, site_sizes$K,
      paste(site_sizes$n_source_sizes, collapse = ","), site_sizes$n_target))
  }

  validate_simulation_params(
    estimand_type      = estimand_type,
    site_allocation    = site_allocation,
    transform_type     = transform_type,
    outcome_type       = outcome_type,
    heterogeneity_type = heterogeneity_type,
    shift_strength     = shift_strength,
    n_folds            = n_folds,
    n_sims             = n_sims,
    n_total            = n_total_vec,
    K                  = K_vec,
    p                  = p_vec,
    config             = configs,
    dgp_type           = dgp_type,
    ate_deviation      = ate_deviation,
    n_deviated_sites   = n_deviated_sites,
    warn_ignored       = TRUE
  )

  # Resolve checkpoint config (NULL disables checkpointing)
  ckpt_file         <- checkpoint_config$file
  ckpt_interval     <- checkpoint_config$setting_interval %||% 1L
  sim_ckpt_interval <- checkpoint_config$sim_interval %||% n_sims
  preempt_file      <- checkpoint_config$preempt_signal
  saved_file        <- checkpoint_config$saved_signal
  use_checkpoint    <- !is.null(ckpt_file)

  # Create parameter grid
  param_grid <- expand.grid(
    n_total = n_total_vec,
    K = K_vec,
    p = p_vec,
    config = configs,
    stringsAsFactors = FALSE
  )
  total_settings <- nrow(param_grid)

  cat(sprintf("\n[RoCE] Simulation study: %d settings, %d sims each, %d cores\n",
              total_settings, n_sims, n_cores))
  cat(sprintf("  n_total: %s | K: %s | p: %s | configs: %s\n",
              paste(n_total_vec, collapse = ","),
              paste(K_vec, collapse = ","),
              paste(p_vec, collapse = ","),
              paste(configs, collapse = ",")))
  cat(sprintf("  estimand: %s | outcome: %s | heterogeneity: %s\n",
              estimand_type, outcome_type, heterogeneity_type))
  cat(sprintf("  checkpoint: %s (setting every %d, sim every %d)\n",
              if (use_checkpoint) "ON" else "OFF", ckpt_interval, sim_ckpt_interval))

  # Try to load checkpoint
  checkpoint <- if (use_checkpoint) load_checkpoint(ckpt_file) else NULL

  # Variables for simulation-level resume
  resume_sim_results <- NULL
  resume_sim_idx <- 0

  if (!is.null(checkpoint)) {
    all_results <- checkpoint$all_results

    if (!is.null(checkpoint$sim_results) && !is.null(checkpoint$current_sim_idx)) {
      start_idx <- checkpoint$current_setting_idx
      resume_sim_results <- checkpoint$sim_results
      resume_sim_idx <- checkpoint$current_sim_idx
      cat(sprintf("  Resuming setting %d from simulation %d\n",
                  start_idx, resume_sim_idx + 1))
    } else {
      start_idx <- checkpoint$current_setting_idx + 1
      cat(sprintf("  Resuming from setting %d/%d\n", start_idx, total_settings))
    }

    if (checkpoint$current_setting_idx >= total_settings &&
        is.null(checkpoint$sim_results)) {
      cat("  Checkpoint indicates all settings complete. Finalizing.\n")
      final_results <- .bind_sim_result_list(all_results)
      if (use_checkpoint) cleanup_checkpoint(ckpt_file, preempt_file, saved_file)
      return(final_results)
    }
  } else {
    all_results <- list()
    start_idx <- 1
    cat("  Starting fresh (no checkpoint found)\n")
  }

  # Set up parallel processing
  cl <- NULL
  n_cores_internal_parallel <- 1
  n_cores_outer <- n_cores

  if (n_cores > 1) {
    if (nested_parallel) {
      max_k <- max(K_vec)

      if (parallel_strategy == "outer_only" || max_k <= 1) {
        outer_cores <- n_cores
        inner_cores <- 1L
      } else if (parallel_strategy == "outer_priority") {
        # High-core default: prioritize simulation-level throughput.
        # Keep inner parallelism minimal (<=2) to avoid nested overhead.
        inner_cores <- if (n_cores >= 16 && max_k >= 2) min(2L, max_k) else 1L
        outer_cores <- max(1L, floor(n_cores / inner_cores))
      } else {
        # balanced
        inner_cores <- min(max_k, n_cores)
        outer_cores <- max(2L, floor(n_cores / inner_cores))
        inner_cores <- max(1L, floor(n_cores / outer_cores))
        inner_cores <- min(inner_cores, max_k)
      }

      if (outer_cores < 2L) {
        outer_cores <- n_cores
        inner_cores <- 1L
      }

      n_cores_outer <- outer_cores
      n_cores_internal_parallel <- inner_cores

      cat(sprintf("  Nested parallel (%s): %d outer x %d inner (of %d total)\n",
                  parallel_strategy, outer_cores, inner_cores, n_cores))
    } else {
      n_cores_outer <- n_cores
      cat(sprintf("  Using %d cores for simulation-level parallelization\n", n_cores))
    }

    cl <- parallel::makeCluster(n_cores_outer)
    on.exit({
      if (!is.null(cl)) parallel::stopCluster(cl)
    }, add = TRUE)
    doParallel::registerDoParallel(cl)

    parallel::clusterEvalQ(cl, {
      if (requireNamespace("RoCE", quietly = TRUE)) {
        library(RoCE)
      } else if (requireNamespace("devtools", quietly = TRUE)) {
        devtools::load_all(".")
      } else {
        stop("Workers require either installed RoCE package or devtools.")
      }
    })

    parallel::clusterExport(cl, c(
      "run_single_simulation",
      ".bind_sim_result_list", ".abort_with_context",
      ".warn_with_context", ".require_method_result",
      "n_cores_internal_parallel", "nlambda_init",
      "nuisance_lambda_rule", "estimand_type",
      "site_allocation", "transform_type", "outcome_type",
      "heterogeneity_type", "shift_strength", "n_folds",
      "aggregation_lambda",
      "n_weight_bootstrap",
      "estimate_ate", "dgp_type", "ate_deviation", "n_deviated_sites"
    ), envir = environment())
  }

  # Main loop over parameter settings
  for (i in start_idx:total_settings) {

    # Check for preemption signal
    if (use_checkpoint && check_preempt_signal(preempt_file)) {
      cat(sprintf("\n[WARN] Preemption signal at setting %d/%d. Saving checkpoint...\n",
                  i, total_settings))
      state <- list(
        all_results = all_results,
        current_setting_idx = i - 1,
        total_settings = total_settings,
        timestamp = Sys.time()
      )
      checkpoint_saved <- save_checkpoint(state, ckpt_file)
      if (checkpoint_saved) signal_checkpoint_saved(saved_file)
      if (!is.null(cl)) {
        parallel::stopCluster(cl)
        cl <- NULL
      }
      cat("  Checkpoint saved. Exiting for requeue.\n")
      quit(save = "no", status = 0)
    }

    params <- param_grid[i, ]
    cat(sprintf("\n[%s] Setting %d/%d: n=%d, K=%d, p=%d, config=%s\n",
                format(Sys.time(), "%H:%M:%S"), i, total_settings,
                params$n_total, params$K, params$p, params$config))

    methods <- c("two_round_crossfit", "one_round_crossfit",
                 "target_only", "sample_size", "inverse_variance",
                 "federated_dr", "pooled_dr", "tilted_aipw",
                 "oracle_dr")

    # Resume data if available
    if (i == start_idx && !is.null(resume_sim_results)) {
      sim_results <- resume_sim_results
      sim_start_idx <- resume_sim_idx + 1
      cat(sprintf("  Resuming with %d completed simulations\n", resume_sim_idx))
      resume_sim_results <- NULL
      resume_sim_idx <- 0
    } else {
      sim_results <- list()
      sim_start_idx <- 1
    }

    setting_start_time <- Sys.time()
    completed_sims_timing <- numeric(0)

    if (n_cores_outer > 1) {
        # --- Parallel execution (batched) ---
        batch_size <- n_cores_outer * 2
        first_batch <- ceiling(sim_start_idx / batch_size)
        n_batches <- ceiling(n_sims / batch_size)

        for (batch in first_batch:n_batches) {
          batch_start_time <- Sys.time()
          start_sim <- max((batch - 1) * batch_size + 1, sim_start_idx)
          end_sim <- min(batch * batch_size, n_sims)
          batch_ids <- start_sim:end_sim

          elapsed_secs <- as.numeric(difftime(Sys.time(), setting_start_time, units = "secs"))
          sims_done <- start_sim - sim_start_idx
          eta_str <- if (sims_done > 0) {
            sprintf("ETA: %.1f min",
                    elapsed_secs / sims_done * (n_sims - start_sim + 1) / 60)
          } else {
            "ETA: calculating..."
          }

          pct <- end_sim / n_sims * 100
          bar_width <- 20
          filled <- floor(bar_width * pct / 100)
          bar <- paste0("[", strrep("\u2588", filled),
                        strrep("\u2591", bar_width - filled), "]")
          cat(sprintf("  [%s] %s %.0f%% | Sims %d-%d/%d | Batch %d/%d | %s\n",
                      format(Sys.time(), "%H:%M:%S"), bar, pct,
                      start_sim, end_sim, n_sims, batch, n_batches, eta_str))

          batch_results <- parallel::parLapply(cl, batch_ids, function(sim_id) {
            run_single_simulation(
              sim_id, params$n_total, params$K, params$p,
              params$config, methods,
              verbose = FALSE,
              n_cores_internal = n_cores_internal_parallel,
              nlambda_init = nlambda_init,
              nuisance_lambda_rule = nuisance_lambda_rule,
              estimand_type = estimand_type,
              site_allocation = site_allocation,
              transform_type = transform_type,
              outcome_type = outcome_type,
              heterogeneity_type = heterogeneity_type,
              shift_strength = shift_strength,
              n_folds = n_folds,
              use_lambda_cache = use_lambda_cache,
              aggregation_lambda = aggregation_lambda,
              estimate_ate = estimate_ate,
              dgp_type = dgp_type,
              ate_deviation = ate_deviation,
              n_deviated_sites = n_deviated_sites,
              deviation_mechanism = deviation_mechanism,
              n_target = n_target,
              n_source_sizes = n_source_sizes,
              n_weight_bootstrap = n_weight_bootstrap)
          })

          batch_elapsed <- as.numeric(
            difftime(Sys.time(), batch_start_time, units = "secs"))
          n_success <- sum(vapply(batch_results,
                                  function(r) is.data.frame(r) && nrow(r) > 0,
                                  logical(1)))
          cat(sprintf("    Batch done: %.1fs (%.2fs/sim) | %d/%d succeeded\n",
                      batch_elapsed, batch_elapsed / length(batch_ids),
                      n_success, length(batch_ids)))

          for (j in seq_along(batch_results)) {
            sim_results[[batch_ids[j]]] <- batch_results[[j]]
          }

          # Sim-level checkpoint
          if (use_checkpoint &&
              end_sim %% sim_ckpt_interval == 0 && end_sim < n_sims) {
            state <- list(
              all_results = all_results, current_setting_idx = i,
              total_settings = total_settings, sim_results = sim_results,
              current_sim_idx = end_sim, n_sims = n_sims,
              timestamp = Sys.time()
            )
            save_checkpoint(state, ckpt_file, sim_level = TRUE)
          }
        }
    } else {
        # --- Sequential execution ---
        n_cores_internal <- min(params$K,
                               max(1, parallel::detectCores() - 1))
        if (n_cores_internal > 1) {
          cat(sprintf("  Using %d internal cores for source-site parallelization\n",
                      n_cores_internal))
        }

        for (sim_id in sim_start_idx:n_sims) {
          sim_start <- Sys.time()

          sims_done <- sim_id - sim_start_idx
          eta_str <- if (sims_done > 0 && length(completed_sims_timing) > 0) {
            sprintf("ETA: %.1f min",
                    mean(completed_sims_timing) * (n_sims - sim_id + 1) / 60)
          } else {
            "ETA: calculating..."
          }

          if (sim_id %% 10 == 1 || sim_id == sim_start_idx) {
            cat(sprintf("  [%s] Progress: sim %d/%d (%.0f%%) | %s\n",
                        format(Sys.time(), "%H:%M:%S"),
                        sim_id, n_sims, sim_id / n_sims * 100, eta_str))
          }

          sim_results[[sim_id]] <- run_single_simulation(
            sim_id, params$n_total, params$K, params$p,
            params$config, methods,
            verbose = (sim_id == sim_start_idx) || ((sim_id - sim_start_idx) %% max(1L, as.integer(verbose_every)) == 0L),
            n_cores_internal = n_cores_internal,
            nlambda_init = nlambda_init,
            nuisance_lambda_rule = nuisance_lambda_rule,
            estimand_type = estimand_type,
            site_allocation = site_allocation,
            transform_type = transform_type,
            outcome_type = outcome_type,
            heterogeneity_type = heterogeneity_type,
            shift_strength = shift_strength,
            n_folds = n_folds,
            use_lambda_cache = use_lambda_cache,
            aggregation_lambda = aggregation_lambda,
            estimate_ate = estimate_ate,
            dgp_type = dgp_type,
            ate_deviation = ate_deviation,
            n_deviated_sites = n_deviated_sites,
            deviation_mechanism = deviation_mechanism,
            n_target = n_target,
            n_source_sizes = n_source_sizes,
            n_weight_bootstrap = n_weight_bootstrap
          )

          sim_elapsed <- as.numeric(difftime(Sys.time(), sim_start, units = "secs"))
          completed_sims_timing <- c(completed_sims_timing, sim_elapsed)
          if (length(completed_sims_timing) > 10) {
            completed_sims_timing <- utils::tail(completed_sims_timing, 10)
          }

          # Sim-level checkpoint
          if (use_checkpoint &&
              sim_id %% sim_ckpt_interval == 0 && sim_id < n_sims) {
            state <- list(
              all_results = all_results, current_setting_idx = i,
              total_settings = total_settings, sim_results = sim_results,
              current_sim_idx = sim_id, n_sims = n_sims,
              timestamp = Sys.time()
            )
            save_checkpoint(state, ckpt_file, sim_level = TRUE)
          }
        }
    }

    setting_results <- .bind_sim_result_list(sim_results)
    all_results[[i]] <- setting_results

    setting_elapsed <- as.numeric(
      difftime(Sys.time(), setting_start_time, units = "secs"))
    cat(sprintf("  Setting %d/%d done: %d rows, %.1f min (%.2f s/sim)\n",
                i, total_settings, nrow(setting_results),
                setting_elapsed / 60, setting_elapsed / n_sims))

    # Setting-level checkpoint
    if (use_checkpoint && (i %% ckpt_interval == 0 || i == total_settings)) {
      state <- list(
        all_results = all_results, current_setting_idx = i,
        total_settings = total_settings, timestamp = Sys.time()
      )
      save_checkpoint(state, ckpt_file)
    }
  }

  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }

  final_results <- .bind_sim_result_list(all_results)

  if (is.null(final_results) || nrow(final_results) == 0) {
    cat("\n[WARN] No simulation results collected!\n")
    final_results <- data.frame()
  }

  if (use_checkpoint) cleanup_checkpoint(ckpt_file, preempt_file, saved_file)

  total_rows <- if (is.data.frame(final_results)) nrow(final_results) else 0
  cat(sprintf("\nSimulation study completed: %d rows, %d settings\n",
              total_rows, total_settings))

  return(final_results)
}
