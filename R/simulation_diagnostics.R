# Simulation result diagnostics ------------------------------------------------
#
# These helpers deliberately operate on replicate-level output.  Summary-only
# checks cannot distinguish a genuinely incomplete setting from non-finite
# values that were silently removed by na.rm = TRUE.

.simulation_diagnostic_group_columns <- function(results) {
  preferred <- c(
    "experiment", "dgp_type", "outcome_family", "config", "heterogeneity_type",
    "estimand_type", "p", "K", "rho", "cutoff",
    "n_site", "n_folds", "nlambda_init", "nuisance_lambda_rule",
    "n_bootstrap", "n_weight_bootstrap", "M_tau", "M_tau_inference",
    "estimand_scope", "method"
  )
  intersect(preferred, names(results))
}

.tate_benchmark_methods <- function() {
  c(
    "target_only_ate", "sample_size_ate", "inverse_variance_ate",
    "federated_dr_ate", "pooled_dr_ate"
  )
}

# Frozen production method-row set of one replicate (HISTORY #0009): the
# treated-mean rows, the arm-wise diagnostic, the primary common-weight TATE,
# the gated sensitivity/diagnostic rules, and the TATE benchmarks. Audits and
# tests compare a task's rows against this set.
.tate_production_method_rows <- function(include_hard_threshold = FALSE,
                                         include_quadratic_bias = TRUE) {
  c(
    "one_round_crossfit", "target_only", "sample_size", "inverse_variance",
    "federated_dr", "pooled_dr",
    "one_round_crossfit_ate_armwise", "one_round_crossfit_ate",
    if (isTRUE(include_quadratic_bias)) "one_round_crossfit_ate_quadratic_bias",
    if (isTRUE(include_hard_threshold)) "one_round_crossfit_ate_hard_threshold",
    .tate_benchmark_methods()
  )
}

.face_deviation_mechanisms <- function() c("treated_arm", "both_arms")

# Columns that must carry numeric storage in a simulation result frame.
# .read_simulation_result_files() restores this storage after a CSV round trip
# and the production validator enforces it.
.simulation_numeric_columns <- function() {
  unique(c(
    "estimate", "bias", "se", "truth", "ci_lower", "ci_upper",
    "ci_width", "min_site_arm_outcome_cell_n",
    "min_target_arm_outcome_cell_n",
    "n_site_arm_outcome_cells_below_8", "dr_weight_n",
    "dr_weight_n_clipped", "dr_weight_fraction_clipped",
    "dr_weight_max_site_fraction_clipped",
    "dr_weight_min_before_clipping", "dr_weight_max_before_clipping",
    "comparison_se_analytic", "comparison_se_bootstrap",
    "comparison_n_bootstrap",
    "n_weight_bootstrap", "variance_weight_relearn_bootstrap",
    "se_weight_relearn_bootstrap", "se_fixed_weight_bootstrap",
    "weight_uncertainty_ratio", "weight_relearn_n_bootstrap",
    "weight_bootstrap_seed", "weight_bootstrap_failures",
    .direct_tate_aggregation_diagnostic_columns()
  ))
}

.direct_tate_aggregation_diagnostic_columns <- function() {
  c(
    "target_anchor_weight", "mean_abs_source_weight",
    "max_abs_source_weight", "max_wald_statistic",
    "mean_wald_statistic", "penalized_source_fold_fraction",
    "max_weight_optimizer_iterations", "max_weight_psd_ridge",
    "weight_psd_ridge_fold_fraction",
    "se_fixed_weights", "weight_layer_indirect_variance",
    "weight_layer_cross_term", "weight_layer_kink_cells",
    "inference_logit_truncated", "inference_logit_truncation_fraction",
    "inference_max_abs_logit", "inference_safety_clip_count"
  )
}

.validate_face_production_scientific_metadata <- function(results,
                                                          expected_n_folds = 5L) {
  if (length(expected_n_folds) != 1L || !is.numeric(expected_n_folds) ||
      !is.finite(expected_n_folds) || expected_n_folds < 2 ||
      expected_n_folds != as.integer(expected_n_folds)) {
    stop("expected_n_folds must be one integer of at least 2.", call. = FALSE)
  }
  required <- c(
    "dgp_type", "outcome_family", "heterogeneity_type", "estimand_type",
    "config", "rho", "deviation_mechanism", "p", "K", "n_site", "n_folds",
    "estimand_scope", "method"
  )
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop(
      "FACE production results are missing scientific metadata: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyNA(results[required])) {
    stop("FACE production scientific metadata must be nonmissing.",
         call. = FALSE)
  }
  numeric_metadata <- c("rho", "p", "K", "n_site", "n_folds")
  if (!all(vapply(results[numeric_metadata], is.numeric, logical(1L)))) {
    stop("FACE production numeric design metadata have invalid storage.",
         call. = FALSE)
  }
  mechanism <- as.character(results$deviation_mechanism)
  expected_heterogeneity <- vapply(
    seq_len(nrow(results)),
    function(i) {
      rho <- results$rho[[i]]
      if (!is.finite(rho) || rho < 0 ||
          !mechanism[[i]] %in% .face_deviation_mechanisms()) {
        return(NA_character_)
      }
      .face_heterogeneity_type(
        ate_deviation = rho,
        n_deviated_sites = as.integer(rho > 0),
        deviation_mechanism = mechanism[[i]]
      )
    },
    character(1L)
  )
  expected_scope <- ifelse(
    .is_tate_method(results$method), "tate", "treated_mean"
  )
  valid <- is.finite(results$rho) &
    results$rho %in% c(0, 0.5, 1, 1.5, 2, 2.5) &
    is.finite(results$p) & is.finite(results$K) &
    is.finite(results$n_site) & is.finite(results$n_folds) &
    results$p == 100L & results$dgp_type == "face" &
    results$config %in% c("C1", "C2", "C3") &
    results$K %in% c(2L, 4L, 8L) & results$n_site == 1000L &
    results$n_folds == as.integer(expected_n_folds) &
    results$outcome_family == "binomial" &
    results$estimand_type == "superpopulation" &
    results$estimand_scope == expected_scope &
    !is.na(expected_heterogeneity) &
    results$heterogeneity_type == expected_heterogeneity
  if (!all(valid)) {
    stop(
      paste0(
        "FACE production scientific metadata must use the audited p=100 ",
        "design grid, the binary FACE DGP, the superpopulation estimand, ",
        "method-consistent estimand scopes, a known deviation mechanism, ",
        "and rho-consistent heterogeneity labels."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.simulation_coverage_as_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  if (is.numeric(x)) {
    result <- rep(NA, length(x))
    result[!is.na(x) & x == 1] <- TRUE
    result[!is.na(x) & x == 0] <- FALSE
    return(result)
  }
  normalized <- toupper(trimws(as.character(x)))
  result <- rep(NA, length(normalized))
  result[normalized %in% c("TRUE", "T", "1")] <- TRUE
  result[normalized %in% c("FALSE", "F", "0")] <- FALSE
  result
}

.collapse_diagnostic_flags <- function(flags) {
  active <- names(flags)[vapply(flags, isTRUE, logical(1L))]
  if (length(active) == 0L) "ok" else paste(active, collapse = ";")
}

.simulation_group_key <- function(data, columns) {
  grouping_data <- data[columns]
  grouping_data[] <- lapply(grouping_data, function(values) {
    values <- as.character(values)
    missing <- is.na(values) | !nzchar(trimws(values))
    values[missing] <- "<unspecified>"
    values
  })
  interaction(grouping_data, drop = TRUE, lex.order = TRUE)
}

.bind_simulation_result_frames <- function(parts) {
  if (!is.list(parts) || length(parts) == 0L ||
      any(!vapply(parts, is.data.frame, logical(1L)))) {
    stop("parts must be a non-empty list of data frames.", call. = FALSE)
  }
  all_columns <- unique(unlist(lapply(parts, names), use.names = FALSE))
  aligned <- lapply(parts, function(part) {
    missing_columns <- setdiff(all_columns, names(part))
    for (column in missing_columns) {
      part[[column]] <- rep(NA, nrow(part))
    }
    part[all_columns]
  })
  result <- do.call(rbind, aligned)
  rownames(result) <- NULL
  result
}

.read_simulation_result_files <- function(files) {
  if (!is.character(files) || length(files) == 0L || any(!file.exists(files))) {
    stop("files must name one or more existing simulation CSVs.", call. = FALSE)
  }
  results <- .bind_simulation_result_frames(lapply(
    files, utils::read.csv, stringsAsFactors = FALSE
  ))
  # A numeric column that is NA in every row (a disabled diagnostic, such as
  # the weight-relearn bootstrap when n_weight_bootstrap = 0) is written as
  # "NA" and read back as logical.  Restore the declared storage; a column with
  # any non-numeric value is left alone so the production validator still
  # rejects it.
  for (column in intersect(.simulation_numeric_columns(), names(results))) {
    values <- results[[column]]
    if (is.logical(values) && all(is.na(values))) {
      results[[column]] <- as.numeric(values)
    }
  }
  results
}

.summarize_group_nuisance_diagnostics <- function(group) {
  primary_count_columns <- c(
    "face_initial_dr_nonconverged",
    "face_calibrated_dr_nonconverged",
    "face_calibrated_outcome_nonconverged"
  )
  count_columns <- c(
    primary_count_columns,
    "face_initial_dr_line_search_failures",
    "face_calibrated_dr_line_search_failures",
    "face_calibrated_outcome_line_search_failures",
    "face_initial_dr_cv_invalid_fold_fits",
    "face_initial_dr_cv_invalid_lambdas",
    "face_initial_dr_cv_path_tail_skipped_fold_fits",
    "face_calibrated_dr_cv_invalid_fold_fits",
    "face_calibrated_dr_cv_invalid_lambdas",
    "face_calibrated_dr_cv_path_tail_skipped_fold_fits",
    "face_calibrated_outcome_cv_invalid_fold_fits",
    "face_calibrated_outcome_cv_invalid_lambdas",
    "face_calibrated_outcome_cv_path_tail_skipped_fold_fits",
    "face_initial_dr_support_floor_applied",
    "face_calibrated_dr_support_floor_applied",
    "face_mu1_initial_dr_nonconverged",
    "face_mu0_initial_dr_nonconverged",
    "face_mu1_calibrated_dr_nonconverged",
    "face_mu0_calibrated_dr_nonconverged",
    "face_mu1_calibrated_outcome_nonconverged",
    "face_mu0_calibrated_outcome_nonconverged",
    "face_mu1_initial_dr_line_search_failures",
    "face_mu0_initial_dr_line_search_failures",
    "face_mu1_calibrated_dr_line_search_failures",
    "face_mu0_calibrated_dr_line_search_failures",
    "face_mu1_calibrated_outcome_line_search_failures",
    "face_mu0_calibrated_outcome_line_search_failures",
    "face_mu1_initial_dr_cv_invalid_fold_fits",
    "face_mu0_initial_dr_cv_invalid_fold_fits",
    "face_mu1_initial_dr_cv_invalid_lambdas",
    "face_mu0_initial_dr_cv_invalid_lambdas",
    "face_mu1_initial_dr_cv_path_tail_skipped_fold_fits",
    "face_mu0_initial_dr_cv_path_tail_skipped_fold_fits",
    "face_mu1_calibrated_dr_cv_invalid_fold_fits",
    "face_mu0_calibrated_dr_cv_invalid_fold_fits",
    "face_mu1_calibrated_dr_cv_invalid_lambdas",
    "face_mu0_calibrated_dr_cv_invalid_lambdas",
    "face_mu1_calibrated_dr_cv_path_tail_skipped_fold_fits",
    "face_mu0_calibrated_dr_cv_path_tail_skipped_fold_fits",
    "face_mu1_calibrated_outcome_cv_invalid_fold_fits",
    "face_mu0_calibrated_outcome_cv_invalid_fold_fits",
    "face_mu1_calibrated_outcome_cv_invalid_lambdas",
    "face_mu0_calibrated_outcome_cv_invalid_lambdas",
    "face_mu1_calibrated_outcome_cv_path_tail_skipped_fold_fits",
    "face_mu0_calibrated_outcome_cv_path_tail_skipped_fold_fits",
    "face_mu1_initial_dr_support_floor_applied",
    "face_mu0_initial_dr_support_floor_applied",
    "face_mu1_calibrated_dr_support_floor_applied",
    "face_mu0_calibrated_dr_support_floor_applied"
  )
  iteration_columns <- c(
    "face_max_initial_dr_iterations",
    "face_max_calibrated_dr_iterations",
    "face_max_calibrated_outcome_iterations"
  )
  continuous_max_columns <- c(
    "face_max_initial_dr_update_ratio",
    "face_max_initial_dr_abs_coefficient",
    "face_max_calibrated_dr_update_ratio",
    "face_max_calibrated_dr_abs_coefficient",
    "face_max_calibrated_outcome_update_ratio",
    "face_max_calibrated_outcome_abs_coefficient",
    "face_max_initial_dr_support_floor",
    "face_max_calibrated_dr_support_floor"
  )
  available_counts <- intersect(count_columns, names(group))
  available_iterations <- intersect(iteration_columns, names(group))
  available_continuous <- intersect(continuous_max_columns, names(group))
  method_values <- if ("method" %in% names(group)) {
    unique(as.character(group$method))
  } else {
    character(0)
  }
  direct_tate_diagnostics_required <-
    length(method_values) == 1L &&
    method_values %in% c(
      "one_round_crossfit_ate", "two_round_crossfit_ate",
      "one_round_crossfit_ate_hard_threshold",
      "two_round_crossfit_ate_hard_threshold",
      "one_round_crossfit_ate_quadratic_bias",
      "two_round_crossfit_ate_quadratic_bias"
    )

  count_matrix <- if (length(available_counts) > 0L) {
    as.matrix(data.frame(lapply(
      group[available_counts], function(values) suppressWarnings(as.numeric(values))
    )))
  } else {
    matrix(numeric(0), nrow = nrow(group), ncol = 0L)
  }
  iteration_matrix <- if (length(available_iterations) > 0L) {
    as.matrix(data.frame(lapply(
      group[available_iterations],
      function(values) suppressWarnings(as.numeric(values))
    )))
  } else {
    matrix(numeric(0), nrow = nrow(group), ncol = 0L)
  }
  continuous_matrix <- if (length(available_continuous) > 0L) {
    as.matrix(data.frame(lapply(
      group[available_continuous],
      function(values) suppressWarnings(as.numeric(values))
    )))
  } else {
    matrix(numeric(0), nrow = nrow(group), ncol = 0L)
  }

  observed_count_rows <- if (ncol(count_matrix) > 0L) {
    rowSums(!is.na(count_matrix)) > 0L
  } else {
    rep(FALSE, nrow(group))
  }
  invalid_counts <- if (any(observed_count_rows)) {
    observed <- count_matrix[!is.na(count_matrix)]
    any(!is.finite(observed) | observed < 0 |
          observed != floor(observed))
  } else {
    FALSE
  }
  complete_primary_counts <-
    all(primary_count_columns %in% available_counts) &&
    all(vapply(group[primary_count_columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      all(is.finite(values) & values >= 0 & values == floor(values))
    }, logical(1L)))
  invalid_counts <- invalid_counts ||
    (direct_tate_diagnostics_required && !complete_primary_counts)
  invalid_iterations <- if (ncol(iteration_matrix) > 0L &&
      any(!is.na(iteration_matrix))) {
    observed <- iteration_matrix[!is.na(iteration_matrix)]
    any(!is.finite(observed) | observed < 0 |
          observed != floor(observed))
  } else {
    FALSE
  }
  invalid_continuous <- if (ncol(continuous_matrix) > 0L &&
      any(!is.na(continuous_matrix))) {
    observed <- continuous_matrix[!is.na(continuous_matrix)]
    any(!is.finite(observed) | observed < 0)
  } else {
    FALSE
  }

  count_totals <- setNames(rep(NA_real_, length(count_columns)), count_columns)
  for (column in available_counts) {
    values <- count_matrix[, match(column, available_counts)]
    if (any(!is.na(values))) {
      count_totals[[column]] <- sum(values[is.finite(values)])
    }
  }
  iteration_maxima <- setNames(
    rep(NA_real_, length(iteration_columns)), iteration_columns
  )
  for (column in available_iterations) {
    values <- iteration_matrix[, match(column, available_iterations)]
    if (any(is.finite(values))) {
      iteration_maxima[[column]] <- max(values[is.finite(values)])
    }
  }
  continuous_maxima <- setNames(
    rep(NA_real_, length(continuous_max_columns)), continuous_max_columns
  )
  for (column in available_continuous) {
    values <- continuous_matrix[, match(column, available_continuous)]
    if (any(is.finite(values))) {
      continuous_maxima[[column]] <- max(values[is.finite(values)])
    }
  }

  available_primary_counts <- intersect(
    primary_count_columns, available_counts
  )
  row_nonconvergence <- if (length(available_primary_counts) > 0L) {
    primary_matrix <- count_matrix[
      , match(available_primary_counts, available_counts), drop = FALSE
    ]
    rowSums(replace(primary_matrix, !is.finite(primary_matrix), 0)) > 0
  } else {
    rep(FALSE, nrow(group))
  }
  primary_count_totals <- count_totals[primary_count_columns]
  total_nonconvergence <- sum(
    primary_count_totals[is.finite(primary_count_totals)]
  )

  c(
    nuisance_diagnostic_replications = sum(observed_count_rows),
    nuisance_nonconverged_replications = sum(row_nonconvergence),
    nuisance_nonconverged_fits_total = total_nonconvergence,
    setNames(
      count_totals,
      paste0(names(count_totals), "_total")
    ),
    setNames(
      iteration_maxima,
      paste0(names(iteration_maxima), "_observed")
    ),
    setNames(
      continuous_maxima,
      paste0(names(continuous_maxima), "_observed")
    ),
    invalid_nuisance_diagnostics =
      invalid_counts || invalid_iterations || invalid_continuous
  )
}

.summarize_weight_relearn_bootstrap_diagnostics <- function(
    group, method_values, identity_tolerance) {
  numeric_fields <- c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
    "weight_bootstrap_failures"
  )
  character_fields <- c(
    "weight_bootstrap_multiplier", "weight_bootstrap_screening_rule"
  )
  logical_field <- "weight_bootstrap_relearn_weights"
  configured <- if ("n_weight_bootstrap" %in% names(group)) {
    suppressWarnings(as.numeric(group$n_weight_bootstrap))
  } else {
    rep(0, nrow(group))
  }
  valid_config <- all(is.finite(configured)) &&
    all(configured >= 0 & configured == floor(configured)) &&
    all(configured != 1) && length(unique(configured)) == 1L
  B <- if (valid_config) unique(configured) else NA_real_
  bootstrap_primary <- length(method_values) == 1L &&
    identical(method_values, "one_round_crossfit_ate")

  schema_fields <- c(numeric_fields, character_fields, logical_field)
  present_fields <- intersect(schema_fields, names(group))
  any_schema <- length(present_fields) > 0L
  has_schema <- length(present_fields) == length(schema_fields)
  observed <- rep(FALSE, nrow(group))
  if (any_schema) {
    observed <- Reduce(`|`, lapply(
      group[present_fields], function(value) !is.na(value)
    ))
  }
  valid_rows <- rep(FALSE, nrow(group))
  variance_identity_error <- NA_real_
  ratio_identity_error <- NA_real_
  if (has_schema && is.finite(B) && B >= 2 && bootstrap_primary) {
    numeric_data <- lapply(
      group[numeric_fields], function(value) suppressWarnings(as.numeric(value))
    )
    valid_logical <- if (is.logical(group[[logical_field]])) {
      !is.na(group[[logical_field]]) & group[[logical_field]]
    } else {
      rep(FALSE, nrow(group))
    }
    valid_character <- all(vapply(
      group[character_fields], is.character, logical(1L)
    ))
    valid_rows <- Reduce(`&`, lapply(numeric_data, is.finite)) &
      numeric_data$variance_weight_relearn_bootstrap > 0 &
      numeric_data$se_weight_relearn_bootstrap > 0 &
      numeric_data$se_fixed_weight_bootstrap > 0 &
      numeric_data$weight_uncertainty_ratio > 0 &
      numeric_data$weight_relearn_n_bootstrap == B &
      numeric_data$weight_bootstrap_seed >= 1 &
      numeric_data$weight_bootstrap_seed ==
        floor(numeric_data$weight_bootstrap_seed) &
      numeric_data$weight_bootstrap_failures >= 0 &
      numeric_data$weight_bootstrap_failures <= B &
      numeric_data$weight_bootstrap_failures ==
        floor(numeric_data$weight_bootstrap_failures) &
      as.character(group$weight_bootstrap_multiplier) %in%
        c("exponential", "site_stratified_exponential") &
      as.character(group$weight_bootstrap_screening_rule) == "soft_penalty" &
      valid_logical
    if (!valid_character) valid_rows[] <- FALSE
    if (any(valid_rows)) {
      variance_identity_error <- max(abs(
        numeric_data$variance_weight_relearn_bootstrap[valid_rows] -
          numeric_data$se_weight_relearn_bootstrap[valid_rows]^2
      ))
      ratio_identity_error <- max(abs(
        numeric_data$weight_uncertainty_ratio[valid_rows] -
          numeric_data$se_weight_relearn_bootstrap[valid_rows] /
            numeric_data$se_fixed_weight_bootstrap[valid_rows]
      ))
    }
  }

  invalid <- !valid_config || (any_schema && !has_schema) ||
    (is.finite(B) && B == 0 && any(observed)) ||
    (is.finite(B) && B >= 2 && bootstrap_primary &&
       (!has_schema || !all(valid_rows) ||
          !is.finite(variance_identity_error) ||
          variance_identity_error > identity_tolerance ||
          !is.finite(ratio_identity_error) ||
          ratio_identity_error > identity_tolerance)) ||
    (is.finite(B) && B >= 2 && !bootstrap_primary && any(observed))

  relearn_se <- if (has_schema) {
    suppressWarnings(as.numeric(group$se_weight_relearn_bootstrap))
  } else numeric(0)
  fixed_se <- if (has_schema) {
    suppressWarnings(as.numeric(group$se_fixed_weight_bootstrap))
  } else numeric(0)
  ratios <- if (has_schema) {
    suppressWarnings(as.numeric(group$weight_uncertainty_ratio))
  } else numeric(0)
  failures <- if (has_schema) {
    suppressWarnings(as.numeric(group$weight_bootstrap_failures))
  } else numeric(0)
  finite_mean <- function(value) {
    value <- value[is.finite(value)]
    if (length(value) > 0L) mean(value) else NA_real_
  }
  list(
    invalid_weight_relearn_bootstrap_diagnostics = invalid,
    weight_relearn_bootstrap_replications = sum(valid_rows),
    mean_se_weight_relearn_bootstrap = finite_mean(relearn_se[valid_rows]),
    mean_se_fixed_weight_bootstrap = finite_mean(fixed_se[valid_rows]),
    mean_weight_uncertainty_ratio = finite_mean(ratios[valid_rows]),
    weight_bootstrap_failures_total = sum(failures[valid_rows], na.rm = TRUE),
    weight_bootstrap_variance_identity_error = variance_identity_error,
    weight_bootstrap_ratio_identity_error = ratio_identity_error
  )
}

.summarize_simulation_diagnostic_group <- function(
    group, expected_replications, nominal_coverage, mc_level,
    se_ratio_limits, identity_tolerance) {
  coverage <- .simulation_coverage_as_logical(group$coverage)
  finite_estimate <- is.finite(group$estimate)
  finite_bias <- is.finite(group$bias)
  finite_se <- is.finite(group$se)
  finite_coverage <- !is.na(coverage)
  has_truth <- "truth" %in% names(group)
  finite_truth <- if (has_truth) {
    is.finite(group$truth)
  } else {
    rep(TRUE, nrow(group))
  }
  has_ci <- all(c("ci_lower", "ci_upper") %in% names(group))
  finite_ci <- if (has_ci) {
    is.finite(group$ci_lower) & is.finite(group$ci_upper)
  } else {
    rep(TRUE, nrow(group))
  }
  ordered_ci <- if (has_ci) {
    finite_ci & group$ci_lower <= group$ci_upper
  } else {
    rep(TRUE, nrow(group))
  }
  valid_replication_id <- if ("sim_id" %in% names(group)) {
    is.finite(group$sim_id) & group$sim_id >= 1 &
      group$sim_id == floor(group$sim_id)
  } else {
    rep(TRUE, nrow(group))
  }
  complete <- finite_estimate & finite_bias & finite_se & finite_coverage &
    finite_truth & finite_ci & ordered_ci & valid_replication_id

  bias <- group$bias[complete]
  se <- group$se[complete]
  coverage <- coverage[complete]
  n_complete <- length(bias)
  n_unique <- if ("sim_id" %in% names(group)) {
    length(unique(group$sim_id[!is.na(group$sim_id)]))
  } else {
    nrow(group)
  }
  n_duplicate <- if ("sim_id" %in% names(group)) {
    sum(duplicated(group$sim_id[!is.na(group$sim_id)]))
  } else {
    0L
  }
  n_invalid_replication_id <- sum(!valid_replication_id)

  bias_mean <- if (n_complete > 0L) mean(bias) else NA_real_
  empirical_sd <- if (n_complete > 1L) stats::sd(bias) else NA_real_
  rmse <- if (n_complete > 0L) sqrt(mean(bias^2)) else NA_real_
  mean_se <- if (n_complete > 0L) mean(se) else NA_real_
  coverage_rate <- if (n_complete > 0L) mean(coverage) else NA_real_
  truth_below_interval_fraction <- if (
      has_ci && has_truth && n_complete > 0L) {
    mean(group$truth[complete] < group$ci_lower[complete])
  } else {
    NA_real_
  }
  truth_above_interval_fraction <- if (
      has_ci && has_truth && n_complete > 0L) {
    mean(group$truth[complete] > group$ci_upper[complete])
  } else {
    NA_real_
  }
  bias_mcse <- if (n_complete > 1L) empirical_sd / sqrt(n_complete) else NA_real_
  coverage_mcse <- if (n_complete > 0L) {
    sqrt(coverage_rate * (1 - coverage_rate) / n_complete)
  } else {
    NA_real_
  }
  nominal_mcse <- if (n_complete > 0L) {
    sqrt(nominal_coverage * (1 - nominal_coverage) / n_complete)
  } else {
    NA_real_
  }
  mc_multiplier <- stats::qnorm(1 - (1 - mc_level) / 2)
  coverage_mc_lower <- max(0, nominal_coverage - mc_multiplier * nominal_mcse)
  coverage_mc_upper <- min(1, nominal_coverage + mc_multiplier * nominal_mcse)
  se_to_empirical_sd <- if (is.finite(empirical_sd) && empirical_sd > 0) {
    mean_se / empirical_sd
  } else {
    NA_real_
  }
  centered_bias <- bias - bias_mean
  second_central_moment <- if (n_complete > 1L) {
    mean(centered_bias^2)
  } else {
    NA_real_
  }
  bias_skewness <- if (is.finite(second_central_moment) &&
      second_central_moment > 0) {
    mean(centered_bias^3) / second_central_moment^(3 / 2)
  } else {
    NA_real_
  }
  bias_excess_kurtosis <- if (is.finite(second_central_moment) &&
      second_central_moment > 0) {
    mean(centered_bias^4) / second_central_moment^2 - 3
  } else {
    NA_real_
  }
  normal_reference_coverage <- if (
      is.finite(empirical_sd) && empirical_sd > 0 && is.finite(mean_se)
  ) {
    z <- stats::qnorm(0.5 + nominal_coverage / 2)
    stats::pnorm((z * mean_se - bias_mean) / empirical_sd) -
      stats::pnorm((-z * mean_se - bias_mean) / empirical_sd)
  } else {
    NA_real_
  }
  rmse_from_moments <- if (n_complete > 1L) {
    sqrt(bias_mean^2 + ((n_complete - 1) / n_complete) * empirical_sd^2)
  } else if (n_complete == 1L) {
    abs(bias_mean)
  } else {
    NA_real_
  }
  rmse_identity_error <- abs(rmse - rmse_from_moments)
  nuisance_diagnostics <- .summarize_group_nuisance_diagnostics(group)
  method_values <- if ("method" %in% names(group)) {
    unique(as.character(group$method))
  } else {
    character(0)
  }
  weight_relearn_bootstrap_diagnostics <-
    .summarize_weight_relearn_bootstrap_diagnostics(
      group, method_values, identity_tolerance
    )
  direct_tate_diagnostics_required <-
    length(method_values) == 1L &&
    method_values %in% c(
      "one_round_crossfit_ate", "two_round_crossfit_ate",
      "one_round_crossfit_ate_hard_threshold",
      "two_round_crossfit_ate_hard_threshold",
      "one_round_crossfit_ate_quadratic_bias",
      "two_round_crossfit_ate_quadratic_bias"
    )
  direct_aggregation_columns <-
    .direct_tate_aggregation_diagnostic_columns()
  complete_direct_aggregation_diagnostics <-
    all(direct_aggregation_columns %in% names(group)) &&
    all(vapply(group[direct_aggregation_columns], function(values) {
      values <- suppressWarnings(as.numeric(values))
      length(values) == nrow(group) && all(is.finite(values))
    }, logical(1L)))
  invalid_direct_tate_diagnostics <-
    direct_tate_diagnostics_required &&
    !complete_direct_aggregation_diagnostics
  if (direct_tate_diagnostics_required &&
      complete_direct_aggregation_diagnostics) {
    invalid_direct_tate_diagnostics <-
      any(group$mean_abs_source_weight < 0) ||
      any(group$max_abs_source_weight < group$mean_abs_source_weight) ||
      any(group$max_wald_statistic < group$mean_wald_statistic) ||
      any(group$mean_wald_statistic < 0) ||
      any(group$penalized_source_fold_fraction < 0 |
            group$penalized_source_fold_fraction > 1) ||
      any(group$max_weight_optimizer_iterations < 1 |
            group$max_weight_optimizer_iterations !=
              floor(group$max_weight_optimizer_iterations)) ||
      any(group$max_weight_psd_ridge < 0) ||
      any(group$weight_psd_ridge_fold_fraction < 0 |
            group$weight_psd_ridge_fold_fraction > 1) ||
      any(group$se_fixed_weights < 0) ||
      any(group$weight_layer_indirect_variance < 0) ||
      any(group$weight_layer_kink_cells < 0 |
            group$weight_layer_kink_cells !=
              floor(group$weight_layer_kink_cells)) ||
      any(group$inference_logit_truncated < 0 |
            group$inference_logit_truncated !=
              floor(group$inference_logit_truncated)) ||
      any(group$inference_logit_truncation_fraction < 0 |
            group$inference_logit_truncation_fraction > 1) ||
      any(group$inference_max_abs_logit < 0) ||
      any(group$inference_safety_clip_count < 0 |
            group$inference_safety_clip_count !=
              floor(group$inference_safety_clip_count))

    hard_threshold_row <- grepl("_ate_hard_threshold$", method_values)
    if (hard_threshold_row) {
      inclusion_columns <- grep(
        "^source_.+_inclusion_fraction$", names(group), value = TRUE
      )
      fold_inclusion_columns <- grep(
        "^source_.+_fold_[0-9]+_included$", names(group), value = TRUE
      )
      fold_weight_columns <- grep(
        "^source_.+_fold_[0-9]+_weight$", names(group), value = TRUE
      )
      source_fold_key <- function(column, suffix) sub(suffix, "", column)
      inclusion_keys <- source_fold_key(
        fold_inclusion_columns, "_included$"
      )
      weight_keys <- source_fold_key(fold_weight_columns, "_weight$")
      hard_schema_ok <- length(inclusion_columns) > 0L &&
        length(fold_inclusion_columns) > 0L &&
        setequal(inclusion_keys, weight_keys)
      if (hard_schema_ok) {
        inclusion_values <- suppressWarnings(as.numeric(unlist(
          group[c(inclusion_columns, fold_inclusion_columns)],
          use.names = FALSE
        )))
        hard_schema_ok <- all(is.finite(inclusion_values)) &&
          all(inclusion_values >= 0 & inclusion_values <= 1)
      }
      excluded_weight_ok <- hard_schema_ok
      if (excluded_weight_ok) {
        for (column in fold_inclusion_columns) {
          key <- source_fold_key(column, "_included$")
          weight_column <- paste0(key, "_weight")
          included <- suppressWarnings(as.numeric(group[[column]]))
          weights <- suppressWarnings(as.numeric(group[[weight_column]]))
          if (any(included == 0 & weights != 0)) {
            excluded_weight_ok <- FALSE
            break
          }
        }
      }
      invalid_direct_tate_diagnostics <-
        invalid_direct_tate_diagnostics || !hard_schema_ok ||
        !excluded_weight_ok
    }
  }
  max_abs_source_weight <- if ("max_abs_source_weight" %in% names(group)) {
    values <- suppressWarnings(as.numeric(group$max_abs_source_weight))
    values <- values[is.finite(values)]
    if (length(values) > 0L) max(values) else NA_real_
  } else {
    NA_real_
  }
  observed_weight_iterations <- if (
    "max_weight_optimizer_iterations" %in% names(group)
  ) {
    suppressWarnings(as.numeric(group$max_weight_optimizer_iterations))
  } else {
    numeric(0)
  }
  observed_weight_ridges <- if ("max_weight_psd_ridge" %in% names(group)) {
    suppressWarnings(as.numeric(group$max_weight_psd_ridge))
  } else {
    numeric(0)
  }
  finite_weight_iterations <- observed_weight_iterations[
    is.finite(observed_weight_iterations)
  ]
  finite_weight_ridges <- observed_weight_ridges[is.finite(observed_weight_ridges)]
  max_weight_optimizer_iterations <- if (length(finite_weight_iterations) > 0L) {
    max(finite_weight_iterations)
  } else {
    NA_real_
  }
  max_weight_psd_ridge <- if (length(finite_weight_ridges) > 0L) {
    max(finite_weight_ridges)
  } else {
    NA_real_
  }
  weight_psd_ridge_replications <- if (length(finite_weight_ridges) > 0L) {
    sum(finite_weight_ridges > 0)
  } else {
    NA_integer_
  }
  invalid_weight_optimizer_diagnostics <-
    any(!is.na(observed_weight_iterations) &
          (!is.finite(observed_weight_iterations) | observed_weight_iterations < 1)) ||
    any(!is.na(observed_weight_ridges) &
          (!is.finite(observed_weight_ridges) | observed_weight_ridges < 0))
  truncation_fraction <- if (
    "inference_logit_truncation_fraction" %in% names(group)
  ) {
    values <- suppressWarnings(as.numeric(
      group$inference_logit_truncation_fraction
    ))
    values <- values[is.finite(values)]
    if (length(values) > 0L) mean(values) else NA_real_
  } else {
    NA_real_
  }
  max_inference_abs_logit <- if (
    "inference_max_abs_logit" %in% names(group)
  ) {
    values <- suppressWarnings(as.numeric(group$inference_max_abs_logit))
    values <- values[is.finite(values)]
    if (length(values) > 0L) max(values) else NA_real_
  } else {
    NA_real_
  }
  inference_safety_clip_total <- if (
    "inference_safety_clip_count" %in% names(group)
  ) {
    values <- suppressWarnings(as.numeric(group$inference_safety_clip_count))
    values <- values[is.finite(values)]
    if (length(values) > 0L) sum(values) else NA_real_
  } else {
    NA_real_
  }
  numeric_column <- function(column) {
    # Preserve one value per replicate when optional metadata are absent.
    # Returning numeric(0) here would let method-specific fail-closed checks
    # silently collapse to length zero (notably for DR clipping metadata).
    if (!column %in% names(group)) return(rep(NA_real_, nrow(group)))
    suppressWarnings(as.numeric(group[[column]]))
  }
  finite_column <- function(column) {
    values <- numeric_column(column)
    values[is.finite(values)]
  }
  min_site_cells <- finite_column("min_site_arm_outcome_cell_n")
  min_target_cells <- finite_column("min_target_arm_outcome_cell_n")
  sparse_cell_counts <- finite_column("n_site_arm_outcome_cells_below_8")
  min_site_arm_outcome_cell_n_observed <- if (length(min_site_cells)) {
    min(min_site_cells)
  } else {
    NA_real_
  }
  min_target_arm_outcome_cell_n_observed <- if (length(min_target_cells)) {
    min(min_target_cells)
  } else {
    NA_real_
  }
  sparse_binary_cell_replications <- if (length(sparse_cell_counts)) {
    sum(sparse_cell_counts > 0)
  } else {
    NA_integer_
  }
  cell_diagnostic_columns <- c(
    "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
    "n_site_arm_outcome_cells_below_8"
  )
  cell_diagnostic_matrix <- as.matrix(data.frame(lapply(
    cell_diagnostic_columns,
    function(column) {
      values <- numeric_column(column)
      if (length(values) == 0L) rep(NA_real_, nrow(group)) else values
    }
  )))
  cell_diagnostic_values <- as.numeric(cell_diagnostic_matrix)
  observed_cell_diagnostics <- !is.na(cell_diagnostic_values)
  outcome_families <- if ("outcome_family" %in% names(group)) {
    unique(as.character(group$outcome_family))
  } else {
    character(0)
  }
  binary_cell_diagnostics_required <-
    length(outcome_families) == 1L && identical(outcome_families, "binomial")
  max_sparse_cells <- if (
      "K" %in% names(group) && length(unique(group$K)) == 1L &&
      is.finite(unique(group$K))) {
    4 * (as.integer(unique(group$K)) + 1L)
  } else {
    NA_integer_
  }
  invalid_cell_diagnostics <-
    (binary_cell_diagnostics_required && any(!observed_cell_diagnostics)) ||
    any(
    observed_cell_diagnostics &
      (!is.finite(cell_diagnostic_values) |
         cell_diagnostic_values < 0 |
         cell_diagnostic_values != floor(cell_diagnostic_values))
    ) ||
    any(
      is.finite(cell_diagnostic_matrix[, 1L]) &
        is.finite(cell_diagnostic_matrix[, 2L]) &
        cell_diagnostic_matrix[, 1L] > cell_diagnostic_matrix[, 2L]
    ) ||
    (is.finite(max_sparse_cells) && any(
      is.finite(cell_diagnostic_matrix[, 3L]) &
        cell_diagnostic_matrix[, 3L] > max_sparse_cells
    ))

  dr_fractions <- numeric_column("dr_weight_fraction_clipped")
  dr_max_site_fractions <- numeric_column(
    "dr_weight_max_site_fraction_clipped"
  )
  dr_n <- numeric_column("dr_weight_n")
  dr_n_clipped <- numeric_column("dr_weight_n_clipped")
  dr_min_before_clipping_raw <- numeric_column(
    "dr_weight_min_before_clipping"
  )
  dr_max_before_clipping_raw <- numeric_column(
    "dr_weight_max_before_clipping"
  )
  dr_min_before_clipping <- dr_min_before_clipping_raw[
    is.finite(dr_min_before_clipping_raw)
  ]
  dr_max_before_clipping <- dr_max_before_clipping_raw[
    is.finite(dr_max_before_clipping_raw)
  ]
  dr_diagnostic_matrix <- cbind(
    dr_n, dr_n_clipped, dr_fractions, dr_max_site_fractions,
    dr_min_before_clipping_raw, dr_max_before_clipping_raw
  )
  dr_methods <- c(
    "federated_dr", "pooled_dr", "federated_dr_ate", "pooled_dr_ate"
  )
  dr_diagnostics_required <-
    length(method_values) == 1L && method_values %in% dr_methods
  dr_rows_with_any_metadata <- rowSums(!is.na(dr_diagnostic_matrix)) > 0L
  dr_rows_requiring_complete_metadata <-
    dr_rows_with_any_metadata | dr_diagnostics_required
  valid_dr_rows <- is.finite(dr_fractions)
  mean_dr_weight_fraction_clipped <- if (any(valid_dr_rows)) {
    mean(dr_fractions[valid_dr_rows])
  } else {
    NA_real_
  }
  max_dr_weight_site_fraction_clipped <- if (
      any(is.finite(dr_max_site_fractions))) {
    max(dr_max_site_fractions[is.finite(dr_max_site_fractions)])
  } else {
    NA_real_
  }
  dr_weight_clipped_replications <- if (any(valid_dr_rows)) {
    sum(dr_fractions[valid_dr_rows] > 0)
  } else {
    NA_integer_
  }
  min_dr_weight_before_clipping <- if (length(dr_min_before_clipping)) {
    min(dr_min_before_clipping)
  } else {
    NA_real_
  }
  max_dr_weight_before_clipping <- if (length(dr_max_before_clipping)) {
    max(dr_max_before_clipping)
  } else {
    NA_real_
  }
  invalid_observed <- function(values, predicate) {
    observed <- !is.na(values)
    any(observed & (!is.finite(values) | predicate(values)))
  }
  paired_dr_counts <- !is.na(dr_n) & !is.na(dr_n_clipped)
  paired_dr_fractions <- !is.na(dr_fractions) &
    !is.na(dr_max_site_fractions)
  invalid_dr_weight_diagnostics <-
    any(dr_rows_requiring_complete_metadata &
          rowSums(is.finite(dr_diagnostic_matrix)) !=
            ncol(dr_diagnostic_matrix)) ||
    invalid_observed(dr_fractions, function(values) {
      values < 0 | values > 1
    }) ||
    invalid_observed(dr_max_site_fractions, function(values) {
      values < 0 | values > 1
    }) ||
    invalid_observed(dr_n, function(values) {
      values < 1 | values != floor(values)
    }) ||
    invalid_observed(dr_n_clipped, function(values) {
      values < 0 | values != floor(values)
    }) ||
    any(paired_dr_counts & dr_n_clipped > dr_n) ||
    any(paired_dr_fractions & dr_fractions > dr_max_site_fractions) ||
    any(
      paired_dr_counts & !is.na(dr_fractions) &
        abs(dr_fractions - dr_n_clipped / dr_n) > identity_tolerance
    ) ||
    invalid_observed(dr_min_before_clipping_raw, function(values) {
      values <= 0
    }) ||
    invalid_observed(dr_max_before_clipping_raw, function(values) {
      values <= 0
    }) ||
    any(
      !is.na(dr_min_before_clipping_raw) &
        !is.na(dr_max_before_clipping_raw) &
        dr_min_before_clipping_raw > dr_max_before_clipping_raw
    )
  bootstrap_diagnostic_columns <- c(
    "comparison_variance_method", "comparison_se_analytic",
    "comparison_se_bootstrap", "comparison_n_bootstrap"
  )
  has_bootstrap_diagnostics <- all(
    bootstrap_diagnostic_columns %in% names(group)
  )
  bootstrap_rows <- if (has_bootstrap_diagnostics) {
    is.finite(group$comparison_se_analytic) |
      is.finite(group$comparison_se_bootstrap) |
      is.finite(group$comparison_n_bootstrap)
  } else {
    rep(FALSE, nrow(group))
  }
  bootstrap_methods <- c(
    "sample_size", "inverse_variance", "federated_dr", "pooled_dr",
    "sample_size_ate", "inverse_variance_ate", "federated_dr_ate",
    "pooled_dr_ate"
  )
  bootstrap_required <-
    length(method_values) == 1L && method_values %in% bootstrap_methods
  valid_bootstrap_rows <- if (has_bootstrap_diagnostics) {
    bootstrap_rows &
      is.finite(group$comparison_se_analytic) &
      is.finite(group$comparison_se_bootstrap) &
      is.finite(group$comparison_n_bootstrap) &
      group$comparison_se_analytic > 0 &
      group$comparison_se_bootstrap > 0 &
      group$comparison_n_bootstrap >= 2 &
      group$comparison_n_bootstrap == floor(group$comparison_n_bootstrap) &
      !is.na(group$comparison_variance_method) &
      group$comparison_variance_method == "bootstrap"
  } else {
    rep(FALSE, nrow(group))
  }
  mean_comparison_se_analytic <- if (any(valid_bootstrap_rows)) {
    mean(group$comparison_se_analytic[valid_bootstrap_rows])
  } else {
    NA_real_
  }
  mean_comparison_se_bootstrap <- if (any(valid_bootstrap_rows)) {
    mean(group$comparison_se_bootstrap[valid_bootstrap_rows])
  } else {
    NA_real_
  }
  bootstrap_to_analytic_se <- if (
      is.finite(mean_comparison_se_analytic) &&
      mean_comparison_se_analytic > 0) {
    mean_comparison_se_bootstrap / mean_comparison_se_analytic
  } else {
    NA_real_
  }
  bootstrap_report_error <- if (
      any(valid_bootstrap_rows) &&
      "comparison_variance_method" %in% names(group)) {
    reported_bootstrap <- valid_bootstrap_rows &
      group$comparison_variance_method == "bootstrap"
    if (any(reported_bootstrap)) {
      max(abs(
        group$se[reported_bootstrap] -
          group$comparison_se_bootstrap[reported_bootstrap]
      ))
    } else {
      0
    }
  } else {
    NA_real_
  }
  invalid_bootstrap_diagnostics <-
    (bootstrap_required && (!has_bootstrap_diagnostics ||
       !all(valid_bootstrap_rows))) ||
    (any(bootstrap_rows) &&
       (!all(valid_bootstrap_rows[bootstrap_rows]) ||
          (is.finite(bootstrap_report_error) &&
             bootstrap_report_error > identity_tolerance)))

  bias_identity_error <- NA_real_
  if ("truth" %in% names(group)) {
    valid_bias_identity <- finite_estimate & finite_bias & finite_truth
    if (any(valid_bias_identity)) {
      bias_identity_error <- max(abs(
        group$bias[valid_bias_identity] -
          (group$estimate[valid_bias_identity] -
             group$truth[valid_bias_identity])
      ))
    }
  }
  ci_arithmetic_error <- NA_real_
  if (has_ci) {
    valid_ci_arithmetic <- finite_estimate & finite_se & finite_ci
    if (any(valid_ci_arithmetic)) {
      z <- stats::qnorm(0.5 + nominal_coverage / 2)
      ci_errors <- c(
        group$ci_lower[valid_ci_arithmetic] -
          (group$estimate[valid_ci_arithmetic] -
             z * group$se[valid_ci_arithmetic]),
        group$ci_upper[valid_ci_arithmetic] -
          (group$estimate[valid_ci_arithmetic] +
             z * group$se[valid_ci_arithmetic])
      )
      if ("ci_width" %in% names(group)) {
        valid_width <- valid_ci_arithmetic & is.finite(group$ci_width)
        if (any(valid_width)) {
          ci_errors <- c(
            ci_errors,
            group$ci_width[valid_width] -
              (group$ci_upper[valid_width] - group$ci_lower[valid_width])
          )
        }
      }
      ci_arithmetic_error <- max(abs(ci_errors))
    }
  }

  coverage_inconsistent <- FALSE
  if (has_ci && "truth" %in% names(group)) {
    valid_ci <- finite_coverage & finite_truth & finite_ci & ordered_ci
    if (any(valid_ci)) {
      recomputed <- group$truth[valid_ci] >= group$ci_lower[valid_ci] &
        group$truth[valid_ci] <= group$ci_upper[valid_ci]
      coverage_inconsistent <- any(
        recomputed != .simulation_coverage_as_logical(group$coverage[valid_ci])
      )
    }
  }

  flags <- c(
    incomplete_replications = n_unique < expected_replications,
    excess_replications = n_unique > expected_replications,
    duplicated_replications = n_duplicate > 0L,
    invalid_replication_id = n_invalid_replication_id > 0L,
    nonfinite_values = n_complete < nrow(group),
    nonpositive_standard_error = any(finite_se & group$se <= 0),
    invalid_confidence_interval = has_ci && any(!ordered_ci),
    bias_truth_inconsistent = is.finite(bias_identity_error) &&
      bias_identity_error > identity_tolerance,
    ci_arithmetic_inconsistent = is.finite(ci_arithmetic_error) &&
      ci_arithmetic_error > identity_tolerance,
    coverage_ci_inconsistent = coverage_inconsistent,
    rmse_identity_failed = is.finite(rmse_identity_error) &&
      rmse_identity_error > identity_tolerance,
    invalid_nuisance_diagnostics =
      isTRUE(as.logical(nuisance_diagnostics[[
        "invalid_nuisance_diagnostics"
      ]])),
    nuisance_nonconvergence_detected =
      nuisance_diagnostics[["nuisance_nonconverged_fits_total"]] > 0,
    nuisance_line_search_failures_detected =
      sum(c(
        nuisance_diagnostics[[
          "face_initial_dr_line_search_failures_total"
        ]],
        nuisance_diagnostics[[
          "face_calibrated_dr_line_search_failures_total"
        ]],
        nuisance_diagnostics[[
          "face_calibrated_outcome_line_search_failures_total"
        ]]
      ), na.rm = TRUE) > 0,
    nuisance_cv_candidates_excluded =
      sum(c(
        nuisance_diagnostics[[
          "face_initial_dr_cv_invalid_fold_fits_total"
        ]],
        nuisance_diagnostics[[
          "face_calibrated_dr_cv_invalid_fold_fits_total"
        ]],
        nuisance_diagnostics[[
          "face_calibrated_outcome_cv_invalid_fold_fits_total"
        ]]
      ), na.rm = TRUE) > 0,
    invalid_weight_optimizer_diagnostics =
      invalid_weight_optimizer_diagnostics,
    invalid_direct_tate_diagnostics = invalid_direct_tate_diagnostics,
    invalid_bootstrap_diagnostics = invalid_bootstrap_diagnostics,
    invalid_weight_relearn_bootstrap_diagnostics =
      isTRUE(weight_relearn_bootstrap_diagnostics[[
        "invalid_weight_relearn_bootstrap_diagnostics"
      ]]),
    weight_relearn_bootstrap_failures_detected =
      weight_relearn_bootstrap_diagnostics[[
        "weight_bootstrap_failures_total"
      ]] > 0,
    invalid_cell_diagnostics = invalid_cell_diagnostics,
    invalid_dr_weight_diagnostics = invalid_dr_weight_diagnostics,
    sparse_binary_cells_detected =
      is.finite(sparse_binary_cell_replications) &&
      sparse_binary_cell_replications > 0,
    density_ratio_clipping_detected =
      is.finite(dr_weight_clipped_replications) &&
      dr_weight_clipped_replications > 0,
    inference_safety_clipping_detected =
      is.finite(inference_safety_clip_total) &&
      inference_safety_clip_total > 0,
    extreme_aggregation_weight = is.finite(max_abs_source_weight) &&
      max_abs_source_weight > 10,
    coverage_below_mc_band = is.finite(coverage_rate) &&
      coverage_rate < coverage_mc_lower,
    coverage_above_mc_band = is.finite(coverage_rate) &&
      coverage_rate > coverage_mc_upper,
    se_empirical_ratio_low = is.finite(se_to_empirical_sd) &&
      se_to_empirical_sd < se_ratio_limits[[1L]],
    se_empirical_ratio_high = is.finite(se_to_empirical_sd) &&
      se_to_empirical_sd > se_ratio_limits[[2L]],
    bias_exceeds_2_mcse = is.finite(bias_mcse) && bias_mcse > 0 &&
      abs(bias_mean) > 2 * bias_mcse
  )
  coverage_diagnosis <- if (!is.finite(coverage_rate)) {
    "not_evaluable"
  } else if (isTRUE(flags[["coverage_below_mc_band"]])) {
    biased <- isTRUE(flags[["bias_exceeds_2_mcse"]])
    se_low <- isTRUE(flags[["se_empirical_ratio_low"]])
    se_high <- isTRUE(flags[["se_empirical_ratio_high"]])
    if (biased && se_low) {
      "bias_and_se_underestimation"
    } else if (biased && se_high) {
      "centering_bias_despite_conservative_se"
    } else if (biased) {
      "primarily_centering_bias"
    } else if (se_low) {
      "primarily_se_underestimation"
    } else {
      "unresolved_finite_sample_shape"
    }
  } else if (isTRUE(flags[["coverage_above_mc_band"]])) {
    if (isTRUE(flags[["se_empirical_ratio_high"]])) {
      "primarily_se_overestimation"
    } else {
      "unresolved_finite_sample_shape"
    }
  } else {
    "within_nominal_mc_band"
  }

  data.frame(
    n_rows = nrow(group),
    n_unique_replications = n_unique,
    n_complete_replications = n_complete,
    n_duplicate_replications = n_duplicate,
    n_invalid_replication_ids = n_invalid_replication_id,
    expected_replications = expected_replications,
    bias = bias_mean,
    bias_mcse = bias_mcse,
    empirical_sd = empirical_sd,
    bias_skewness = bias_skewness,
    bias_excess_kurtosis = bias_excess_kurtosis,
    mean_se = mean_se,
    se_to_empirical_sd = se_to_empirical_sd,
    rmse = rmse,
    rmse_from_moments = rmse_from_moments,
    rmse_identity_error = rmse_identity_error,
    bias_identity_error = bias_identity_error,
    ci_arithmetic_error = ci_arithmetic_error,
    coverage = coverage_rate,
    truth_below_interval_fraction = truth_below_interval_fraction,
    truth_above_interval_fraction = truth_above_interval_fraction,
    coverage_mcse = coverage_mcse,
    coverage_mc_lower = coverage_mc_lower,
    coverage_mc_upper = coverage_mc_upper,
    normal_reference_coverage = normal_reference_coverage,
    max_abs_source_weight_observed = max_abs_source_weight,
    max_weight_optimizer_iterations_observed = max_weight_optimizer_iterations,
    max_weight_psd_ridge_observed = max_weight_psd_ridge,
    weight_psd_ridge_replications = weight_psd_ridge_replications,
    mean_inference_logit_truncation_fraction = truncation_fraction,
    max_inference_abs_logit_observed = max_inference_abs_logit,
    inference_safety_clip_total = inference_safety_clip_total,
    min_site_arm_outcome_cell_n_observed =
      min_site_arm_outcome_cell_n_observed,
    min_target_arm_outcome_cell_n_observed =
      min_target_arm_outcome_cell_n_observed,
    sparse_binary_cell_replications = sparse_binary_cell_replications,
    mean_dr_weight_fraction_clipped = mean_dr_weight_fraction_clipped,
    max_dr_weight_site_fraction_clipped =
      max_dr_weight_site_fraction_clipped,
    dr_weight_clipped_replications = dr_weight_clipped_replications,
    min_dr_weight_before_clipping = min_dr_weight_before_clipping,
    max_dr_weight_before_clipping = max_dr_weight_before_clipping,
    mean_comparison_se_analytic = mean_comparison_se_analytic,
    mean_comparison_se_bootstrap = mean_comparison_se_bootstrap,
    bootstrap_to_analytic_se = bootstrap_to_analytic_se,
    bootstrap_report_identity_error = bootstrap_report_error,
    weight_relearn_bootstrap_replications =
      weight_relearn_bootstrap_diagnostics[[
        "weight_relearn_bootstrap_replications"
      ]],
    mean_se_weight_relearn_bootstrap =
      weight_relearn_bootstrap_diagnostics[[
        "mean_se_weight_relearn_bootstrap"
      ]],
    mean_se_fixed_weight_bootstrap =
      weight_relearn_bootstrap_diagnostics[[
        "mean_se_fixed_weight_bootstrap"
      ]],
    mean_weight_uncertainty_ratio =
      weight_relearn_bootstrap_diagnostics[[
        "mean_weight_uncertainty_ratio"
      ]],
    weight_bootstrap_failures_total =
      weight_relearn_bootstrap_diagnostics[[
        "weight_bootstrap_failures_total"
      ]],
    weight_bootstrap_variance_identity_error =
      weight_relearn_bootstrap_diagnostics[[
        "weight_bootstrap_variance_identity_error"
      ]],
    weight_bootstrap_ratio_identity_error =
      weight_relearn_bootstrap_diagnostics[[
        "weight_bootstrap_ratio_identity_error"
      ]],
    nuisance_diagnostic_replications =
      nuisance_diagnostics[["nuisance_diagnostic_replications"]],
    nuisance_nonconverged_replications =
      nuisance_diagnostics[["nuisance_nonconverged_replications"]],
    nuisance_nonconverged_fits_total =
      nuisance_diagnostics[["nuisance_nonconverged_fits_total"]],
    face_initial_dr_nonconverged_total =
      nuisance_diagnostics[["face_initial_dr_nonconverged_total"]],
    face_calibrated_dr_nonconverged_total =
      nuisance_diagnostics[["face_calibrated_dr_nonconverged_total"]],
    face_calibrated_outcome_nonconverged_total =
      nuisance_diagnostics[["face_calibrated_outcome_nonconverged_total"]],
    face_initial_dr_cv_invalid_fold_fits_total =
      nuisance_diagnostics[[
        "face_initial_dr_cv_invalid_fold_fits_total"
      ]],
    face_initial_dr_cv_invalid_lambdas_total =
      nuisance_diagnostics[["face_initial_dr_cv_invalid_lambdas_total"]],
    face_initial_dr_cv_path_tail_skipped_fold_fits_total =
      nuisance_diagnostics[[
        "face_initial_dr_cv_path_tail_skipped_fold_fits_total"
      ]],
    face_calibrated_dr_cv_invalid_fold_fits_total =
      nuisance_diagnostics[[
        "face_calibrated_dr_cv_invalid_fold_fits_total"
      ]],
    face_calibrated_dr_cv_invalid_lambdas_total =
      nuisance_diagnostics[["face_calibrated_dr_cv_invalid_lambdas_total"]],
    face_calibrated_dr_cv_path_tail_skipped_fold_fits_total =
      nuisance_diagnostics[[
        "face_calibrated_dr_cv_path_tail_skipped_fold_fits_total"
      ]],
    face_calibrated_outcome_cv_invalid_fold_fits_total =
      nuisance_diagnostics[[
        "face_calibrated_outcome_cv_invalid_fold_fits_total"
      ]],
    face_calibrated_outcome_cv_invalid_lambdas_total =
      nuisance_diagnostics[[
        "face_calibrated_outcome_cv_invalid_lambdas_total"
      ]],
    face_calibrated_outcome_cv_path_tail_skipped_fold_fits_total =
      nuisance_diagnostics[[
        "face_calibrated_outcome_cv_path_tail_skipped_fold_fits_total"
      ]],
    face_initial_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_initial_dr_line_search_failures_total"
      ]],
    face_calibrated_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_calibrated_dr_line_search_failures_total"
      ]],
    face_calibrated_outcome_line_search_failures_total =
      nuisance_diagnostics[[
        "face_calibrated_outcome_line_search_failures_total"
      ]],
    face_initial_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_initial_dr_support_floor_applied_total"
      ]],
    face_calibrated_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_calibrated_dr_support_floor_applied_total"
      ]],
    face_mu1_initial_dr_nonconverged_total =
      nuisance_diagnostics[["face_mu1_initial_dr_nonconverged_total"]],
    face_mu0_initial_dr_nonconverged_total =
      nuisance_diagnostics[["face_mu0_initial_dr_nonconverged_total"]],
    face_mu1_calibrated_dr_nonconverged_total =
      nuisance_diagnostics[["face_mu1_calibrated_dr_nonconverged_total"]],
    face_mu0_calibrated_dr_nonconverged_total =
      nuisance_diagnostics[["face_mu0_calibrated_dr_nonconverged_total"]],
    face_mu1_calibrated_outcome_nonconverged_total =
      nuisance_diagnostics[[
        "face_mu1_calibrated_outcome_nonconverged_total"
      ]],
    face_mu0_calibrated_outcome_nonconverged_total =
      nuisance_diagnostics[[
        "face_mu0_calibrated_outcome_nonconverged_total"
      ]],
    face_mu1_initial_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu1_initial_dr_line_search_failures_total"
      ]],
    face_mu0_initial_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu0_initial_dr_line_search_failures_total"
      ]],
    face_mu1_calibrated_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu1_calibrated_dr_line_search_failures_total"
      ]],
    face_mu0_calibrated_dr_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu0_calibrated_dr_line_search_failures_total"
      ]],
    face_mu1_calibrated_outcome_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu1_calibrated_outcome_line_search_failures_total"
      ]],
    face_mu0_calibrated_outcome_line_search_failures_total =
      nuisance_diagnostics[[
        "face_mu0_calibrated_outcome_line_search_failures_total"
      ]],
    face_mu1_initial_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_mu1_initial_dr_support_floor_applied_total"
      ]],
    face_mu0_initial_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_mu0_initial_dr_support_floor_applied_total"
      ]],
    face_mu1_calibrated_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_mu1_calibrated_dr_support_floor_applied_total"
      ]],
    face_mu0_calibrated_dr_support_floor_applied_total =
      nuisance_diagnostics[[
        "face_mu0_calibrated_dr_support_floor_applied_total"
      ]],
    face_max_initial_dr_iterations_observed =
      nuisance_diagnostics[["face_max_initial_dr_iterations_observed"]],
    face_max_calibrated_dr_iterations_observed =
      nuisance_diagnostics[["face_max_calibrated_dr_iterations_observed"]],
    face_max_calibrated_outcome_iterations_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_outcome_iterations_observed"
      ]],
    face_max_initial_dr_update_ratio_observed =
      nuisance_diagnostics[[
        "face_max_initial_dr_update_ratio_observed"
      ]],
    face_max_initial_dr_abs_coefficient_observed =
      nuisance_diagnostics[[
        "face_max_initial_dr_abs_coefficient_observed"
      ]],
    face_max_calibrated_dr_update_ratio_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_dr_update_ratio_observed"
      ]],
    face_max_calibrated_dr_abs_coefficient_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_dr_abs_coefficient_observed"
      ]],
    face_max_calibrated_outcome_update_ratio_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_outcome_update_ratio_observed"
      ]],
    face_max_calibrated_outcome_abs_coefficient_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_outcome_abs_coefficient_observed"
      ]],
    face_max_initial_dr_support_floor_observed =
      nuisance_diagnostics[[
        "face_max_initial_dr_support_floor_observed"
      ]],
    face_max_calibrated_dr_support_floor_observed =
      nuisance_diagnostics[[
        "face_max_calibrated_dr_support_floor_observed"
      ]],
    coverage_diagnosis = coverage_diagnosis,
    diagnostic_status = .collapse_diagnostic_flags(flags),
    stringsAsFactors = FALSE
  )
}

#' Diagnose replicate-level simulation results
#'
#' Internal quality-control helper used by the MSI aggregation scripts.  It
#' checks every setting-method combination without dropping non-finite rows.
#'
#' @param results Replicate-level simulation result data frame.
#' @param expected_replications Intended number of Monte Carlo replications.
#' @param nominal_coverage Nominal confidence interval coverage.
#' @param mc_level Confidence level of the nominal Monte Carlo band.
#' @param se_ratio_limits Lower and upper review thresholds for mean SE divided
#'   by the empirical Monte Carlo standard deviation.
#' @param identity_tolerance Numerical tolerance for the RMSE moment identity.
#' @return One row per setting and method with diagnostic metrics and flags.
diagnose_simulation_results <- function(
    results, expected_replications = 200L, nominal_coverage = 0.95,
    mc_level = 0.95, se_ratio_limits = c(0.85, 1.15),
    identity_tolerance = 1e-10) {
  required <- c("method", "estimate", "bias", "se", "coverage")
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop(
      "simulation results are missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(expected_replications) != 1L ||
      !is.numeric(expected_replications) ||
      is.na(expected_replications) || !is.finite(expected_replications) ||
      expected_replications < 1L ||
      expected_replications != floor(expected_replications)) {
    stop("expected_replications must be one positive integer.", call. = FALSE)
  }
  if (length(nominal_coverage) != 1L || !is.finite(nominal_coverage) ||
      nominal_coverage <= 0 || nominal_coverage >= 1) {
    stop("nominal_coverage must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(mc_level) != 1L || !is.finite(mc_level) ||
      mc_level <= 0 || mc_level >= 1) {
    stop("mc_level must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(se_ratio_limits) != 2L ||
      any(!is.finite(se_ratio_limits)) ||
      se_ratio_limits[[1L]] <= 0 ||
      se_ratio_limits[[1L]] >= se_ratio_limits[[2L]]) {
    stop("se_ratio_limits must be two increasing positive values.",
         call. = FALSE)
  }
  numeric_columns <- intersect(.simulation_numeric_columns(), names(results))
  nonnumeric_columns <- numeric_columns[!vapply(
    results[numeric_columns], is.numeric, logical(1L)
  )]
  if (length(nonnumeric_columns) > 0L) {
    stop(
      "simulation numeric columns have nonnumeric storage: ",
      paste(nonnumeric_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(identity_tolerance) != 1L ||
      !is.numeric(identity_tolerance) || !is.finite(identity_tolerance) ||
      identity_tolerance < 0) {
    stop("identity_tolerance must be one finite non-negative number.",
         call. = FALSE)
  }

  # Older task CSVs predate the explicit estimand label.  Infer it from the
  # stable method suffix so retained results remain diagnosable and are never
  # silently dropped by interaction() because of a missing/NA group value.
  inferred_scope <- ifelse(
    .is_tate_method(results$method),
    "tate",
    "treated_mean"
  )
  if (!"estimand_scope" %in% names(results)) {
    results$estimand_scope <- inferred_scope
  } else {
    missing_scope <- is.na(results$estimand_scope) |
      !nzchar(trimws(as.character(results$estimand_scope)))
    results$estimand_scope[missing_scope] <- inferred_scope[missing_scope]
  }

  group_columns <- .simulation_diagnostic_group_columns(results)
  if (!"method" %in% group_columns) {
    stop("method must be a simulation diagnostic grouping column.",
         call. = FALSE)
  }
  group_key <- .simulation_group_key(results, group_columns)
  diagnostics <- lapply(split(results, group_key), function(group) {
    identifiers <- group[1L, group_columns, drop = FALSE]
    metrics <- .summarize_simulation_diagnostic_group(
      group = group,
      expected_replications = as.integer(expected_replications),
      nominal_coverage = nominal_coverage,
      mc_level = mc_level,
      se_ratio_limits = se_ratio_limits,
      identity_tolerance = identity_tolerance
    )
    cbind(identifiers, metrics)
  })
  result <- do.call(rbind, diagnostics)
  rownames(result) <- NULL
  result
}

#' Add within-setting RMSE comparisons to simulation diagnostics
#'
#' @param diagnostics Output from \code{diagnose_simulation_results}.
#' @param target_method Optional method label used for the target-only
#'   benchmark.  By default, TATE rows use \code{target_only_ate} and retained
#'   potential-outcome-mean rows use \code{target_only}.
#' @return Diagnostics augmented with RMSE rank and target-relative RMSE.
add_simulation_rmse_comparisons <- function(
    diagnostics, target_method = NULL) {
  required <- c("method", "rmse")
  missing <- setdiff(required, names(diagnostics))
  if (length(missing) > 0L) {
    stop(
      "diagnostics are missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  setting_columns <- setdiff(
    .simulation_diagnostic_group_columns(diagnostics), "method"
  )
  setting_key <- if (length(setting_columns) > 0L) {
    .simulation_group_key(diagnostics, setting_columns)
  } else {
    factor(rep("all", nrow(diagnostics)))
  }
  augmented <- lapply(split(diagnostics, setting_key), function(setting) {
    setting_target_method <- target_method
    if (is.null(setting_target_method)) {
      setting_target_method <- if (
          "estimand_scope" %in% names(setting) &&
          identical(as.character(setting$estimand_scope[[1L]]), "treated_mean")
      ) {
        "target_only"
      } else {
        "target_only_ate"
      }
    }
    target_rmse <- setting$rmse[setting$method == setting_target_method]
    target_rmse <- if (length(target_rmse) == 1L) target_rmse else NA_real_
    setting$rmse_rank <- rank(
      setting$rmse, ties.method = "min", na.last = "keep"
    )
    setting$rmse_relative_to_target <- setting$rmse / target_rmse
    setting
  })
  result <- do.call(rbind, augmented)
  rownames(result) <- NULL
  result
}

#' Summarize paired squared-error differences across Monte Carlo replicates
#'
#' Internal checkpoint helper. Unlike a comparison of two separately rounded
#' RMSE values, this calculation preserves the common-random-number pairing and
#' supplies a Monte Carlo standard error for the MSE difference.
#'
#' @param results Replicate-level simulation results.
#' @param method Method whose squared error is compared with the benchmark.
#' @param benchmark_method Benchmark method label.
#' @return One row per simulation setting with the paired mean squared-error
#'   difference (method minus benchmark), its Monte Carlo standard error and
#'   standardized ratio, and the fraction of replicates on which \code{method}
#'   has smaller squared error.
summarize_paired_mse_differences <- function(
    results, method, benchmark_method) {
  required <- c("sim_id", "method", "bias")
  missing <- setdiff(required, names(results))
  if (length(missing) > 0L) {
    stop(
      "results are missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  labels <- c(method, benchmark_method)
  if (!is.character(method) || length(method) != 1L ||
      !is.character(benchmark_method) || length(benchmark_method) != 1L ||
      anyNA(labels) || any(!nzchar(trimws(labels))) ||
      any(labels != trimws(labels)) ||
      identical(method, benchmark_method)) {
    stop("method and benchmark_method must be distinct nonempty labels.",
         call. = FALSE)
  }
  if (!is.numeric(results$bias)) {
    stop("results$bias must be numeric.", call. = FALSE)
  }

  paired_rows <- results[results$method %in% labels, , drop = FALSE]
  if (!all(labels %in% paired_rows$method)) {
    stop("both requested methods must be present in results.", call. = FALSE)
  }
  if (anyNA(paired_rows$sim_id)) {
    stop("paired sim_id values must be nonmissing.", call. = FALSE)
  }
  if (!is.numeric(paired_rows$sim_id) ||
      any(!is.finite(paired_rows$sim_id)) ||
      any(paired_rows$sim_id < 1) ||
      any(paired_rows$sim_id != floor(paired_rows$sim_id))) {
    stop("paired sim_id values must be positive finite integers.",
         call. = FALSE)
  }
  setting_columns <- setdiff(
    .simulation_diagnostic_group_columns(paired_rows), "method"
  )
  setting_key <- if (length(setting_columns) > 0L) {
    .simulation_group_key(paired_rows, setting_columns)
  } else {
    factor(rep("all", nrow(paired_rows)))
  }

  summaries <- lapply(split(paired_rows, setting_key), function(setting) {
    focal <- setting[setting$method == method, , drop = FALSE]
    benchmark <- setting[
      setting$method == benchmark_method, , drop = FALSE
    ]
    if (nrow(focal) == 0L || nrow(benchmark) == 0L ||
        anyDuplicated(focal$sim_id) || anyDuplicated(benchmark$sim_id)) {
      stop(
        "each setting must contain one unique row per method and sim_id.",
        call. = FALSE
      )
    }
    benchmark_index <- match(focal$sim_id, benchmark$sim_id)
    if (nrow(focal) != nrow(benchmark) || anyNA(benchmark_index)) {
      stop("method and benchmark sim_id sets must match exactly.",
           call. = FALSE)
    }
    focal_bias <- focal$bias
    benchmark_bias <- benchmark$bias[benchmark_index]
    if (any(!is.finite(c(focal_bias, benchmark_bias)))) {
      stop("paired squared-error inputs must be finite.", call. = FALSE)
    }
    difference <- focal_bias^2 - benchmark_bias^2
    identifiers <- if (length(setting_columns) > 0L) {
      focal[1L, setting_columns, drop = FALSE]
    } else {
      data.frame(.setting = "all", stringsAsFactors = FALSE)
    }
    difference_mcse <- if (length(difference) > 1L) {
      stats::sd(difference) / sqrt(length(difference))
    } else {
      NA_real_
    }
    cbind(
      identifiers,
      data.frame(
        method = method,
        benchmark_method = benchmark_method,
        n_paired_replications = length(difference),
        mean_squared_error_difference = mean(difference),
        squared_error_difference_mcse = difference_mcse,
        squared_error_difference_z = if (
            is.finite(difference_mcse) && difference_mcse > 0) {
          mean(difference) / difference_mcse
        } else {
          NA_real_
        },
        method_lower_squared_error_fraction = mean(difference < 0),
        stringsAsFactors = FALSE
      )
    )
  })
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL
  result
}

#' Summarize TATE paired MSE against all standard benchmarks
#'
#' @param results Replicate-level simulation results.
#' @param method TATE method label.
#' @return Bound output from \code{summarize_paired_mse_differences}, with one
#'   set of setting rows for every standard TATE benchmark.
summarize_tate_paired_mse_benchmarks <- function(
    results, method = "one_round_crossfit_ate") {
  summaries <- lapply(.tate_benchmark_methods(), function(benchmark) {
    summarize_paired_mse_differences(
      results,
      method = method,
      benchmark_method = benchmark
    )
  })
  result <- do.call(rbind, summaries)
  rownames(result) <- NULL
  result
}
