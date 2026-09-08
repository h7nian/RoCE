#!/usr/bin/env Rscript

# Shared policy for classifying simulation diagnostics.  The package-level
# diagnostic table deliberately records numerical events without deciding
# whether they invalidate a result.  Production aggregation and checkpoint
# audits source this file so that the classification is defined once.

roce_append_diagnostic_status <- function(status, flag, selected) {
  if (length(status) != length(selected)) {
    stop("status and selected must have the same length.", call. = FALSE)
  }
  selected <- !is.na(selected) & selected
  status <- as.character(status)
  already_present <- grepl(
    paste0("(^|;)", flag, "($|;)"), status
  )
  already_present[is.na(already_present)] <- FALSE
  selected <- selected & !already_present
  empty <- is.na(status) | !nzchar(status) | status == "ok"
  status[selected & empty] <- flag
  status[selected & !empty] <- paste0(status[selected & !empty], ";", flag)
  status
}

roce_augment_simulation_qc <- function(
    diagnostics, parameter_bound, boundary_fraction = 0.9) {
  if (!is.data.frame(diagnostics)) {
    stop("diagnostics must be a data frame.", call. = FALSE)
  }
  if (length(parameter_bound) != 1L || !is.finite(parameter_bound) ||
      parameter_bound <= 0) {
    stop("parameter_bound must be one finite positive number.", call. = FALSE)
  }
  if (length(boundary_fraction) != 1L || !is.finite(boundary_fraction) ||
      boundary_fraction <= 0 || boundary_fraction >= 1) {
    stop("boundary_fraction must lie strictly between zero and one.",
         call. = FALSE)
  }

  density_ratio_coefficient_columns <- c(
    "face_max_initial_dr_abs_coefficient_observed",
    "face_max_calibrated_dr_abs_coefficient_observed"
  )
  outcome_coefficient_column <-
    "face_max_calibrated_outcome_abs_coefficient_observed"
  coefficient_columns <- c(
    density_ratio_coefficient_columns, outcome_coefficient_column
  )
  missing_columns <- setdiff(
    c("diagnostic_status", coefficient_columns), names(diagnostics)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "simulation diagnostics are missing QC columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  threshold <- boundary_fraction * parameter_bound
  coefficient_matrix <- as.matrix(data.frame(lapply(
    diagnostics[coefficient_columns],
    function(values) suppressWarnings(as.numeric(values))
  )))
  density_ratio_boundary_detected <- apply(
    coefficient_matrix[, density_ratio_coefficient_columns, drop = FALSE], 1L,
    function(values) any(is.finite(values) & values >= threshold)
  )
  outcome_boundary_detected <- apply(
    coefficient_matrix[, outcome_coefficient_column, drop = FALSE], 1L,
    function(values) any(is.finite(values) & values >= threshold)
  )
  diagnostics$nuisance_density_ratio_boundary_detected <-
    density_ratio_boundary_detected
  diagnostics$nuisance_density_ratio_boundary_threshold <- threshold
  diagnostics$diagnostic_status <- roce_append_diagnostic_status(
    diagnostics$diagnostic_status,
    "nuisance_density_ratio_boundary_detected",
    density_ratio_boundary_detected
  )
  diagnostics$nuisance_outcome_boundary_detected <-
    outcome_boundary_detected
  diagnostics$nuisance_outcome_boundary_threshold <- threshold
  diagnostics$diagnostic_status <- roce_append_diagnostic_status(
    diagnostics$diagnostic_status,
    "nuisance_outcome_boundary_detected",
    outcome_boundary_detected
  )
  diagnostics
}

roce_classify_simulation_qc <- function(diagnostics) {
  if (!is.data.frame(diagnostics) ||
      !"diagnostic_status" %in% names(diagnostics)) {
    stop("diagnostics must contain diagnostic_status.", call. = FALSE)
  }
  implementation_patterns <- c(
    "incomplete_replications", "excess_replications",
    "duplicated_replications", "nonfinite_values",
    "invalid_replication_id", "invalid_confidence_interval",
    "bias_truth_inconsistent", "ci_arithmetic_inconsistent",
    "nonpositive_standard_error", "coverage_ci_inconsistent",
    "rmse_identity_failed", "invalid_nuisance_diagnostics",
    "nuisance_nonconvergence_detected",
    "nuisance_density_ratio_boundary_detected",
    "nuisance_outcome_boundary_detected",
    "invalid_weight_optimizer_diagnostics",
    "invalid_direct_tate_diagnostics",
    "invalid_bootstrap_diagnostics",
    "invalid_weight_relearn_bootstrap_diagnostics",
    "weight_relearn_bootstrap_failures_detected",
    "invalid_cell_diagnostics",
    "invalid_dr_weight_diagnostics",
    "inference_safety_clipping_detected"
  )
  review_patterns <- c(
    "coverage_below_mc_band", "coverage_above_mc_band",
    "se_empirical_ratio_low", "se_empirical_ratio_high",
    "bias_exceeds_2_mcse", "extreme_aggregation_weight",
    "sparse_binary_cells_detected", "density_ratio_clipping_detected",
    "nuisance_cv_candidates_excluded",
    # This counter records rejected intermediate coordinate steps. A
    # converged fit necessarily ended on a complete scan with no unrecovered
    # line-search failure, so the historical count is a review signal.
    "nuisance_line_search_failures_detected"
  )
  matches_any <- function(status, patterns) {
    matched <- grepl(paste(patterns, collapse = "|"), status)
    matched[is.na(matched)] <- FALSE
    matched
  }
  list(
    implementation_failed = matches_any(
      diagnostics$diagnostic_status, implementation_patterns
    ),
    statistical_review_required = matches_any(
      diagnostics$diagnostic_status, review_patterns
    ),
    implementation_patterns = implementation_patterns,
    review_patterns = review_patterns
  )
}
