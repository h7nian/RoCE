#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) args[[1L]] else {
  "results/direct_tate_mc500_b5000/rhc_supported_primary"
}
result_path <- file.path(output_root, "rhc_direct_tate.rds")
if (!file.exists(result_path)) {
  stop("RHC result not found: ", result_path, call. = FALSE)
}

result <- readRDS(result_path)
fit <- result$tate_fit

scalar_or_na <- function(value, mode = c("numeric", "integer", "logical")) {
  mode <- match.arg(mode)
  if (is.null(value) || length(value) == 0L) {
    return(switch(mode, numeric = NA_real_, integer = NA_integer_, logical = NA))
  }
  switch(
    mode,
    numeric = as.numeric(value[[1L]]),
    integer = as.integer(value[[1L]]),
    logical = as.logical(value[[1L]])
  )
}

parameter_summary <- function(value) {
  line_search_failures <- suppressWarnings(
    as.integer(attr(value, "line_search_failures"))
  )
  lambda_min <- suppressWarnings(as.numeric(attr(value, "lambda_min")))
  lambda_1se <- suppressWarnings(as.numeric(attr(value, "lambda_1se")))
  lambda_before_floor <- suppressWarnings(as.numeric(
    attr(value, "lambda_selected_before_support_floor")
  ))
  lambda_selected_on_path <- suppressWarnings(as.numeric(
    attr(value, "lambda_selected_on_path")
  ))
  cv_invalid_fold_fits <- suppressWarnings(as.integer(
    attr(value, "cv_invalid_fold_fits")
  ))
  cv_invalid_lambdas <- suppressWarnings(as.integer(
    attr(value, "cv_invalid_lambdas")
  ))
  support_floor <- suppressWarnings(as.numeric(
    attr(value, "support_penalty_floor")
  ))
  support_constant_count <- suppressWarnings(as.integer(
    attr(value, "support_constant_feature_count")
  ))
  support_floor_applied <-
    isTRUE(attr(value, "support_penalty_floor_applied"))
  support_path_floor_applied <-
    isTRUE(attr(value, "support_path_floor_applied"))
  value <- as.numeric(value)
  finite <- is.finite(value)
  c(
    parameter_length = length(value),
    parameter_nonfinite = sum(!finite),
    parameter_max_abs = if (any(finite)) max(abs(value[finite])) else NA_real_,
    line_search_failures = if (length(line_search_failures) == 1L) {
      line_search_failures
    } else {
      NA_integer_
    },
    lambda_min = if (length(lambda_min) == 1L) lambda_min else NA_real_,
    lambda_1se = if (length(lambda_1se) == 1L) lambda_1se else NA_real_,
    lambda_selected_on_path = if (
      length(lambda_selected_on_path) == 1L
    ) lambda_selected_on_path else NA_real_,
    lambda_selected_before_support_floor = if (
      length(lambda_before_floor) == 1L
    ) lambda_before_floor else NA_real_,
    cv_invalid_fold_fits = if (
      length(cv_invalid_fold_fits) == 1L
    ) cv_invalid_fold_fits else NA_integer_,
    cv_invalid_lambdas = if (
      length(cv_invalid_lambdas) == 1L
    ) cv_invalid_lambdas else NA_integer_,
    support_penalty_floor = if (
      length(support_floor) == 1L
    ) support_floor else NA_real_,
    support_penalty_floor_applied = support_floor_applied,
    support_path_floor_applied = support_path_floor_applied,
    support_constant_feature_count = if (
      length(support_constant_count) == 1L
    ) support_constant_count else NA_integer_
  )
}

diagnostic_rows <- list()
append_diagnostic <- function(arm, treatment, outer_fold, source, inner_fold,
                              stage, converged, iterations, lambda,
                              update_ratio, parameters) {
  parameter <- parameter_summary(parameters)
  diagnostic_rows[[length(diagnostic_rows) + 1L]] <<- data.frame(
    arm = arm,
    treatment = treatment,
    outer_fold = outer_fold,
    source = source,
    inner_fold = inner_fold,
    stage = stage,
    converged = converged,
    iterations = iterations,
    lambda = lambda,
    update_to_threshold_ratio = update_ratio,
    parameter_length = parameter[["parameter_length"]],
    parameter_nonfinite = parameter[["parameter_nonfinite"]],
    parameter_max_abs = parameter[["parameter_max_abs"]],
    line_search_failures = parameter[["line_search_failures"]],
    lambda_min = parameter[["lambda_min"]],
    lambda_1se = parameter[["lambda_1se"]],
    lambda_selected_on_path = parameter[["lambda_selected_on_path"]],
    lambda_selected_before_support_floor =
      parameter[["lambda_selected_before_support_floor"]],
    cv_invalid_fold_fits = parameter[["cv_invalid_fold_fits"]],
    cv_invalid_lambdas = parameter[["cv_invalid_lambdas"]],
    support_penalty_floor = parameter[["support_penalty_floor"]],
    support_penalty_floor_applied =
      as.logical(parameter[["support_penalty_floor_applied"]]),
    support_path_floor_applied =
      as.logical(parameter[["support_path_floor_applied"]]),
    support_constant_feature_count =
      as.integer(parameter[["support_constant_feature_count"]]),
    stringsAsFactors = FALSE
  )
}

for (arm in intersect(c("mu1", "mu0"), names(fit$arm_results))) {
  arm_result <- fit$arm_results[[arm]]
  treatment <- if (identical(arm, "mu1")) 1L else 0L
  for (outer_fold in seq_along(arm_result$fold_results)) {
    source_results <- arm_result$fold_results[[outer_fold]]$source_results
    for (source in names(source_results)) {
      source_result <- source_results[[source]]
      diagnostic <- source_result$nuisance_fit_diagnostics
      initial_parameters <- source_result$per_k2_gamma
      n_inner <- max(
        length(initial_parameters),
        length(diagnostic$initial_density_ratio_converged),
        length(diagnostic$initial_density_ratio_iterations),
        length(diagnostic$initial_density_ratio_lambdas),
        length(diagnostic$initial_density_ratio_update_ratios)
      )
      if (n_inner > 0L) {
        for (inner_index in seq_len(n_inner)) {
          inner_name <- names(initial_parameters)[inner_index]
          if (is.null(inner_name) || is.na(inner_name) || !nzchar(inner_name)) {
            inner_name <- as.character(inner_index)
          }
          append_diagnostic(
            arm, treatment, outer_fold, source, inner_name,
            "initial_density_ratio",
            scalar_or_na(diagnostic$initial_density_ratio_converged[inner_index], "logical"),
            scalar_or_na(diagnostic$initial_density_ratio_iterations[inner_index], "integer"),
            scalar_or_na(diagnostic$initial_density_ratio_lambdas[inner_index], "numeric"),
            scalar_or_na(diagnostic$initial_density_ratio_update_ratios[inner_index], "numeric"),
            initial_parameters[[inner_index]]
          )
        }
      }
      append_diagnostic(
        arm, treatment, outer_fold, source, NA_character_,
        "calibrated_density_ratio",
        scalar_or_na(diagnostic$calibrated_density_ratio_converged, "logical"),
        scalar_or_na(diagnostic$calibrated_density_ratio_iterations, "integer"),
        scalar_or_na(diagnostic$calibrated_density_ratio_lambda, "numeric"),
        scalar_or_na(diagnostic$calibrated_density_ratio_update_ratio, "numeric"),
        source_result$gamma_s
      )
      append_diagnostic(
        arm, treatment, outer_fold, source, NA_character_,
        "calibrated_outcome",
        scalar_or_na(diagnostic$calibrated_outcome_converged, "logical"),
        scalar_or_na(diagnostic$calibrated_outcome_iterations, "integer"),
        scalar_or_na(diagnostic$calibrated_outcome_lambda, "numeric"),
        scalar_or_na(
          diagnostic$calibrated_outcome_update_ratio, "numeric"
        ),
        source_result$alpha_ts
      )
    }
  }
}

nuisance_detail <- if (length(diagnostic_rows) == 0L) data.frame() else {
  do.call(rbind, diagnostic_rows)
}
write.csv(
  nuisance_detail,
  file.path(output_root, "rhc_direct_tate_nuisance_detail.csv"),
  row.names = FALSE
)

reconstruct_target_phi <- function(fold_info) {
  if (is.null(fold_info) || length(fold_info) == 0L) return(NULL)
  max_index <- max(unlist(lapply(fold_info, `[[`, "target_idx")))
  values <- rep(NA_real_, max_index)
  counts <- integer(max_index)
  for (info in fold_info) {
    indices <- as.integer(info$target_idx)
    values[indices] <- as.numeric(info$varphi_ot) + info$fold_target_estimate
    counts[indices] <- counts[indices] + 1L
  }
  attr(values, "assignment_counts") <- counts
  values
}

pseudo_objects <- list(
  direct_tate = fit$all_phi_tau,
  target_tate_stored = fit$intermediates$target_only_phi,
  target_tate_reconstructed = reconstruct_target_phi(fit$intermediates$fold_info),
  target_mu1_stored = fit$arm_results$mu1$intermediates$target_only_phi,
  target_mu0_stored = fit$arm_results$mu0$intermediates$target_only_phi
)
pseudo_summary <- do.call(rbind, lapply(names(pseudo_objects), function(name) {
  value <- pseudo_objects[[name]]
  finite <- is.finite(value)
  counts <- attr(value, "assignment_counts")
  data.frame(
    component = name,
    class = paste(class(value), collapse = "/"),
    typeof = typeof(value),
    length = length(value),
    finite_count = sum(finite),
    nonfinite_count = sum(!finite),
    min = if (any(finite)) min(value[finite]) else NA_real_,
    max = if (any(finite)) max(value[finite]) else NA_real_,
    mean = if (any(finite)) mean(value[finite]) else NA_real_,
    assignment_min = if (length(counts)) min(counts) else NA_integer_,
    assignment_max = if (length(counts)) max(counts) else NA_integer_,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  pseudo_summary,
  file.path(output_root, "rhc_direct_tate_pseudovalue_audit.csv"),
  row.names = FALSE
)

cat("Nuisance rows:", nrow(nuisance_detail), "\n")
if (nrow(nuisance_detail) > 0L) {
  failed <- nuisance_detail[
    !is.na(nuisance_detail$converged) & !nuisance_detail$converged,
    , drop = FALSE
  ]
  cat("Nonconverged rows:", nrow(failed), "\n")
  if (nrow(failed) > 0L) print(failed, row.names = FALSE)
}
print(pseudo_summary, row.names = FALSE)
