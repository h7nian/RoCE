#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    paste0(
      "usage: compare_nlambda_smokes.R GRID_A_RESULT.csv ",
      "GRID_B_RESULT.csv [OUTPUT.csv]"
    ),
    call. = FALSE
  )
}

paths <- vapply(args[1:2], normalizePath, character(1L), mustWork = TRUE)
output_path <- if (length(args) == 3L) args[[3L]] else NA_character_
results <- lapply(paths, read.csv, stringsAsFactors = FALSE)

required <- c(
  "sim_id", "method", "estimate", "se", "truth", "config", "p", "K",
  "rho", "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
  "M_tau", "M_tau_inference", "task_elapsed_seconds"
)
for (i in seq_along(results)) {
  missing <- setdiff(required, names(results[[i]]))
  if (length(missing) > 0L) {
    stop(
      sprintf("result %d is missing columns: %s", i, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (anyDuplicated(results[[i]]$method)) {
    stop(sprintf("result %d contains duplicated methods.", i), call. = FALSE)
  }
}

metadata_fields <- c(
  "sim_id", "config", "p", "K", "rho", "cutoff", "n_site", "n_folds",
  "n_bootstrap", "M_tau", "M_tau_inference"
)
metadata <- lapply(results, function(result) {
  unique(result[metadata_fields])
})
if (nrow(metadata[[1L]]) != 1L || nrow(metadata[[2L]]) != 1L ||
    !identical(metadata[[1L]], metadata[[2L]])) {
  stop("the two smoke results do not have identical non-grid metadata.", call. = FALSE)
}

grid_sizes <- vapply(results, function(result) {
  values <- unique(result$nlambda_init)
  if (length(values) != 1L || !is.finite(values) || values < 2L) {
    stop("each smoke result must contain exactly one valid nlambda_init.", call. = FALSE)
  }
  as.integer(values)
}, integer(1L))
if (grid_sizes[[1L]] == grid_sizes[[2L]]) {
  stop("the two smoke results use the same nlambda_init.", call. = FALSE)
}

methods <- sort(results[[1L]]$method)
if (!identical(methods, sort(results[[2L]]$method))) {
  stop("the two smoke results do not contain the same method set.", call. = FALSE)
}

ordered <- lapply(results, function(result) {
  result[match(methods, result$method), , drop = FALSE]
})
if (!isTRUE(all.equal(
  ordered[[1L]]$truth, ordered[[2L]]$truth,
  tolerance = 0, check.attributes = FALSE
))) {
  stop("the two smoke results do not use identical estimand truths.", call. = FALSE)
}

comparison <- data.frame(
  method = methods,
  nlambda_a = grid_sizes[[1L]],
  nlambda_b = grid_sizes[[2L]],
  estimate_a = ordered[[1L]]$estimate,
  estimate_b = ordered[[2L]]$estimate,
  estimate_difference = ordered[[1L]]$estimate - ordered[[2L]]$estimate,
  se_a = ordered[[1L]]$se,
  se_b = ordered[[2L]]$se,
  se_ratio_a_over_b = ordered[[1L]]$se / ordered[[2L]]$se,
  difference_in_b_se = (
    ordered[[1L]]$estimate - ordered[[2L]]$estimate
  ) / ordered[[2L]]$se,
  stringsAsFactors = FALSE
)

diagnostic_fields <- c(
  "target_anchor_weight", "mean_abs_source_weight", "max_abs_source_weight",
  "mean_wald_statistic", "max_wald_statistic",
  "penalized_source_fold_fraction", "face_initial_dr_nonconverged",
  "face_calibrated_dr_nonconverged", "face_calibrated_outcome_nonconverged",
  "face_initial_dr_line_search_failures",
  "face_calibrated_dr_line_search_failures",
  "face_initial_dr_support_floor_applied",
  "face_calibrated_dr_support_floor_applied", "max_weight_psd_ridge",
  "inference_safety_clip_count"
)
missing_diagnostics <- setdiff(
  diagnostic_fields, intersect(names(results[[1L]]), names(results[[2L]]))
)
if (length(missing_diagnostics) > 0L) {
  stop(
    paste0(
      "paired TATE diagnostics are incomplete: ",
      paste(missing_diagnostics, collapse = ", ")
    ),
    call. = FALSE
  )
}

direct_method <- "one_round_crossfit_ate"
direct <- lapply(results, function(result) {
  row <- result[result$method == direct_method, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("each result must contain exactly one TATE row.", call. = FALSE)
  }
  row
})
diagnostics <- data.frame(
  diagnostic = diagnostic_fields,
  value_a = unlist(direct[[1L]][diagnostic_fields], use.names = FALSE),
  value_b = unlist(direct[[2L]][diagnostic_fields], use.names = FALSE),
  stringsAsFactors = FALSE
)
diagnostics$difference_a_minus_b <- diagnostics$value_a - diagnostics$value_b

elapsed <- vapply(results, function(result) {
  values <- unique(result$task_elapsed_seconds)
  if (length(values) != 1L || !is.finite(values) || values <= 0) {
    stop("each result must contain exactly one positive task runtime.", call. = FALSE)
  }
  values
}, numeric(1L))
runtime <- data.frame(
  nlambda_a = grid_sizes[[1L]],
  nlambda_b = grid_sizes[[2L]],
  seconds_a = elapsed[[1L]],
  seconds_b = elapsed[[2L]],
  runtime_ratio_a_over_b = elapsed[[1L]] / elapsed[[2L]],
  stringsAsFactors = FALSE
)

cat("Paired nuisance-grid comparison (same seed and DGP)\n")
print(comparison, row.names = FALSE, digits = 6)
cat("\nTATE diagnostics\n")
print(diagnostics, row.names = FALSE, digits = 6)
cat("\nRuntime\n")
print(runtime, row.names = FALSE, digits = 6)
cat(
  "\nInterpretation: this one-seed comparison is an implementation and runtime ",
  "diagnostic; it is not evidence of Monte Carlo equivalence.\n",
  sep = ""
)

if (!is.na(output_path)) {
  output_dir <- dirname(output_path)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output <- rbind(
    transform(
      comparison,
      record_type = "method",
      name = method,
      value_a = estimate_a,
      value_b = estimate_b,
      difference = estimate_difference,
      ratio = se_ratio_a_over_b
    )[c(
      "record_type", "name", "nlambda_a", "nlambda_b", "value_a",
      "value_b", "difference", "ratio"
    )],
    transform(
      diagnostics,
      record_type = "diagnostic",
      name = diagnostic,
      nlambda_a = grid_sizes[[1L]],
      nlambda_b = grid_sizes[[2L]],
      difference = difference_a_minus_b,
      ratio = ifelse(value_b == 0, NA_real_, value_a / value_b)
    )[c(
      "record_type", "name", "nlambda_a", "nlambda_b", "value_a",
      "value_b", "difference", "ratio"
    )],
    data.frame(
      record_type = "runtime",
      name = "task_elapsed_seconds",
      nlambda_a = grid_sizes[[1L]],
      nlambda_b = grid_sizes[[2L]],
      value_a = elapsed[[1L]],
      value_b = elapsed[[2L]],
      difference = elapsed[[1L]] - elapsed[[2L]],
      ratio = elapsed[[1L]] / elapsed[[2L]],
      stringsAsFactors = FALSE
    )
  )
  temporary_path <- tempfile(
    pattern = ".nlambda_comparison_", tmpdir = output_dir, fileext = ".csv"
  )
  write.csv(output, temporary_path, row.names = FALSE)
  if (!file.rename(temporary_path, output_path)) {
    unlink(temporary_path)
    stop("failed to atomically write ", output_path, call. = FALSE)
  }
  message("[done] wrote ", output_path)
}
