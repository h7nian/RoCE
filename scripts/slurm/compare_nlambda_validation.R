#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: compare_nlambda_validation.R GRID_A_RAW GRID_B_RAW OUTPUT_DIR",
    call. = FALSE
  )
}

grid_directories <- vapply(
  args[1:2], normalizePath, character(1L), mustWork = TRUE
)
output_directory <- args[[3L]]
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
direct_method <- "one_round_crossfit_ate"
target_method <- "target_only_ate"
expected_methods <- c(
  "one_round_crossfit", "target_only",
  "one_round_crossfit_ate_armwise", direct_method, target_method
)
cv_diagnostic_columns <- c(
  "face_initial_dr_cv_invalid_fold_fits",
  "face_initial_dr_cv_invalid_lambdas",
  "face_calibrated_dr_cv_invalid_fold_fits",
  "face_calibrated_dr_cv_invalid_lambdas",
  "face_calibrated_outcome_cv_invalid_fold_fits",
  "face_calibrated_outcome_cv_invalid_lambdas"
)

read_grid <- function(directory) {
  files <- list.files(
    directory, pattern = "^task_[0-9]+[.]csv$", full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("no task results found in ", directory, call. = FALSE)
  }
  parts <- lapply(files, read.csv, stringsAsFactors = FALSE)
  columns <- unique(unlist(lapply(parts, names), use.names = FALSE))
  parts <- lapply(parts, function(part) {
    for (column in setdiff(columns, names(part))) {
      part[[column]] <- NA
    }
    part[columns]
  })
  result <- do.call(rbind, parts)
  rownames(result) <- NULL
  result
}

grids <- lapply(grid_directories, read_grid)
expected_cv_threads <- roce_read_positive_integer_env(
  "ROCE_EXPECT_NUISANCE_CV_THREADS", 5L, 5L
)
required_columns <- c(
  "task_id", "sim_id", "method", "estimate", "se", "bias", "coverage",
  "truth",
  "experiment", "config", "p",
  "K", "rho", "cutoff", "n_site", "n_folds", "nlambda_init",
  "target_anchor_nlambda",
  "n_bootstrap", "M_tau", "M_tau_inference", "task_elapsed_seconds",
  "package_library", "package_fingerprint",
  "workflow_fingerprint", "manifest_fingerprint",
  roce_resource_metadata_columns(),
  "target_anchor_weight", "mean_abs_source_weight",
  "max_abs_source_weight", "mean_wald_statistic", "max_wald_statistic",
  "penalized_source_fold_fraction", "face_initial_dr_nonconverged",
  "face_calibrated_dr_nonconverged",
  "face_calibrated_outcome_nonconverged",
  "face_initial_dr_line_search_failures",
  "face_calibrated_dr_line_search_failures",
  "face_max_initial_dr_abs_coefficient",
  "face_max_calibrated_dr_abs_coefficient", cv_diagnostic_columns,
  "inference_safety_clip_count"
)

provenance <- vector("list", length(grids))
workflow_fingerprints <- character(length(grids))
manifest_fingerprints <- character(length(grids))
for (index in seq_along(grids)) {
  missing_columns <- setdiff(required_columns, names(grids[[index]]))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "grid %d is missing columns: %s", index,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  core_numeric_columns <- c("estimate", "se", "bias", "truth")
  if (any(!is.finite(as.matrix(grids[[index]][core_numeric_columns]))) ||
      any(grids[[index]]$se <= 0)) {
    stop(sprintf("grid %d has non-finite core values or non-positive SEs.",
                 index), call. = FALSE)
  }
  expected_metadata <-
    grids[[index]]$experiment == "nlambda_validation" &
    grids[[index]]$config == "C3" & grids[[index]]$p == 100L &
    grids[[index]]$K == 4L & grids[[index]]$rho == 0 &
    grids[[index]]$cutoff == 2 & grids[[index]]$n_site == 1000L &
    grids[[index]]$n_folds == 5L & grids[[index]]$n_bootstrap == 5000L &
    grids[[index]]$M_tau == 5 & grids[[index]]$M_tau_inference == 5 &
    grids[[index]]$target_anchor_nlambda == 100L
  if (!all(expected_metadata)) {
    stop(sprintf("grid %d does not match the bounded p=100 C3 validation design.",
                 index), call. = FALSE)
  }
  roce_validate_result_resource_metadata(
    grids[[index]], expected_threads = expected_cv_threads
  )
  method_sets <- split(grids[[index]]$method, grids[[index]]$sim_id)
  valid_method_sets <- vapply(method_sets, function(methods) {
    length(methods) == length(expected_methods) && !anyDuplicated(methods) &&
      identical(sort(as.character(methods)), sort(expected_methods))
  }, logical(1L))
  if (!all(valid_method_sets)) {
    stop(
      sprintf(
        "grid %d has an invalid five-method schema for sim_id: %s", index,
        paste(names(valid_method_sets)[!valid_method_sets], collapse = ",")
      ),
      call. = FALSE
    )
  }
  task_seed_map <- unique(grids[[index]][c("task_id", "sim_id")])
  task_seed_map <- task_seed_map[order(task_seed_map$sim_id), , drop = FALSE]
  rownames(task_seed_map) <- NULL
  expected_seeds <- sort(unique(as.integer(grids[[index]]$sim_id)))
  if (anyNA(expected_seeds) || !identical(expected_seeds, seq_along(expected_seeds)) ||
      nrow(task_seed_map) != length(expected_seeds) ||
      !identical(as.integer(task_seed_map$task_id), expected_seeds) ||
      !identical(as.integer(task_seed_map$sim_id), expected_seeds)) {
    stop(
      sprintf("grid %d does not contain an exact task_id-to-seed 1:n mapping.",
              index),
      call. = FALSE
    )
  }
  provenance[[index]] <- roce_result_provenance(
    grids[[index]], project_library
  )
  if (!isTRUE(provenance[[index]]$passed)) {
    stop(
      sprintf(
        "grid %d failed tested-package provenance validation: %s",
        index, provenance[[index]]$detail
      ),
      call. = FALSE
    )
  }
  workflow_values <- unique(tolower(trimws(
    as.character(grids[[index]]$workflow_fingerprint)
  )))
  manifest_values <- unique(tolower(trimws(
    as.character(grids[[index]]$manifest_fingerprint)
  )))
  expected_workflow <- tolower(trimws(Sys.getenv(
    "ROCE_WORKFLOW_FINGERPRINT", ""
  )))
  if (length(workflow_values) != 1L ||
      !grepl("^[0-9a-f]{64}$", workflow_values) ||
      !identical(workflow_values, expected_workflow)) {
    stop(sprintf("grid %d has invalid or unexpected workflow provenance.", index),
         call. = FALSE)
  }
  if (length(manifest_values) != 1L ||
      !grepl("^[0-9a-f]{64}$", manifest_values)) {
    stop(sprintf("grid %d has invalid or mixed manifest provenance.", index),
         call. = FALSE)
  }
  workflow_fingerprints[[index]] <- workflow_values
  manifest_fingerprints[[index]] <- manifest_values
}
if (length(unique(vapply(
  provenance, `[[`, character(1L), "library"
))) != 1L || length(unique(vapply(
  provenance, `[[`, character(1L), "fingerprint"
))) != 1L) {
  stop(
    "paired grids do not share one tested package installation.",
    call. = FALSE
  )
}
if (length(unique(workflow_fingerprints)) != 1L) {
  stop("paired grids do not share one simulation workflow fingerprint.",
       call. = FALSE)
}
resource_metadata <- lapply(grids, function(grid) {
  unique(grid[roce_resource_metadata_columns()])
})
if (!identical(resource_metadata[[1L]], resource_metadata[[2L]])) {
  stop("paired grids do not share one Slurm resource topology.", call. = FALSE)
}

grid_sizes <- vapply(grids, function(grid) {
  values <- unique(grid$nlambda_init)
  if (length(values) != 1L || !is.finite(values) || values < 2L) {
    stop("each grid must contain one valid nlambda_init.", call. = FALSE)
  }
  as.integer(values)
}, integer(1L))
if (grid_sizes[[1L]] == grid_sizes[[2L]]) {
  stop("the paired grids must use different nlambda_init values.", call. = FALSE)
}

metadata_columns <- c(
  "sim_id", "config", "p", "K", "rho", "cutoff", "n_site", "n_folds",
  "n_bootstrap", "M_tau", "M_tau_inference"
)
metadata <- lapply(grids, function(grid) {
  unique(grid[grid$method == direct_method, metadata_columns, drop = FALSE])
})
metadata <- lapply(metadata, function(values) {
  values[order(values$sim_id), , drop = FALSE]
})
if (!identical(metadata[[1L]], metadata[[2L]])) {
  stop("paired grids do not have identical seed and DGP metadata.", call. = FALSE)
}

select_method <- function(grid, method) {
  rows <- grid[grid$method == method, , drop = FALSE]
  rows[order(rows$sim_id), , drop = FALSE]
}
direct <- lapply(grids, select_method, method = direct_method)
target <- lapply(grids, select_method, method = target_method)
if (!identical(direct[[1L]]$truth, direct[[2L]]$truth)) {
  stop("paired grids use different TATE truths.", call. = FALSE)
}
if (!isTRUE(all.equal(
  target[[1L]][c("estimate", "se", "truth")],
  target[[2L]][c("estimate", "se", "truth")],
  tolerance = 0, check.attributes = FALSE
))) {
  stop(
    "the fold-aligned target anchor changed across nuisance-grid sizes.",
    call. = FALSE
  )
}

diagnostic_columns <- c(
  "target_anchor_weight", "mean_abs_source_weight",
  "max_abs_source_weight", "mean_wald_statistic", "max_wald_statistic",
  "penalized_source_fold_fraction", "face_initial_dr_nonconverged",
  "face_calibrated_dr_nonconverged",
  "face_calibrated_outcome_nonconverged",
  "face_initial_dr_line_search_failures",
  "face_calibrated_dr_line_search_failures",
  "face_max_initial_dr_abs_coefficient",
  "face_max_calibrated_dr_abs_coefficient", cv_diagnostic_columns,
  "inference_safety_clip_count"
)
for (index in seq_along(direct)) {
  if (any(!is.finite(as.matrix(direct[[index]][diagnostic_columns])))) {
    stop(sprintf("grid %d has non-finite TATE diagnostics.", index),
         call. = FALSE)
  }
}

paired <- data.frame(
  sim_id = direct[[1L]]$sim_id,
  nlambda_a = grid_sizes[[1L]],
  nlambda_b = grid_sizes[[2L]],
  estimate_a = direct[[1L]]$estimate,
  estimate_b = direct[[2L]]$estimate,
  estimate_difference = direct[[1L]]$estimate - direct[[2L]]$estimate,
  se_a = direct[[1L]]$se,
  se_b = direct[[2L]]$se,
  se_ratio_a_over_b = direct[[1L]]$se / direct[[2L]]$se,
  difference_in_b_se =
    (direct[[1L]]$estimate - direct[[2L]]$estimate) / direct[[2L]]$se,
  runtime_a_seconds = direct[[1L]]$task_elapsed_seconds,
  runtime_b_seconds = direct[[2L]]$task_elapsed_seconds,
  runtime_ratio_a_over_b =
    direct[[1L]]$task_elapsed_seconds / direct[[2L]]$task_elapsed_seconds,
  stringsAsFactors = FALSE
)
for (column in diagnostic_columns) {
  paired[[paste0(column, "_a")]] <- direct[[1L]][[column]]
  paired[[paste0(column, "_b")]] <- direct[[2L]][[column]]
  paired[[paste0(column, "_difference")]] <-
    direct[[1L]][[column]] - direct[[2L]][[column]]
}

summary <- data.frame(
  package_library = provenance[[1L]]$library,
  package_fingerprint = provenance[[1L]]$fingerprint,
  workflow_fingerprint = workflow_fingerprints[[1L]],
  manifest_fingerprint_a = manifest_fingerprints[[1L]],
  manifest_fingerprint_b = manifest_fingerprints[[2L]],
  n_paired_seeds = nrow(paired),
  nlambda_a = grid_sizes[[1L]],
  nlambda_b = grid_sizes[[2L]],
  mean_estimate_difference = mean(paired$estimate_difference),
  max_abs_estimate_difference = max(abs(paired$estimate_difference)),
  rms_estimate_difference = sqrt(mean(paired$estimate_difference^2)),
  mean_abs_difference_in_b_se = mean(abs(paired$difference_in_b_se)),
  max_abs_difference_in_b_se = max(abs(paired$difference_in_b_se)),
  mean_se_ratio_a_over_b = mean(paired$se_ratio_a_over_b),
  mean_runtime_ratio_a_over_b = mean(paired$runtime_ratio_a_over_b),
  median_runtime_ratio_a_over_b = stats::median(
    paired$runtime_ratio_a_over_b
  ),
  any_nuisance_nonconvergence = any(
    paired$face_initial_dr_nonconverged_a > 0 |
      paired$face_initial_dr_nonconverged_b > 0 |
      paired$face_calibrated_dr_nonconverged_a > 0 |
      paired$face_calibrated_dr_nonconverged_b > 0 |
      paired$face_calibrated_outcome_nonconverged_a > 0 |
      paired$face_calibrated_outcome_nonconverged_b > 0
  ),
  any_density_ratio_line_search_failure = any(
    paired$face_initial_dr_line_search_failures_a > 0 |
      paired$face_initial_dr_line_search_failures_b > 0 |
      paired$face_calibrated_dr_line_search_failures_a > 0 |
      paired$face_calibrated_dr_line_search_failures_b > 0
  ),
  any_density_ratio_boundary = any(
    paired$face_max_initial_dr_abs_coefficient_a >= 90 |
      paired$face_max_initial_dr_abs_coefficient_b >= 90 |
      paired$face_max_calibrated_dr_abs_coefficient_a >= 90 |
      paired$face_max_calibrated_dr_abs_coefficient_b >= 90
  ),
  any_cv_candidate_excluded = any(vapply(
    c(
      paste0(cv_diagnostic_columns, "_a"),
      paste0(cv_diagnostic_columns, "_b")
    ),
    function(column) any(paired[[column]] > 0),
    logical(1L)
  )),
  invalid_cv_fold_fits_a = sum(unlist(paired[
    paste0(cv_diagnostic_columns[grepl("fold_fits$", cv_diagnostic_columns)],
           "_a")
  ])),
  invalid_cv_fold_fits_b = sum(unlist(paired[
    paste0(cv_diagnostic_columns[grepl("fold_fits$", cv_diagnostic_columns)],
           "_b")
  ])),
  any_inference_safety_clip = any(
    paired$inference_safety_clip_count_a > 0 |
      paired$inference_safety_clip_count_b > 0
  ),
  stringsAsFactors = FALSE
)
summary$accelerated_grid_eligible <-
  summary$n_paired_seeds >= 5L &&
  summary$nlambda_a < summary$nlambda_b &&
  summary$mean_abs_difference_in_b_se <= 0.10 &&
  summary$max_abs_difference_in_b_se <= 0.25 &&
  summary$mean_se_ratio_a_over_b >= 0.98 &&
  summary$mean_se_ratio_a_over_b <= 1.02 &&
  summary$median_runtime_ratio_a_over_b <= 0.80 &&
  !summary$any_nuisance_nonconvergence &&
  !summary$any_density_ratio_line_search_failure &&
  !summary$any_density_ratio_boundary &&
  !summary$any_inference_safety_clip

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write_atomic_csv <- function(value, path) {
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "_"),
    tmpdir = dirname(path), fileext = ".tmp"
  )
  write.csv(value, temporary_path, row.names = FALSE)
  if (!file.rename(temporary_path, path)) {
    unlink(temporary_path)
    stop("failed to atomically write ", path, call. = FALSE)
  }
}
write_atomic_csv(paired, file.path(output_directory, "paired_seed_results.csv"))
write_atomic_csv(summary, file.path(output_directory, "paired_summary.csv"))

if (summary$any_nuisance_nonconvergence ||
    summary$any_density_ratio_line_search_failure ||
    summary$any_density_ratio_boundary ||
    summary$any_inference_safety_clip) {
  stop(
    "paired nlambda validation contains an implementation-QC failure; ",
    "inspect paired_seed_results.csv before proceeding.",
    call. = FALSE
  )
}

sentinel_path <- file.path(
  output_directory, "implementation_audit_passed.txt"
)
sentinel_temporary <- tempfile(
  pattern = ".implementation_audit_passed_",
  tmpdir = output_directory, fileext = ".tmp"
)
writeLines(
  c(
    paste0("n_paired_seeds=", summary$n_paired_seeds),
    paste0("package_fingerprint=", summary$package_fingerprint),
    paste0("workflow_fingerprint=", summary$workflow_fingerprint),
    "implementation_qc=passed"
  ),
  sentinel_temporary
)
if (!file.rename(sentinel_temporary, sentinel_path)) {
  unlink(sentinel_temporary)
  stop("failed to atomically write paired-audit sentinel.", call. = FALSE)
}

print(summary, row.names = FALSE, digits = 6)
message(
  "[done] paired nuisance-grid comparison written to ", output_directory,
  ". Accelerated-grid eligibility requires five paired seeds and the ",
  "pre-specified numerical/SE/runtime thresholds. This bounded study is not ",
  "a coverage experiment."
)
