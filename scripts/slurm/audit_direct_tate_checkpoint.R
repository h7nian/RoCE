#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6L || length(args) > 7L) {
  stop(paste(
    "usage: audit_direct_tate_checkpoint.R ROOT MANIFEST CONFIG K RHO",
    "N_REPLICATIONS [CUTOFF]"
  ), call. = FALSE)
}

root <- args[[1L]]
manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
config <- args[[3L]]
source_count <- as.integer(args[[4L]])
rho <- as.numeric(args[[5L]])
expected_replications <- as.integer(args[[6L]])
cutoff <- if (length(args) == 7L) {
  as.numeric(args[[7L]])
} else {
  suppressWarnings(as.numeric(Sys.getenv("ROCE_PRIMARY_CUTOFF", "1")))
}

if (!config %in% c("C1", "C2", "C3")) {
  stop("CONFIG must be C1, C2, or C3.", call. = FALSE)
}
if (is.na(source_count) || !source_count %in% c(2L, 4L, 8L)) {
  stop("K must be 2, 4, or 8.", call. = FALSE)
}
if (!is.finite(rho) || !rho %in% c(0, 0.5, 1, 1.5, 2, 2.5)) {
  stop("RHO is not in the production grid.", call. = FALSE)
}
if (is.na(expected_replications) || expected_replications < 1L ||
    expected_replications > 500L) {
  stop("N_REPLICATIONS must be an integer in 1:500.", call. = FALSE)
}
if (!is.finite(cutoff) || cutoff <= 0) {
  stop("CUTOFF must be positive.", call. = FALSE)
}

default_project_library <- file.path(
  normalizePath(".", mustWork = TRUE),
  "results", "direct_tate_mc500_b5000", "Rlib_current"
)
project_library <- Sys.getenv(
  "ROCE_PROJECT_LIB",
  if (dir.exists(default_project_library)) default_project_library else ""
)
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
package_primary_cutoff <- 1 / get(
  "AGG_WALD_LAMBDA", envir = asNamespace("RoCE")
)
if (!isTRUE(all.equal(cutoff, package_primary_cutoff, tolerance = 1e-12))) {
  stop(
    "requested checkpoint cutoff does not match the tested package default.",
    call. = FALSE
  )
}
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
source(file.path("scripts", "slurm", "simulation_qc_policy.R"))

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
required_manifest_columns <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
  "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap", "M_tau",
  "M_tau_inference", "methods"
)
missing_manifest_columns <- setdiff(required_manifest_columns, names(manifest))
if (length(missing_manifest_columns) > 0L) {
  stop(
    "manifest is missing columns: ",
    paste(missing_manifest_columns, collapse = ", "),
    call. = FALSE
  )
}

setting <- manifest[
  manifest$experiment == "negative_transfer" &
    manifest$config == config & manifest$p == 100L &
    manifest$K == source_count & manifest$rho == rho &
    manifest$cutoff == cutoff,
  , drop = FALSE
]
setting <- setting[order(setting$sim_id), , drop = FALSE]
expected_ids <- seq_len(expected_replications)
setting <- setting[setting$sim_id %in% expected_ids, , drop = FALSE]
if (nrow(setting) != expected_replications ||
    !identical(as.integer(setting$sim_id), expected_ids)) {
  stop(
    "manifest does not contain exactly sim_id 1:", expected_replications,
    " for the requested setting.",
    call. = FALSE
  )
}

raw_directory <- file.path(root, "raw")
task_files <- file.path(
  raw_directory,
  sprintf("task_%06d.csv", as.integer(setting$task_id))
)
missing_files <- task_files[!file.exists(task_files)]
if (length(missing_files) > 0L) {
  stop(
    sprintf(
      "checkpoint is incomplete: %d/%d expected task files are missing; first=%s",
      length(missing_files), expected_replications, missing_files[[1L]]
    ),
    call. = FALSE
  )
}

raw <- RoCE:::.read_simulation_result_files(task_files)
required_result_columns <- c(
  "task_id", "sim_id", "experiment", "method", "estimand_scope",
  "dgp_type", "outcome_family", "heterogeneity_type", "estimand_type",
  "estimate", "se", "bias", "coverage", "truth", "ci_lower", "ci_upper",
  "p", "K", "config", "rho", "cutoff", "n_site", "n_folds",
  "aggregation_lambda", "aggregation_cutoff",
  "primary_cutoff",
  "nlambda_init", "nuisance_lambda_rule", "n_bootstrap", "M_tau",
  "M_tau_inference",
  "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
  "n_site_arm_outcome_cells_below_8",
  "dr_weight_n", "dr_weight_n_clipped",
  "dr_weight_fraction_clipped", "dr_weight_max_site_fraction_clipped",
  "dr_weight_min_before_clipping", "dr_weight_max_before_clipping",
  "comparison_variance_method", "comparison_se_analytic",
  "comparison_se_bootstrap", "comparison_n_bootstrap",
  RoCE:::.direct_tate_aggregation_diagnostic_columns(),
  "target_anchor_nlambda", "package_library", "package_fingerprint",
  "workflow_fingerprint", "manifest_fingerprint",
  roce_resource_metadata_columns(), roce_scheduler_provenance_columns()
)
missing_result_columns <- setdiff(required_result_columns, names(raw))
if (length(missing_result_columns) > 0L) {
  stop(
    "task results are missing columns: ",
    paste(missing_result_columns, collapse = ", "),
    call. = FALSE
  )
}
expected_methods <- c(
  "one_round_crossfit", "target_only", "sample_size", "inverse_variance",
  "federated_dr", "pooled_dr", "one_round_crossfit_ate_armwise",
  "one_round_crossfit_ate", "target_only_ate", "sample_size_ate",
  "inverse_variance_ate", "federated_dr_ate", "pooled_dr_ate"
)
# The quadratic-bias row is requested by default; results written before the
# column existed fail the method-set check below by design.
quadratic_bias_rule_disabled <-
  "quadratic_bias_rule_requested" %in% names(raw) &&
  all(!is.na(raw$quadratic_bias_rule_requested)) &&
  all(!as.logical(raw$quadratic_bias_rule_requested))
if (!quadratic_bias_rule_disabled) {
  expected_methods <- c(
    expected_methods, "one_round_crossfit_ate_quadratic_bias"
  )
}
hard_threshold_requested <-
  "hard_threshold_diagnostic_requested" %in% names(raw) &&
  all(!is.na(raw$hard_threshold_diagnostic_requested)) &&
  all(as.logical(raw$hard_threshold_diagnostic_requested))
if (hard_threshold_requested) {
  expected_methods <- c(
    expected_methods, "one_round_crossfit_ate_hard_threshold"
  )
}

rows_by_task <- split(raw, raw$task_id)
valid_method_sets <- vapply(rows_by_task, function(task_rows) {
  nrow(task_rows) == length(expected_methods) &&
    !anyDuplicated(task_rows$method) &&
    identical(sort(as.character(task_rows$method)), sort(expected_methods))
}, logical(1L))
if (length(rows_by_task) != expected_replications || !all(valid_method_sets)) {
  bad_tasks <- names(valid_method_sets)[!valid_method_sets]
  stop(
    "one or more tasks do not contain the expected result schema: ",
    paste(head(bad_tasks, 5L), collapse = ","),
    call. = FALSE
  )
}

observed_task_map <- unique(raw[c("task_id", "sim_id")])
observed_task_map <- observed_task_map[order(observed_task_map$sim_id), ]
expected_task_map <- setting[c("task_id", "sim_id")]
expected_nlambda <- setting$nlambda_init[
  match(as.integer(raw$task_id), as.integer(setting$task_id))
]
rownames(observed_task_map) <- NULL
rownames(expected_task_map) <- NULL
provenance <- roce_result_provenance(raw, project_library)
expected_cv_threads <- roce_read_positive_integer_env(
  "ROCE_EXPECT_NUISANCE_CV_THREADS",
  roce_expected_production_cv_threads(source_count),
  5L
)
resource_metadata_ok <- isTRUE(tryCatch(
  roce_validate_result_resource_metadata(
    raw, expected_threads = expected_cv_threads
  ),
  error = function(error) FALSE
))
scheduler_provenance_ok <- isTRUE(tryCatch(
  roce_validate_scheduler_provenance(raw, require_slurm = TRUE),
  error = function(error) FALSE
))
RoCE:::.validate_face_production_scientific_metadata(raw)
expected_workflow_fingerprint <- tolower(trimws(Sys.getenv(
  "ROCE_WORKFLOW_FINGERPRINT", ""
)))
expected_manifest_fingerprint <- tolower(trimws(Sys.getenv(
  "ROCE_MANIFEST_FINGERPRINT", ""
)))
metadata_ok <- all(raw$experiment == "negative_transfer") &&
  all(raw$p == 100L) && all(raw$config == config) &&
  all(raw$K == source_count) && all(raw$rho == rho) &&
  all(raw$cutoff == cutoff) &&
  all(abs(raw$primary_cutoff - cutoff) <= 1e-12) &&
  all(abs(raw$aggregation_cutoff - cutoff) <= 1e-12) &&
  all(abs(raw$aggregation_lambda - 1 / cutoff) <= 1e-12) &&
  all(raw$n_bootstrap == 5000L) &&
  all(raw$nuisance_lambda_rule == "min") &&
  !anyNA(expected_nlambda) &&
  all(as.integer(raw$nlambda_init) == as.integer(expected_nlambda)) &&
  all(as.integer(raw$target_anchor_nlambda) == 100L) &&
  all(raw$n_site == 1000L) && all(raw$n_folds == 5L) &&
  all(raw$M_tau == 5) && all(raw$M_tau_inference == 5) &&
  all(tolower(raw$workflow_fingerprint) == expected_workflow_fingerprint) &&
  all(tolower(raw$manifest_fingerprint) == expected_manifest_fingerprint) &&
  grepl("^[0-9a-f]{64}$", expected_workflow_fingerprint) &&
  grepl("^[0-9a-f]{64}$", expected_manifest_fingerprint) &&
  setequal(as.integer(raw$task_id), as.integer(setting$task_id)) &&
  setequal(as.integer(raw$sim_id), expected_ids) &&
  identical(observed_task_map, expected_task_map) && provenance$passed &&
  resource_metadata_ok && scheduler_provenance_ok
if (!isTRUE(metadata_ok)) {
  stop(
    "replicate metadata/provenance do not match the requested setting: ",
    provenance$detail,
    call. = FALSE
  )
}

diagnostics <- RoCE:::diagnose_simulation_results(
  raw, expected_replications = expected_replications
)
diagnostics <- RoCE:::add_simulation_rmse_comparisons(diagnostics)
diagnostics <- roce_augment_simulation_qc(
  diagnostics,
  parameter_bound = get("PARAM_MAX", envir = asNamespace("RoCE"))
)
qc_classification <- roce_classify_simulation_qc(diagnostics)
implementation_failed <- qc_classification$implementation_failed
review_required <- qc_classification$statistical_review_required
diagnostics$implementation_failed <- implementation_failed
diagnostics$statistical_review_required <- review_required

rho_label <- gsub("[.]", "p", format(rho, scientific = FALSE, trim = TRUE))
checkpoint_name <- sprintf(
  "%s_K%d_rho%s_n%03d", config, source_count, rho_label,
  expected_replications
)
job_label <- Sys.getenv("SLURM_JOB_ID", "local")
run_label <- sprintf(
  "%s_job%s", format(Sys.time(), "%Y%m%d_%H%M%S"),
  gsub("[^A-Za-z0-9_-]", "_", job_label)
)
checkpoint_directory <- file.path(
  root, "checkpoints", checkpoint_name, run_label
)
if (dir.exists(checkpoint_directory)) {
  stop(
    "refusing to overwrite an existing checkpoint audit: ",
    checkpoint_directory,
    call. = FALSE
  )
}
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
write.csv(
  raw,
  file.path(checkpoint_directory, "replicate_results.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics,
  file.path(checkpoint_directory, "setting_method_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics[review_required, , drop = FALSE],
  file.path(checkpoint_directory, "statistical_review_queue.csv"),
  row.names = FALSE
)

tate <- diagnostics[diagnostics$estimand_scope == "tate", , drop = FALSE]
tate <- tate[order(tate$rmse_rank, tate$method), , drop = FALSE]
write.csv(
  tate,
  file.path(checkpoint_directory, "tate_coverage_rmse_audit.csv"),
  row.names = FALSE
)

paired_benchmarks <- RoCE:::.tate_benchmark_methods()
paired_mse <- RoCE:::summarize_tate_paired_mse_benchmarks(raw)
write.csv(
  paired_mse,
  file.path(checkpoint_directory, "direct_tate_paired_mse_audit.csv"),
  row.names = FALSE
)

hard_threshold_paired <- NULL
if (hard_threshold_requested) {
  hard_vs_target <- RoCE:::summarize_paired_mse_differences(
    raw,
    method = "one_round_crossfit_ate_hard_threshold",
    benchmark_method = "target_only_ate"
  )
  hard_vs_soft <- RoCE:::summarize_paired_mse_differences(
    raw,
    method = "one_round_crossfit_ate_hard_threshold",
    benchmark_method = "one_round_crossfit_ate"
  )
  hard_threshold_paired <- rbind(hard_vs_target, hard_vs_soft)
  write.csv(
    hard_threshold_paired,
    file.path(checkpoint_directory, "hard_threshold_paired_mse_audit.csv"),
    row.names = FALSE
  )
}

direct <- tate[tate$method == "one_round_crossfit_ate", , drop = FALSE]
hard_threshold <- tate[
  tate$method == "one_round_crossfit_ate_hard_threshold", , drop = FALSE
]
target <- tate[tate$method == "target_only_ate", , drop = FALSE]
paired_vs_target <- paired_mse[
  paired_mse$benchmark_method == "target_only_ate", , drop = FALSE
]
nonadaptive_benchmarks <- setdiff(paired_benchmarks, "target_only_ate")
comparison_tate <- tate[
  tate$method %in% nonadaptive_benchmarks, , drop = FALSE
]
best_benchmark <- if (nrow(comparison_tate) > 0L) {
  comparison_tate$method[[which.min(comparison_tate$rmse)]]
} else {
  NA_character_
}
paired_vs_best <- paired_mse[
  paired_mse$benchmark_method == best_benchmark, , drop = FALSE
]
report <- c(
  sprintf("Checkpoint: %s", checkpoint_name),
  sprintf("Package provenance: %s", provenance$detail),
  sprintf("Workflow fingerprint: %s", expected_workflow_fingerprint),
  sprintf("Manifest fingerprint: %s", expected_manifest_fingerprint),
  sprintf("Task files audited: %d", expected_replications),
  sprintf("Result rows audited: %d", nrow(raw)),
  sprintf("Implementation failures: %d", sum(implementation_failed)),
  sprintf("Statistical review rows: %d", sum(review_required)),
  sprintf("Hard-threshold diagnostic requested: %s", hard_threshold_requested),
  sprintf(
    "Density-ratio boundary rows: %d",
    sum(diagnostics$nuisance_density_ratio_boundary_detected)
  ),
  sprintf(
    "Outcome-model boundary combinations: %d",
    sum(diagnostics$nuisance_outcome_boundary_detected)
  ),
  sprintf(
    "Sparse binary-cell setting-method rows: %d",
    sum(diagnostics$sparse_binary_cell_replications > 0, na.rm = TRUE)
  ),
  sprintf(
    "Density-ratio-clipped setting-method rows: %d",
    sum(diagnostics$dr_weight_clipped_replications > 0, na.rm = TRUE)
  ),
  if (nrow(direct) == 1L) sprintf(
    paste0(
      "TATE: bias=%.6g (MCSE %.6g), RMSE=%.6g, coverage=%.4f ",
      "(MCSE %.4f), meanSE/empSD=%.4f, target-relative RMSE=%.4f"
    ),
    direct$bias, direct$bias_mcse, direct$rmse, direct$coverage,
    direct$coverage_mcse, direct$se_to_empirical_sd,
    direct$rmse_relative_to_target
  ) else "TATE row missing.",
  if (nrow(hard_threshold) == 1L) sprintf(
    paste0(
      "Hard-threshold TATE: bias=%.6g (MCSE %.6g), RMSE=%.6g, ",
      "coverage=%.4f, meanSE/empSD=%.4f"
    ),
    hard_threshold$bias, hard_threshold$bias_mcse,
    hard_threshold$rmse, hard_threshold$coverage,
    hard_threshold$se_to_empirical_sd
  ) else if (hard_threshold_requested) {
    "Hard-threshold TATE row missing."
  } else {
    "Hard-threshold TATE not requested."
  },
  if (nrow(target) == 1L) sprintf(
    "Target-only TATE: bias=%.6g, RMSE=%.6g, coverage=%.4f, meanSE/empSD=%.4f",
    target$bias, target$rmse, target$coverage, target$se_to_empirical_sd
  ) else "Target-only TATE row missing.",
  if (nrow(paired_vs_target) == 1L) sprintf(
    paste0(
      "Paired MSE difference, direct minus target: %.6g ",
      "(MCSE %.6g, z=%.3f); ",
      "direct lower on %.3f of replicates"
    ),
    paired_vs_target$mean_squared_error_difference,
    paired_vs_target$squared_error_difference_mcse,
    paired_vs_target$squared_error_difference_z,
    paired_vs_target$method_lower_squared_error_fraction
  ) else "Paired direct-vs-target MSE row missing.",
  if (nrow(paired_vs_best) == 1L) sprintf(
    paste0(
      "Paired MSE difference, direct minus observed lowest-RMSE ",
      "nonadaptive benchmark %s: %.6g ",
      "(MCSE %.6g, z=%.3f); ",
      "direct lower on %.3f of replicates"
    ),
    paired_vs_best$benchmark_method,
    paired_vs_best$mean_squared_error_difference,
    paired_vs_best$squared_error_difference_mcse,
    paired_vs_best$squared_error_difference_z,
    paired_vs_best$method_lower_squared_error_fraction
  ) else paste0(
    "Paired direct-vs-observed-lowest-RMSE-nonadaptive-benchmark ",
    "row missing."
  )
)
writeLines(report, file.path(checkpoint_directory, "checkpoint_report.txt"))
message(paste(report, collapse = "\n"))

if (any(implementation_failed)) {
  failed_methods <- diagnostics$method[implementation_failed]
  stop(
    "checkpoint implementation audit failed for: ",
    paste(failed_methods, collapse = ", "),
    call. = FALSE
  )
}

# Write the machine-readable gate only after every implementation check has
# passed.  Downstream bounded-batch submitters can therefore distinguish a
# completed audit from a directory that merely contains partial reports.  The
# timestamped checkpoint directory is immutable, so the gate also records the
# exact package, workflow, manifest, and cumulative replication count that were
# reviewed.
gate <- data.frame(
  checkpoint = checkpoint_name,
  config = config,
  K = source_count,
  rho = rho,
  cutoff = cutoff,
  expected_replications = expected_replications,
  package_library = provenance$library,
  package_fingerprint = provenance$fingerprint,
  workflow_fingerprint = expected_workflow_fingerprint,
  manifest_fingerprint = expected_manifest_fingerprint,
  statistical_review_rows = sum(review_required),
  sparse_binary_cell_rows = sum(
    diagnostics$sparse_binary_cell_replications > 0, na.rm = TRUE
  ),
  density_ratio_clipping_rows = sum(
    diagnostics$dr_weight_clipped_replications > 0, na.rm = TRUE
  ),
  audit_job_id = job_label,
  audit_completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  stringsAsFactors = FALSE
)
gate_path <- file.path(
  checkpoint_directory, "implementation_audit_passed.csv"
)
gate_temporary_path <- tempfile(
  pattern = ".implementation_audit_passed_",
  tmpdir = checkpoint_directory,
  fileext = ".csv"
)
write.csv(gate, gate_temporary_path, row.names = FALSE)
if (!file.rename(gate_temporary_path, gate_path)) {
  unlink(gate_temporary_path)
  stop("failed to atomically write checkpoint pass gate.", call. = FALSE)
}
message(
  "[pass] implementation checks passed. Statistical review signals are ",
  "reported separately and must be inspected before the next batch."
)
