#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000"
}
expected_replications <- if (length(args) >= 2L) {
  as.integer(args[[2L]])
} else {
  500L
}
manifest_path <- if (length(args) >= 3L) args[[3L]] else ""
strict_implementation_qc <- if (length(args) >= 4L) {
  normalized <- tolower(trimws(args[[4L]]))
  if (!normalized %in% c("true", "false")) {
    stop("strict_implementation_qc must be TRUE or FALSE.", call. = FALSE)
  }
  identical(normalized, "true")
} else {
  FALSE
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
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "simulation_qc_policy.R"))

aggregate_summary_directory <- file.path(root, "summary")
raw_path <- file.path(aggregate_summary_directory, "all_raw.csv")
if (file.exists(raw_path)) {
  raw <- read.csv(raw_path, stringsAsFactors = FALSE)
} else {
  raw_directory <- file.path(root, "raw")
  raw_files <- list.files(
    raw_directory, pattern = "^task_[0-9]+\\.csv$", full.names = TRUE
  )
  if (length(raw_files) == 0L) {
    stop("no aggregate or task-level raw results found under ", root)
  }
  raw <- RoCE:::.read_simulation_result_files(raw_files)
}

required_raw_columns <- c(
  "experiment", "dgp_type", "outcome_family", "config", "heterogeneity_type",
  "estimand_type", "p", "K", "rho", "cutoff", "n_site",
  "aggregation_lambda", "aggregation_cutoff",
  "n_folds", "nlambda_init", "nuisance_lambda_rule", "n_bootstrap",
  "M_tau", "M_tau_inference", "target_anchor_nlambda", "sim_id", "method",
  "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
  "n_site_arm_outcome_cells_below_8", "dr_weight_n",
  "dr_weight_n_clipped", "dr_weight_fraction_clipped",
  "dr_weight_max_site_fraction_clipped",
  "dr_weight_min_before_clipping", "dr_weight_max_before_clipping",
  "comparison_variance_method", "comparison_se_analytic",
  "comparison_se_bootstrap", "comparison_n_bootstrap",
  roce_resource_metadata_columns(),
  "workflow_fingerprint", "manifest_fingerprint",
  "estimand_scope", "estimate",
  "se", "bias", "coverage", "truth", "ci_lower", "ci_upper"
)
missing_raw_columns <- setdiff(required_raw_columns, names(raw))
if (length(missing_raw_columns) > 0L) {
  stop(
    "TATE results are missing required columns: ",
    paste(missing_raw_columns, collapse = ", "),
    call. = FALSE
  )
}
if (any(!is.finite(raw$target_anchor_nlambda)) ||
    any(as.integer(raw$target_anchor_nlambda) != 100L)) {
  stop("target_anchor_nlambda must equal the audited 100-point path.",
       call. = FALSE)
}
if (anyNA(raw$nuisance_lambda_rule) ||
    any(raw$nuisance_lambda_rule != "min")) {
  stop("production nuisance_lambda_rule must equal 'min'.", call. = FALSE)
}
expected_n_folds <- suppressWarnings(as.integer(Sys.getenv("ROCE_N_FOLDS", "5")))
if (length(expected_n_folds) != 1L || is.na(expected_n_folds) ||
    expected_n_folds < 2L) {
  stop("ROCE_N_FOLDS must be one integer of at least 2.", call. = FALSE)
}
RoCE:::.validate_face_production_scientific_metadata(
  raw, expected_n_folds = expected_n_folds
)
roce_validate_production_resource_metadata(raw)
provenance <- roce_result_provenance(raw, project_library)
observed_workflows <- unique(tolower(trimws(
  as.character(raw$workflow_fingerprint)
)))
observed_manifests <- unique(tolower(trimws(
  as.character(raw$manifest_fingerprint)
)))
observed_workflows <- observed_workflows[
  !is.na(observed_workflows) & nzchar(observed_workflows)
]
observed_manifests <- observed_manifests[
  !is.na(observed_manifests) & nzchar(observed_manifests)
]
expected_workflow <- tolower(trimws(Sys.getenv(
  "ROCE_WORKFLOW_FINGERPRINT", ""
)))
expected_manifest <- tolower(trimws(Sys.getenv(
  "ROCE_MANIFEST_FINGERPRINT", ""
)))
workflow_provenance_passed <-
  length(observed_workflows) == 1L &&
  grepl("^[0-9a-f]{64}$", observed_workflows) &&
  (!nzchar(expected_workflow) ||
     identical(observed_workflows, expected_workflow))
manifest_provenance_passed <-
  length(observed_manifests) == 1L &&
  grepl("^[0-9a-f]{64}$", observed_manifests) &&
  (!nzchar(expected_manifest) ||
     identical(observed_manifests, expected_manifest)) &&
  (!strict_implementation_qc ||
     grepl("^[0-9a-f]{64}$", expected_manifest))

diagnostics <- RoCE:::diagnose_simulation_results(
  raw, expected_replications = expected_replications
)
diagnostics <- RoCE:::add_simulation_rmse_comparisons(diagnostics)
diagnostics <- roce_augment_simulation_qc(
  diagnostics,
  parameter_bound = get("PARAM_MAX", envir = asNamespace("RoCE"))
)
job_label <- gsub(
  "[^A-Za-z0-9_-]", "_", Sys.getenv("SLURM_JOB_ID", "local")
)
run_label <- sprintf(
  "%s_job%s", format(Sys.time(), "%Y%m%d_%H%M%S"), job_label
)
final_diagnostic_directory <- file.path(
  root, "diagnostic_audits", run_label
)
dir.create(dirname(final_diagnostic_directory), recursive = TRUE,
           showWarnings = FALSE)
diagnostic_directory <- tempfile(
  pattern = paste0(".", run_label, "_"),
  tmpdir = dirname(final_diagnostic_directory)
)
if (!dir.create(diagnostic_directory, showWarnings = FALSE)) {
  stop("failed to create diagnostic staging directory.", call. = FALSE)
}
diagnostic_committed <- FALSE
on.exit({
  if (!diagnostic_committed && dir.exists(diagnostic_directory)) {
    unlink(diagnostic_directory, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)
write.csv(
  diagnostics,
  file.path(diagnostic_directory, "setting_method_diagnostics.csv"),
  row.names = FALSE
)

flagged <- diagnostics[diagnostics$diagnostic_status != "ok", , drop = FALSE]
write.csv(
  flagged,
  file.path(diagnostic_directory, "flagged_setting_methods.csv"),
  row.names = FALSE
)

paired_mse <- RoCE:::summarize_tate_paired_mse_benchmarks(raw)
write.csv(
  paired_mse,
  file.path(diagnostic_directory, "direct_tate_paired_mse_audit.csv"),
  row.names = FALSE
)

rho_zero <- if ("rho" %in% names(diagnostics)) {
  diagnostics[is.finite(diagnostics$rho) & diagnostics$rho == 0, , drop = FALSE]
} else {
  diagnostics[0L, , drop = FALSE]
}
rho_zero <- rho_zero[order(
  rho_zero$config, rho_zero$p, rho_zero$K, rho_zero$cutoff,
  rho_zero$rmse_rank, rho_zero$method
), , drop = FALSE]
write.csv(
  rho_zero,
  file.path(diagnostic_directory, "rho0_rmse_ranking.csv"),
  row.names = FALSE
)

missing_setting_methods <- data.frame()
if (nzchar(manifest_path)) {
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  output_methods <- function(methods) {
    requested <- trimws(strsplit(methods, ",", fixed = TRUE)[[1L]])
    requested <- requested[nzchar(requested)]
    labels <- paste0(requested, "_ate")
    crossfit <- requested %in% c("one_round_crossfit", "two_round_crossfit")
    c(labels, paste0(requested[crossfit], "_ate_armwise"))
  }
  expand_manifest_row <- function(index) {
    row <- manifest[index, , drop = FALSE]
    cutoffs <- if ("cutoffs" %in% names(row)) {
      as.numeric(strsplit(row$cutoffs, ",", fixed = TRUE)[[1L]])
    } else {
      as.numeric(row$cutoff)
    }
    methods <- if ("methods" %in% names(row)) {
      output_methods(row$methods)
    } else {
      c("one_round_crossfit_ate", "target_only_ate")
    }
    expected <- expand.grid(
      experiment = row$experiment,
      config = row$config,
      p = as.integer(row$p),
      K = as.integer(row$K),
      rho = as.numeric(row$rho),
      cutoff = cutoffs,
      n_site = as.integer(row$n_site),
      n_folds = as.integer(row$n_folds),
      method = methods,
      stringsAsFactors = FALSE
    )
    for (column in intersect(
      c("nlambda_init", "n_bootstrap", "M_tau", "M_tau_inference"),
      names(row)
    )) {
      expected[[column]] <- as.numeric(row[[column]])
    }
    expected
  }
  expected_setting_methods <- unique(do.call(
    rbind, lapply(seq_len(nrow(manifest)), expand_manifest_row)
  ))
  comparison_columns <- intersect(
    names(expected_setting_methods), names(raw)
  )
  observed_setting_methods <- unique(raw[comparison_columns])
  make_key <- function(data) {
    do.call(
      paste,
      c(lapply(data, function(x) format(x, scientific = FALSE, trim = TRUE)),
        sep = "\r")
    )
  }
  expected_key <- make_key(expected_setting_methods[comparison_columns])
  observed_key <- make_key(observed_setting_methods[comparison_columns])
  missing_setting_methods <- expected_setting_methods[
    !expected_key %in% observed_key, , drop = FALSE
  ]
  write.csv(
    missing_setting_methods,
    file.path(diagnostic_directory, "missing_setting_methods.csv"),
    row.names = FALSE
  )
}

coverage_flags <- grepl("coverage_(below|above)_mc_band",
                        diagnostics$diagnostic_status)
nuisance_flags <- grepl(
  paste(
    c(
      "nuisance_nonconvergence_detected", "invalid_nuisance_diagnostics",
      "nuisance_cv_candidates_excluded",
      "nuisance_line_search_failures_detected",
      "nuisance_density_ratio_boundary_detected",
      "nuisance_outcome_boundary_detected"
    ),
    collapse = "|"
  ),
  diagnostics$diagnostic_status
)
weight_flags <- grepl(
  "extreme_aggregation_weight",
  diagnostics$diagnostic_status
)
qc_classification <- roce_classify_simulation_qc(diagnostics)
implementation_flags <- qc_classification$implementation_failed
review_flags <- qc_classification$statistical_review_required
incomplete_flags <- grepl(
  "incomplete_replications|excess_replications",
  diagnostics$diagnostic_status
)

report_lines <- c(
  sprintf("Rows read: %d", nrow(raw)),
  sprintf("Package provenance valid: %s", provenance$passed),
  sprintf("Package provenance: %s", provenance$detail),
  sprintf("Workflow provenance valid: %s", workflow_provenance_passed),
  sprintf("Workflow fingerprint: %s", paste(observed_workflows, collapse = ",")),
  sprintf("Manifest provenance valid: %s", manifest_provenance_passed),
  sprintf("Expected manifest fingerprint: %s",
          if (nzchar(expected_manifest)) expected_manifest else "not supplied"),
  sprintf("Distinct manifest fingerprints: %d", length(observed_manifests)),
  sprintf("Setting-method combinations: %d", nrow(diagnostics)),
  sprintf("Paired TATE MSE comparisons: %d", nrow(paired_mse)),
  sprintf("Expected replications per combination: %d", expected_replications),
  sprintf("Incomplete/excess combinations: %d", sum(incomplete_flags)),
  sprintf("Implementation-QC flagged combinations: %d",
          sum(implementation_flags)),
  sprintf("Statistical-review combinations: %d", sum(review_flags)),
  sprintf("Density-ratio boundary combinations: %d",
          sum(diagnostics$nuisance_density_ratio_boundary_detected)),
  sprintf("Outcome-model boundary combinations: %d",
          sum(diagnostics$nuisance_outcome_boundary_detected)),
  sprintf("Sparse binary-cell review combinations: %d",
          sum(diagnostics$sparse_binary_cell_replications > 0,
              na.rm = TRUE)),
  sprintf("Density-ratio clipping review combinations: %d",
          sum(diagnostics$dr_weight_clipped_replications > 0,
              na.rm = TRUE)),
  sprintf("Nuisance-fit flagged combinations: %d", sum(nuisance_flags)),
  sprintf("Extreme-weight review combinations: %d", sum(weight_flags)),
  sprintf("Coverage outside nominal Monte Carlo band: %d",
          sum(coverage_flags)),
  sprintf("Missing setting-method combinations from manifest: %d",
          nrow(missing_setting_methods)),
  "",
  paste0(
    "Coverage flags are statistical review signals, not automatically code ",
    "failures. Inspect bias, SE/empirical-SD, and normal-reference coverage ",
    "together for every flagged setting."
  )
)
writeLines(
  report_lines, file.path(diagnostic_directory, "diagnostic_report.txt")
)
if (!file.rename(diagnostic_directory, final_diagnostic_directory)) {
  stop("failed to atomically commit the diagnostic audit directory.",
       call. = FALSE)
}
diagnostic_committed <- TRUE
message(paste(report_lines, collapse = "\n"))
message("Diagnostic artifacts: ", final_diagnostic_directory)

if (strict_implementation_qc &&
    (!provenance$passed || !workflow_provenance_passed ||
       !manifest_provenance_passed ||
       any(incomplete_flags) || any(implementation_flags) ||
       nrow(missing_setting_methods) > 0L)) {
  stop(
    paste0(
      "strict implementation QC failed; inspect ",
      final_diagnostic_directory,
      "/diagnostic_report.txt and the machine-readable diagnostic CSVs."
    ),
    call. = FALSE
  )
}
}

main()
