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
raw_directory <- file.path(root, "raw")
final_summary_directory <- file.path(root, "summary")
if (dir.exists(final_summary_directory) ||
    file.exists(final_summary_directory)) {
  stop(
    "refusing to overwrite existing aggregate summary: ",
    final_summary_directory, call. = FALSE
  )
}
dir.create(dirname(final_summary_directory), recursive = TRUE,
           showWarnings = FALSE)
summary_directory <- tempfile(
  pattern = ".direct_tate_summary_",
  tmpdir = dirname(final_summary_directory)
)
if (!dir.create(summary_directory, showWarnings = FALSE)) {
  stop("failed to create aggregate summary staging directory.",
       call. = FALSE)
}
summary_committed <- FALSE
on.exit({
  if (!summary_committed && dir.exists(summary_directory)) {
    unlink(summary_directory, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)
files <- list.files(
  raw_directory, pattern = "^task_[0-9]+\\.csv$", full.names = TRUE
)
if (length(files) == 0L) {
  stop("no task CSV files found in ", raw_directory)
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
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
source(file.path("scripts", "slurm", "simulation_qc_policy.R"))
z_alpha_05 <- stats::qnorm(0.975)
raw <- RoCE:::.read_simulation_result_files(files)
RoCE:::.validate_face_production_scientific_metadata(raw)
provenance <- roce_result_provenance(raw, project_library)
if (!provenance$passed) {
  stop("aggregate input failed package provenance: ", provenance$detail,
       call. = FALSE)
}
roce_validate_production_resource_metadata(raw)
roce_validate_scheduler_provenance(raw, require_slurm = TRUE)
expected_primary_cutoff <- 1 / get(
  "AGG_WALD_LAMBDA", envir = asNamespace("RoCE")
)
required_provenance_columns <- c(
  "primary_cutoff", "workflow_fingerprint", "manifest_fingerprint"
)
missing_provenance_columns <- setdiff(
  required_provenance_columns, names(raw)
)
if (length(missing_provenance_columns) > 0L ||
    any(abs(raw$primary_cutoff - expected_primary_cutoff) > 1e-12)) {
  stop("aggregate input has missing or inconsistent primary-cutoff metadata.",
       call. = FALSE)
}
for (column in c("workflow_fingerprint", "manifest_fingerprint")) {
  values <- unique(tolower(trimws(as.character(raw[[column]]))))
  if (length(values) != 1L || !grepl("^[0-9a-f]{64}$", values)) {
    stop("aggregate input has mixed or malformed ", column, ".",
         call. = FALSE)
  }
}
write.csv(raw, file.path(summary_directory, "all_raw.csv"), row.names = FALSE)

group_names <- intersect(c(
  "experiment", "dgp_type", "outcome_family", "config", "heterogeneity_type",
  "estimand_type", "p", "K", "rho", "cutoff", "n_site", "n_folds",
  "nlambda_init", "nuisance_lambda_rule", "n_bootstrap", "M_tau",
  "M_tau_inference",
  roce_resource_metadata_columns()
), names(raw))
group_key <- RoCE:::.simulation_group_key(raw, group_names)
summaries <- lapply(split(raw, group_key), function(group) {
  summary <- summarize_results(group)
  for (name in group_names) {
    summary[[name]] <- group[[name]][[1L]]
  }
  summary$bias_mcse <- summary$bias_sd / sqrt(summary$n_success)
  summary$coverage_mcse <- sqrt(
    summary$coverage * (1 - summary$coverage) / summary$n_success
  )
  summary$se_to_empirical_sd <- ifelse(
    summary$bias_sd > 0,
    summary$se_mean / summary$bias_sd,
    NA_real_
  )
  # Normal-reference coverage using the observed Monte Carlo center and spread.
  # Comparing this with empirical coverage distinguishes center bias from a
  # mis-scaled standard error without requiring raw replicate inspection.
  summary$normal_reference_coverage <- ifelse(
    summary$bias_sd > 0,
    stats::pnorm(
      (z_alpha_05 * summary$se_mean - summary$bias_mean) / summary$bias_sd
    ) - stats::pnorm(
      (-z_alpha_05 * summary$se_mean - summary$bias_mean) / summary$bias_sd
    ),
    NA_real_
  )

  diagnostic_columns <- intersect(
    unique(c(
      "target_anchor_weight", "mean_abs_source_weight",
      "max_abs_source_weight",
      "max_wald_statistic", "mean_wald_statistic",
      "penalized_source_fold_fraction",
      "max_weight_optimizer_iterations", "max_weight_psd_ridge",
      "weight_psd_ridge_fold_fraction", "task_elapsed_seconds",
      "inference_logit_truncated",
      "inference_logit_truncation_fraction",
      "inference_max_abs_logit", "inference_safety_clip_count",
      "comparison_se_analytic", "comparison_se_bootstrap",
      "comparison_n_bootstrap",
      "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
      "n_site_arm_outcome_cells_below_8", "dr_weight_n",
      "dr_weight_n_clipped", "dr_weight_fraction_clipped",
      "dr_weight_max_site_fraction_clipped",
      "dr_weight_min_before_clipping", "dr_weight_max_before_clipping",
      "data_generation_seconds", "one_round_face_wall_seconds",
      "two_round_face_wall_seconds", "comparison_mu1_seconds",
      "tate_contrast_stage_seconds", "simulation_elapsed_seconds",
      # Backward-compatible readers for completed pre-rename smoke artifacts.
      "one_round_mu1_seconds", "two_round_mu1_seconds",
      "direct_tate_seconds",
      grep(
        "^face_.*_(seconds|nonconverged|iterations|update_ratio|coefficient)$",
        names(group), value = TRUE
      )
    )),
    names(group)
  )
  for (column in diagnostic_columns) {
    summary[[paste0(column, "_mean")]] <- vapply(
      as.character(summary$method),
      function(method_name) {
        values <- group[group$method == method_name, column]
        if (length(values) == 0L || all(is.na(values))) {
          return(NA_real_)
        }
        mean(values, na.rm = TRUE)
      },
      numeric(1L)
    )
  }
  summary
})
summary <- do.call(rbind, summaries)
write.csv(
  summary,
  file.path(summary_directory, "all_summary.csv"),
  row.names = FALSE
)

diagnostics <- RoCE:::diagnose_simulation_results(
  raw, expected_replications = expected_replications
)
diagnostics <- RoCE:::add_simulation_rmse_comparisons(diagnostics)
diagnostics <- roce_augment_simulation_qc(
  diagnostics,
  parameter_bound = get("PARAM_MAX", envir = asNamespace("RoCE"))
)
write.csv(
  diagnostics,
  file.path(summary_directory, "setting_method_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  diagnostics[diagnostics$diagnostic_status != "ok", , drop = FALSE],
  file.path(summary_directory, "flagged_setting_methods.csv"),
  row.names = FALSE
)
paired_mse <- RoCE:::summarize_tate_paired_mse_benchmarks(raw)
write.csv(
  paired_mse,
  file.path(summary_directory, "direct_tate_paired_mse_audit.csv"),
  row.names = FALSE
)
if ("rho" %in% names(diagnostics)) {
  rho_zero <- diagnostics[
    is.finite(diagnostics$rho) & diagnostics$rho == 0, , drop = FALSE
  ]
  rho_zero <- rho_zero[order(
    rho_zero$config, rho_zero$p, rho_zero$K, rho_zero$cutoff,
    rho_zero$rmse_rank, rho_zero$method
  ), , drop = FALSE]
  write.csv(
    rho_zero,
    file.path(summary_directory, "rho0_rmse_ranking.csv"),
    row.names = FALSE
  )
}

for (setting in split(summary, interaction(
  summary$experiment, summary$config, summary$p, summary$K,
  drop = TRUE, lex.order = TRUE
))) {
  filename <- sprintf(
    "negtransfer_%s_p%d_K%d.csv",
    setting$config[[1L]], setting$p[[1L]], setting$K[[1L]]
  )
  write.csv(setting, file.path(summary_directory, filename), row.names = FALSE)
}
if (!file.rename(summary_directory, final_summary_directory)) {
  stop("failed to atomically commit the aggregate summary directory.",
       call. = FALSE)
}
summary_committed <- TRUE
message(
  "aggregated ", length(files), " task files into ",
  final_summary_directory
)
}

main()
