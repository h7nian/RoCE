#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
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
if (is.na(expected_replications) || expected_replications < 1L ||
    expected_replications > 500L) {
  stop("expected_replications must be an integer in 1:500.", call. = FALSE)
}

default_library <- file.path(
  normalizePath(".", mustWork = TRUE),
  "results", "direct_tate_mc500_b5000", "Rlib_current"
)
project_library <- Sys.getenv(
  "ROCE_PROJECT_LIB",
  if (dir.exists(default_library)) default_library else ""
)
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "simulation_qc_policy.R"))
source(file.path("scripts", "slurm", "atomic_output.R"))

read_task_directory <- function(path, label) {
  files <- list.files(
    path, pattern = "^task_[0-9]+[.]csv$", full.names = TRUE
  )
  if (length(files) == 0L) {
    stop("no ", label, " task files found in ", path, call. = FALSE)
  }
  RoCE:::.read_simulation_result_files(files)
}

sidecar <- read_task_directory(
  file.path(root, "reused_sensitivity_raw"), "reused-sensitivity"
)
fitting <- read_task_directory(
  file.path(root, "fitting_radius_diagnostic", "raw"), "fitting-radius"
)
sensitivity_inputs <- list(sidecar = sidecar, fitting = fitting)
expected_cv_threads <- roce_read_positive_integer_env(
  "ROCE_EXPECT_NUISANCE_CV_THREADS", 5L, 5L
)
invisible(lapply(sensitivity_inputs, function(item) {
  roce_validate_result_resource_metadata(
    item, expected_threads = expected_cv_threads
  )
}))
provenance_records <- lapply(sensitivity_inputs, function(item) {
  provenance <- roce_result_provenance(item, project_library)
  if (!provenance$passed) {
    stop("sensitivity package provenance failed: ", provenance$detail,
         call. = FALSE)
  }
  provenance
})
if (length(unique(vapply(
      provenance_records, `[[`, character(1L), "library"
    ))) != 1L ||
    length(unique(vapply(
      provenance_records, `[[`, character(1L), "fingerprint"
    ))) != 1L) {
  stop(
    "reused and fitting-radius sensitivities use different package builds.",
    call. = FALSE
  )
}
workflow_values <- lapply(sensitivity_inputs, function(item) {
  if (!"workflow_fingerprint" %in% names(item)) character(0) else {
    unique(tolower(trimws(as.character(item$workflow_fingerprint))))
  }
})
if (any(lengths(workflow_values) != 1L) ||
    length(unique(unlist(workflow_values, use.names = FALSE))) != 1L ||
    !grepl("^[0-9a-f]{64}$", workflow_values[[1L]])) {
  stop(
    "reused and fitting-radius sensitivities do not share one valid workflow fingerprint.",
    call. = FALSE
  )
}
if (!"primary_experiment" %in% names(sidecar)) {
  stop(
    "reused-sensitivity rows lack primary_experiment provenance.",
    call. = FALSE
  )
}
unexpected_primary_experiments <- setdiff(
  unique(as.character(sidecar$primary_experiment)), "negative_transfer"
)
if (length(unexpected_primary_experiments) > 0L) {
  stop(
    paste0(
      "production reused-sensitivity directory contains non-production rows: ",
      paste(unexpected_primary_experiments, collapse = ", ")
    ),
    call. = FALSE
  )
}

sidecar_required_metadata <- c(
  "experiment", "config", "p", "K", "rho", "cutoff",
  "aggregation_lambda", "aggregation_cutoff", "dgp_type", "outcome_family",
  "heterogeneity_type", "estimand_type", "n_site", "n_folds",
  "nlambda_init", "nuisance_lambda_rule", "n_bootstrap", "M_tau",
  "M_tau_inference", "primary_cutoff",
  "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
  "n_site_arm_outcome_cells_below_8",
  "reused_primary_task", "nuisance_refit_count", "sensitivity_kind",
  "inference_refresh_count", roce_resource_metadata_columns()
)
fitting_required_metadata <- c(
  "experiment", "config", "p", "K", "rho", "cutoff", "n_site",
  "n_folds", "nlambda_init", "n_bootstrap", "M_tau",
  "M_tau_inference", roce_resource_metadata_columns()
)
missing_sidecar_metadata <- setdiff(sidecar_required_metadata, names(sidecar))
missing_fitting_metadata <- setdiff(fitting_required_metadata, names(fitting))
if (length(missing_sidecar_metadata) > 0L ||
    length(missing_fitting_metadata) > 0L) {
  stop(
    paste0(
      "sensitivity metadata schema is incomplete; sidecar={",
      paste(missing_sidecar_metadata, collapse = ","), "}; fitting={",
      paste(missing_fitting_metadata, collapse = ","), "}."
    ),
    call. = FALSE
  )
}
primary_cutoffs <- unique(suppressWarnings(as.numeric(
  sidecar$primary_cutoff
)))
primary_cutoffs <- primary_cutoffs[is.finite(primary_cutoffs)]
package_primary_cutoff <- 1 / get(
  "AGG_WALD_LAMBDA", envir = asNamespace("RoCE")
)
if (length(primary_cutoffs) != 1L || primary_cutoffs <= 0 ||
    !isTRUE(all.equal(primary_cutoffs, package_primary_cutoff,
                      tolerance = 1e-12))) {
  stop(
    paste0(
      "reused sensitivity must record one primary cutoff matching the ",
      "audited package default."
    ),
    call. = FALSE
  )
}
primary_cutoff <- primary_cutoffs[[1L]]
sidecar_metadata_ok <-
  sidecar$experiment == "c3_reused_sensitivity" &
  sidecar$config == "C3" & sidecar$p == 100L & sidecar$K == 4L &
  sidecar$rho %in% c(0, 2.5) & sidecar$n_site == 1000L &
  sidecar$n_folds == 5L & sidecar$n_bootstrap == 5000L &
  sidecar$M_tau == 5 & sidecar$dgp_type == "face" &
  sidecar$outcome_family == "binomial" &
  sidecar$estimand_type == "superpopulation" &
  sidecar$nuisance_lambda_rule == "min" &
  sidecar$aggregation_cutoff == sidecar$cutoff &
  abs(sidecar$aggregation_lambda - 1 / sidecar$cutoff) <= 1e-12 &
  is.finite(sidecar$min_site_arm_outcome_cell_n) &
  is.finite(sidecar$min_target_arm_outcome_cell_n) &
  is.finite(sidecar$n_site_arm_outcome_cells_below_8) &
  sidecar$reused_primary_task &
  sidecar$nuisance_refit_count == 0L &
  sidecar$inference_refresh_count == 3L &
  abs(sidecar$primary_cutoff - primary_cutoff) <= 1e-12
expected_sidecar_heterogeneity <- ifelse(
  sidecar$rho > 0, "one_deviated_source", "none"
)
sidecar_metadata_ok <- sidecar_metadata_ok &
  sidecar$heterogeneity_type == expected_sidecar_heterogeneity
expected_sensitivity_kind <- ifelse(
  sidecar$M_tau_inference != 5,
  "inference_radius",
  ifelse(
    abs(sidecar$cutoff - primary_cutoff) > 1e-12,
    "aggregation_cutoff",
    "reference_identity"
  )
)
sidecar_metadata_ok <- sidecar_metadata_ok &
  sidecar$sensitivity_kind == expected_sensitivity_kind
if (!all(!is.na(sidecar_metadata_ok) & sidecar_metadata_ok)) {
  stop(
    "reused-sensitivity rows do not match the audited p=100 C3 design/reuse metadata.",
    call. = FALSE
  )
}
fitting_metadata_ok <-
  fitting$experiment == "c3_truncation_diagnostic" &
  fitting$config == "C3" & fitting$p == 100L & fitting$K == 4L &
  fitting$rho %in% c(0, 2.5) &
  abs(fitting$cutoff - primary_cutoff) <= 1e-12 &
  fitting$n_site == 1000L & fitting$n_folds == 5L &
  fitting$n_bootstrap == 5000L & fitting$M_tau %in% c(4, 6) &
  fitting$M_tau_inference == 5
if (!all(!is.na(fitting_metadata_ok) & fitting_metadata_ok)) {
  stop(
    "fitting-radius rows do not match the audited p=100 C3 design.",
    call. = FALSE
  )
}
sidecar_nlambda <- unique(as.integer(sidecar$nlambda_init))
fitting_nlambda <- unique(as.integer(fitting$nlambda_init))
if (length(sidecar_nlambda) != 1L || length(fitting_nlambda) != 1L ||
    !identical(sidecar_nlambda, fitting_nlambda)) {
  stop(
    "reused and fitting-radius sensitivities do not share one nuisance-grid size.",
    call. = FALSE
  )
}

methods <- c("one_round_crossfit_ate", "target_only_ate")
sidecar <- sidecar[
  sidecar$estimand_scope == "tate" & sidecar$method %in% methods,
  , drop = FALSE
]
fitting <- fitting[
  fitting$estimand_scope == "tate" & fitting$method %in% methods,
  , drop = FALSE
]

assert_complete_grid <- function(data, setting_columns, n_settings, label) {
  required <- c(
    "sim_id", "method", "estimate", "se", "truth", "bias", "coverage",
    "ci_lower", "ci_upper", "ci_width", setting_columns
  )
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(label, " is missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  setting_key <- interaction(
    data[setting_columns], drop = TRUE, lex.order = TRUE
  )
  cells <- split(data, interaction(
    setting_key, data$method, drop = TRUE, lex.order = TRUE
  ))
  expected_ids <- seq_len(expected_replications)
  valid <- vapply(cells, function(cell) {
    identical(sort(as.integer(cell$sim_id)), expected_ids) &&
      !anyDuplicated(cell$sim_id)
  }, logical(1L))
  if (length(unique(setting_key)) != n_settings ||
      length(cells) != n_settings * length(methods) || !all(valid)) {
    stop(
      sprintf(
        "%s is incomplete: settings=%d, method-cells=%d, invalid-cells=%d.",
        label, length(unique(setting_key)), length(cells), sum(!valid)
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

cutoff <- sidecar[
  sidecar$M_tau == 5 & sidecar$M_tau_inference == 5,
  , drop = FALSE
]
assert_complete_grid(
  cutoff,
  c("config", "p", "K", "rho", "cutoff", "M_tau", "M_tau_inference"),
  n_settings = 10L,
  label = "cutoff sensitivity"
)

sidecar_truncation <- sidecar[
  abs(sidecar$cutoff - primary_cutoff) <= 1e-12, , drop = FALSE
]
fitting_truncation <- fitting[
  abs(fitting$cutoff - primary_cutoff) <= 1e-12 &
    fitting$M_tau %in% c(4, 6) &
    fitting$M_tau_inference == 5,
  , drop = FALSE
]
truncation <- RoCE:::.bind_sim_result_list(list(
  sidecar_truncation, fitting_truncation
))
truncation$experiment <- "c3_truncation_sensitivity"
assert_complete_grid(
  truncation,
  c("config", "p", "K", "rho", "cutoff", "M_tau", "M_tau_inference"),
  n_settings = 12L,
  label = "combined truncation sensitivity"
)

summarize_sensitivity <- function(data, label) {
  diagnostics <- RoCE:::diagnose_simulation_results(
    data, expected_replications = expected_replications
  )
  diagnostics <- RoCE:::add_simulation_rmse_comparisons(diagnostics)
  diagnostics <- roce_augment_simulation_qc(
    diagnostics,
    parameter_bound = get("PARAM_MAX", envir = asNamespace("RoCE"))
  )
  flagged <- diagnostics[diagnostics$diagnostic_status != "ok", , drop = FALSE]
  qc_classification <- roce_classify_simulation_qc(diagnostics)
  implementation_failed <- qc_classification$implementation_failed
  statistical_review <- qc_classification$statistical_review_required
  report <- c(
    sprintf("Sensitivity: %s", label),
    sprintf("Raw rows: %d", nrow(data)),
    sprintf("Setting-method rows: %d", nrow(diagnostics)),
    sprintf("Implementation failures: %d", sum(implementation_failed)),
    sprintf(
      "Coverage/RMSE review rows: %d",
      sum(statistical_review & !implementation_failed)
    )
  )
  list(
    label = label,
    raw = data,
    diagnostics = diagnostics,
    flagged = flagged,
    report = report,
    implementation_failed = implementation_failed
  )
}

summary_root <- file.path(root, "sensitivity_summary")
evaluations <- list(
  cutoff = summarize_sensitivity(cutoff, "cutoff"),
  truncation = summarize_sensitivity(truncation, "truncation")
)
roce_write_atomic_directory(
  summary_root,
  writer = function(staging_directory) {
    for (name in names(evaluations)) {
      evaluation <- evaluations[[name]]
      output_directory <- file.path(staging_directory, name)
      dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
      write.csv(
        evaluation$raw, file.path(output_directory, "all_raw.csv"),
        row.names = FALSE
      )
      write.csv(
        evaluation$diagnostics,
        file.path(output_directory, "setting_method_diagnostics.csv"),
        row.names = FALSE
      )
      write.csv(
        evaluation$flagged,
        file.path(output_directory, "flagged_setting_methods.csv"),
        row.names = FALSE
      )
      writeLines(
        evaluation$report,
        file.path(output_directory, "diagnostic_report.txt")
      )
    }
  },
  caller = "reused sensitivity aggregation"
)
failed_evaluations <- names(evaluations)[vapply(
  evaluations,
  function(evaluation) any(evaluation$implementation_failed),
  logical(1L)
)]
if (length(failed_evaluations) > 0L) {
  stop(
    "sensitivity implementation QC failed for: ",
    paste(failed_evaluations, collapse = ", "),
    "; inspect ", summary_root, call. = FALSE
  )
}
message(sprintf(
  paste0(
    "sensitivity audit passed: cutoff setting-methods=%d; ",
    "truncation setting-methods=%d"
  ),
  nrow(evaluations$cutoff$diagnostics),
  nrow(evaluations$truncation$diagnostics)
))
