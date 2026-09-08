#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 6L) {
  stop(
    paste(
      "usage: summarize_grouped_cutoff_diagnostic.R",
      paste0(
        "RAW_DIRECTORY MANIFEST.csv EXPECTED_REPLICATIONS ",
        paste0(
          "[PRIMARY_RAW_DIRECTORY] [PRIMARY_MANIFEST.csv] ",
          "[require-production-identity|disjoint-pilot]"
        )
      )
    ),
    call. = FALSE
  )
}
raw_directory <- normalizePath(args[[1L]], mustWork = TRUE)
manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
expected_replications <- suppressWarnings(as.integer(args[[3L]]))
if (is.na(expected_replications) || expected_replications < 1L ||
    expected_replications > 25L) {
  stop("EXPECTED_REPLICATIONS must be an integer in 1:25.", call. = FALSE)
}
results_root <- dirname(dirname(raw_directory))
production_identity_mode <- if (length(args) >= 6L) {
  args[[6L]]
} else {
  "require-production-identity"
}
if (!production_identity_mode %in%
    c("require-production-identity", "disjoint-pilot")) {
  stop(
    paste0(
      "production identity mode must be require-production-identity or ",
      "disjoint-pilot."
    ),
    call. = FALSE
  )
}
primary_raw_directory <- normalizePath(
  if (length(args) >= 4L) args[[4L]] else file.path(results_root, "raw"),
  mustWork = TRUE
)
primary_manifest_path <- normalizePath(
  if (length(args) >= 5L) {
    args[[5L]]
  } else {
    file.path(results_root, "manifest_main.csv")
  },
  mustWork = TRUE
)

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "atomic_output.R"))
source(file.path("scripts", "slurm", "grouped_cutoff_diagnostic_helpers.R"))
primary_cutoff <- 1 / get(
  "AGG_WALD_LAMBDA", envir = asNamespace("RoCE")
)
summary_workflow_fingerprint <- roce_sha256_environment(
  "ROCE_CUTOFF_SUMMARY_WORKFLOW_FINGERPRINT"
)

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
manifest_required <- c(
  "task_id", "sim_id", "config", "p", "K", "rho_values", "cutoffs",
  "n_site", "n_folds", "nlambda_init", "n_bootstrap", "M_tau",
  "M_tau_inference"
)
missing_manifest_columns <- setdiff(manifest_required, names(manifest))
if (length(missing_manifest_columns) > 0L) {
  stop("manifest is missing: ", paste(missing_manifest_columns, collapse = ", "),
       call. = FALSE)
}
manifest <- manifest[order(manifest$sim_id), , drop = FALSE]
if (nrow(manifest) < expected_replications) {
  stop("manifest is shorter than the requested simulation prefix.",
       call. = FALSE)
}
manifest <- manifest[seq_len(expected_replications), , drop = FALSE]
expected_simulation_ids <- seq.int(
  as.integer(manifest$sim_id[[1L]]), length.out = expected_replications
)
if (!identical(as.integer(manifest$sim_id), expected_simulation_ids)) {
  stop("manifest does not contain a contiguous requested simulation prefix.",
       call. = FALSE)
}
if (identical(production_identity_mode, "disjoint-pilot") &&
    any(manifest$sim_id <= 500L)) {
  stop(
    "disjoint-pilot mode requires seeds above the primary 1:500 range.",
    call. = FALSE
  )
}
constant_manifest_columns <- c(
  "config", "p", "K", "rho_values", "cutoffs", "n_site", "n_folds",
  "nlambda_init", "n_bootstrap", "M_tau", "M_tau_inference"
)
if (any(vapply(manifest[constant_manifest_columns], function(value) {
  length(unique(value)) != 1L
}, logical(1L))) ||
    manifest$p[[1L]] != 100L ||
    !manifest$K[[1L]] %in% c(2L, 4L, 8L) ||
    manifest$n_site[[1L]] != 1000L ||
    manifest$n_folds[[1L]] != 5L ||
    manifest$nlambda_init[[1L]] != 100L ||
    manifest$n_bootstrap[[1L]] != 5000L ||
    manifest$M_tau[[1L]] != 5 ||
    manifest$M_tau_inference[[1L]] != 5) {
  stop("manifest mixes settings or departs from the locked p=100 production settings.",
       call. = FALSE)
}
files <- file.path(
  raw_directory, sprintf("task_%06d.csv", as.integer(manifest$task_id))
)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) {
  stop("grouped cutoff prefix is incomplete; first missing file: ",
       missing_files[[1L]], call. = FALSE)
}
raw <- RoCE:::.read_simulation_result_files(files)

required <- c(
  "task_id", "sim_id", "config", "p", "K", "rho", "cutoff", "method",
  "estimand_scope", "estimate", "se", "truth", "bias", "coverage",
  "ci_lower", "ci_upper", "ci_width", "aggregation_lambda",
  "aggregation_cutoff",
  "target_anchor_weight", "source1_mean_weight",
  "source1_mean_abs_weight", "source1_zero_weight_fold_fraction",
  "source1_mean_wald_statistic", "source1_max_wald_statistic",
  "source1_penalty_activation_fraction",
  "informative_sources_mean_abs_weight",
  "source1_estimate", "source1_minus_target_estimate",
  "source1_mean_training_discrepancy",
  "primary_cutoff_identity_checked", "primary_cutoff_identity_exact",
  "n_nuisance_refits_for_cutoff_grid", "package_library",
  "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint"
)
missing_columns <- setdiff(required, names(raw))
if (length(missing_columns) > 0L) {
  stop("diagnostic results are missing: ",
       paste(missing_columns, collapse = ", "), call. = FALSE)
}

finite_columns <- c(
  "task_id", "sim_id", "p", "K", "rho", "cutoff", "estimate", "se",
  "truth", "bias", "ci_lower", "ci_upper", "ci_width",
  "aggregation_lambda", "aggregation_cutoff", "target_anchor_weight"
)
face_finite_columns <- c(
  "source1_mean_weight", "source1_mean_abs_weight",
  "source1_zero_weight_fold_fraction", "source1_mean_wald_statistic",
  "source1_max_wald_statistic", "source1_penalty_activation_fraction",
  "source1_estimate", "source1_minus_target_estimate",
  "source1_mean_training_discrepancy"
)
finite_values <- suppressWarnings(as.matrix(data.frame(
  lapply(raw[finite_columns], as.numeric), check.names = FALSE
)))
face_rows_for_validation <- raw$method == "one_round_crossfit_ate"
face_finite_values <- suppressWarnings(as.matrix(data.frame(
  lapply(raw[face_rows_for_validation, face_finite_columns, drop = FALSE],
         as.numeric),
  check.names = FALSE
)))
coverage <- suppressWarnings(as.numeric(raw$coverage))
tolerance <- 1e-12
if (anyNA(finite_values) || any(!is.finite(finite_values)) ||
    anyNA(face_finite_values) || any(!is.finite(face_finite_values)) ||
    anyNA(coverage) || any(!coverage %in% c(0, 1)) ||
    any(raw$se <= 0) || any(raw$ci_width <= 0) ||
    any(abs(raw$bias - (raw$estimate - raw$truth)) > tolerance) ||
    any(abs(raw$ci_width - (raw$ci_upper - raw$ci_lower)) > tolerance) ||
    any(raw$aggregation_lambda <= 0) ||
    any(abs(raw$aggregation_lambda - 1 / raw$cutoff) > tolerance) ||
    any(abs(raw$aggregation_cutoff - raw$cutoff) > tolerance) ||
    any((raw$estimate >= raw$ci_lower & raw$estimate <= raw$ci_upper) != TRUE) ||
    any(coverage != as.numeric(
      raw$truth >= raw$ci_lower & raw$truth <= raw$ci_upper
    ))) {
  stop("grouped cutoff results fail finite-value or statistical identities.",
       call. = FALSE)
}

expected_rhos <- roce_parse_semicolon_numeric_grid(
  manifest$rho_values[[1L]], "rho_values", require_zero_first = TRUE
)
expected_cutoffs <- roce_parse_semicolon_numeric_grid(
  manifest$cutoffs[[1L]], "cutoffs", require_positive = TRUE
)
if (!any(abs(expected_cutoffs - primary_cutoff) <= 1e-12)) {
  stop("cutoff grid does not contain the tested package primary cutoff.",
       call. = FALSE)
}
if ("primary_cutoff" %in% names(manifest)) {
  manifest_primary_cutoff <- suppressWarnings(as.numeric(
    manifest$primary_cutoff
  ))
  if (anyNA(manifest_primary_cutoff) ||
      any(!is.finite(manifest_primary_cutoff)) ||
      any(abs(manifest_primary_cutoff - primary_cutoff) > 1e-12)) {
    stop("manifest primary cutoff differs from the tested package.",
         call. = FALSE)
  }
}
expected_methods <- c("one_round_crossfit_ate", "target_only_ate")
expected_rows_per_task <- length(expected_rhos) * length(expected_cutoffs) *
  length(expected_methods)
rows_by_task <- split(raw, raw$task_id)
schema_ok <- vapply(rows_by_task, function(rows) {
  manifest_row <- manifest[manifest$task_id == rows$task_id[[1L]], , drop = FALSE]
  nrow(rows) == expected_rows_per_task &&
    nrow(manifest_row) == 1L &&
    all(rows$sim_id == manifest_row$sim_id) &&
    all(rows$config == manifest_row$config) &&
    all(rows$p == manifest_row$p) &&
    all(rows$K == manifest_row$K) &&
    setequal(unique(rows$rho), expected_rhos) &&
    setequal(unique(rows$cutoff), expected_cutoffs) &&
    setequal(unique(rows$method), expected_methods) &&
    !anyDuplicated(rows[c("rho", "cutoff", "method")])
}, logical(1L))
if (length(rows_by_task) != expected_replications || !all(schema_ok)) {
  stop("one or more grouped cutoff tasks failed the exact row schema.",
       call. = FALSE)
}

target_rows <- raw[raw$method == "target_only_ate", , drop = FALSE]
target_groups <- split(
  target_rows, interaction(target_rows$task_id, target_rows$rho, drop = TRUE)
)
target_cutoff_invariant <- vapply(target_groups, function(rows) {
  invariant_columns <- c(
    "estimate", "se", "truth", "bias", "coverage", "ci_lower", "ci_upper",
    "ci_width"
  )
  all(vapply(rows[invariant_columns], function(values) {
    length(unique(values)) == 1L
  }, logical(1L)))
}, logical(1L))
if (length(target_groups) != expected_replications * length(expected_rhos) ||
    !all(target_cutoff_invariant)) {
  stop("target-only results changed across cutoff-only reaggregation.",
       call. = FALSE)
}

provenance <- roce_result_provenance(raw, project_library)
if (!provenance$passed) {
  stop("grouped cutoff provenance failed: ", provenance$detail,
       call. = FALSE)
}
workflow_fingerprints <- unique(tolower(trimws(raw$workflow_fingerprint)))
manifest_fingerprints <- unique(tolower(trimws(raw$manifest_fingerprint)))
if (length(workflow_fingerprints) != 1L ||
    !grepl("^[0-9a-f]{64}$", workflow_fingerprints[[1L]]) ||
    length(manifest_fingerprints) != 1L ||
    !grepl("^[0-9a-f]{64}$", manifest_fingerprints[[1L]])) {
  stop("workflow or manifest provenance is mixed or malformed.", call. = FALSE)
}
expected_manifest_sha256 <- roce_sha256_file(manifest_path)
if (!all(tolower(raw$manifest_fingerprint) == tolower(expected_manifest_sha256))) {
  stop("result manifest fingerprint does not match the audited manifest.",
       call. = FALSE)
}
if (any(raw$n_nuisance_refits_for_cutoff_grid != 0L)) {
  stop("cutoff grid recorded an unexpected nuisance refit.", call. = FALSE)
}
identity_rows <- raw$primary_cutoff_identity_checked %in% TRUE
if (sum(identity_rows) != expected_replications * length(expected_rhos) * 2L ||
    !all(raw$primary_cutoff_identity_exact[identity_rows] %in% TRUE)) {
  stop("primary cutoff identity audit is incomplete or failed.", call. = FALSE)
}

# The within-task identity above proves that cutoff-only reaggregation is an
# exact no-refit operation. This second gate is deliberately external: it
# verifies that the tested package's primary cutoff also reproduces the
# separately committed production result for the same setting and seed.
production_identity_fields <- character(0)
if (identical(
    production_identity_mode, "require-production-identity"
)) {
primary_manifest <- read.csv(primary_manifest_path, stringsAsFactors = FALSE)
primary_manifest_required <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho", "cutoff"
)
missing_primary_manifest_columns <- setdiff(
  primary_manifest_required, names(primary_manifest)
)
if (length(missing_primary_manifest_columns) > 0L) {
  stop(
    "primary manifest is missing: ",
    paste(missing_primary_manifest_columns, collapse = ", "), call. = FALSE
  )
}
primary_manifest <- primary_manifest[
  primary_manifest$experiment == "negative_transfer" &
    primary_manifest$config == manifest$config[[1L]] &
    primary_manifest$p == manifest$p[[1L]] &
    primary_manifest$K == manifest$K[[1L]] &
    primary_manifest$rho %in% expected_rhos &
    abs(primary_manifest$cutoff - primary_cutoff) <= 1e-12 &
    primary_manifest$sim_id %in% manifest$sim_id,
  , drop = FALSE
]
primary_manifest <- primary_manifest[
  order(primary_manifest$sim_id, primary_manifest$rho), , drop = FALSE
]
if (nrow(primary_manifest) != expected_replications * length(expected_rhos) ||
    anyDuplicated(primary_manifest[c("sim_id", "rho")])) {
  stop(
    "primary manifest does not contain the exact primary-cutoff prefix.",
    call. = FALSE
  )
}
primary_files <- file.path(
  primary_raw_directory,
  sprintf("task_%06d.csv", as.integer(primary_manifest$task_id))
)
missing_primary_files <- primary_files[!file.exists(primary_files)]
if (length(missing_primary_files) > 0L) {
  stop(
    "primary-cutoff production prefix is incomplete; first missing file: ",
    missing_primary_files[[1L]], call. = FALSE
  )
}
primary_raw <- RoCE:::.read_simulation_result_files(primary_files)
primary_face <- primary_raw[
  primary_raw$method == "one_round_crossfit_ate" &
    primary_raw$estimand_scope == "tate", , drop = FALSE
]
if (nrow(primary_face) != nrow(primary_manifest) ||
    anyDuplicated(primary_face[c("sim_id", "rho")])) {
  stop("production files do not contain one TATE row per seed and rho.",
       call. = FALSE)
}
diagnostic_face_primary <- raw[
  raw$method == "one_round_crossfit_ate" &
    abs(raw$cutoff - primary_cutoff) <= 1e-12,
  , drop = FALSE
]
production_key <- paste(primary_face$sim_id, primary_face$rho, sep = "|")
diagnostic_key <- paste(
  diagnostic_face_primary$sim_id, diagnostic_face_primary$rho, sep = "|"
)
production_index <- match(diagnostic_key, production_key)
if (nrow(diagnostic_face_primary) != nrow(primary_face) ||
    anyNA(production_index) || anyDuplicated(diagnostic_key)) {
  stop("could not pair diagnostic and production primary-cutoff rows.",
       call. = FALSE)
}
production_face_paired <- primary_face[production_index, , drop = FALSE]
production_identity_fields <- c(
  "estimate", "se", "bias", "coverage", "ci_width", "truth", "ci_lower",
  "ci_upper", "target_anchor_weight", "mean_abs_source_weight",
  "max_abs_source_weight", "max_wald_statistic", "mean_wald_statistic",
  "penalized_source_fold_fraction", "max_weight_optimizer_iterations",
  "max_weight_psd_ridge", "weight_psd_ridge_fold_fraction",
  "inference_logit_truncated", "inference_logit_truncation_fraction",
  "inference_max_abs_logit", "inference_safety_clip_count",
  "face_initial_outcome_degenerate", "face_target_only_outcome_degenerate",
  "face_initial_dr_nonconverged", "face_calibrated_dr_nonconverged",
  "face_calibrated_outcome_nonconverged"
)
missing_identity_fields <- setdiff(
  production_identity_fields,
  intersect(names(diagnostic_face_primary), names(production_face_paired))
)
if (length(missing_identity_fields) > 0L) {
  stop(
    "primary-cutoff external identity fields are missing: ",
    paste(missing_identity_fields, collapse = ", "), call. = FALSE
  )
}
production_identity_exact <- vapply(
  production_identity_fields,
  function(field) {
    identical(
      diagnostic_face_primary[[field]], production_face_paired[[field]]
    )
  },
  logical(1L)
)
if (!all(production_identity_exact)) {
  stop(
    "primary-cutoff diagnostic differs from production in: ",
    paste(names(production_identity_exact)[!production_identity_exact],
          collapse = ", "),
    call. = FALSE
  )
}
}

face <- raw[raw$method == "one_round_crossfit_ate", , drop = FALSE]
target <- raw[raw$method == "target_only_ate", , drop = FALSE]
result_key <- function(value) {
  paste(value$task_id, value$rho, value$cutoff, sep = "|")
}
target_index <- match(result_key(face), result_key(target))
if (anyNA(target_index)) {
  stop("could not pair TATE and target-only diagnostic rows.",
       call. = FALSE)
}
face$target_only_bias <- target$bias[target_index]
group_key <- interaction(face$rho, face$cutoff, drop = TRUE, lex.order = TRUE)
summaries <- lapply(split(face, group_key), function(group) {
  empirical_sd <- stats::sd(group$estimate)
  mean_se <- mean(group$se)
  data.frame(
    config = group$config[[1L]],
    p = group$p[[1L]],
    K = group$K[[1L]],
    rho = group$rho[[1L]],
    cutoff = group$cutoff[[1L]],
    n_replications = nrow(group),
    bias = mean(group$bias),
    bias_mcse = stats::sd(group$bias) / sqrt(nrow(group)),
    empirical_sd = empirical_sd,
    mean_se = mean_se,
    se_to_empirical_sd = if (is.finite(empirical_sd) && empirical_sd > 0) {
      mean_se / empirical_sd
    } else {
      NA_real_
    },
    rmse = sqrt(mean(group$bias^2)),
    target_only_rmse = sqrt(mean(group$target_only_bias^2)),
    mean_paired_squared_error_difference_vs_target =
      mean(group$bias^2 - group$target_only_bias^2),
    paired_squared_error_difference_mcse =
      stats::sd(group$bias^2 - group$target_only_bias^2) /
        sqrt(nrow(group)),
    coverage = mean(group$coverage),
    coverage_mcse = sqrt(
      mean(group$coverage) * (1 - mean(group$coverage)) / nrow(group)
    ),
    mean_target_anchor_weight = mean(group$target_anchor_weight),
    mean_source1_weight = mean(group$source1_mean_weight),
    mean_abs_source1_weight = mean(group$source1_mean_abs_weight),
    mean_source1_zero_weight_fold_fraction =
      mean(group$source1_zero_weight_fold_fraction),
    mean_source1_wald_statistic = mean(group$source1_mean_wald_statistic),
    max_source1_wald_statistic = max(group$source1_max_wald_statistic),
    mean_source1_penalty_activation_fraction =
      mean(group$source1_penalty_activation_fraction),
    mean_informative_source_abs_weight =
      mean(group$informative_sources_mean_abs_weight),
    mean_source1_estimate = mean(group$source1_estimate),
    mean_source1_minus_target_estimate =
      mean(group$source1_minus_target_estimate),
    mean_source1_training_discrepancy =
      mean(group$source1_mean_training_discrepancy),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summaries)
summary <- summary[order(summary$rho, summary$cutoff), , drop = FALSE]
rownames(summary) <- NULL

output_directory <- file.path(
  dirname(raw_directory), sprintf("summary_n%03d", expected_replications)
)
roce_write_atomic_directory(
  output_directory,
  writer = function(staging_directory) {
    write.csv(raw, file.path(staging_directory, "replicate_results.csv"),
              row.names = FALSE)
    write.csv(summary, file.path(staging_directory, "cutoff_summary.csv"),
              row.names = FALSE)
    replicate_results_sha256 <- roce_sha256_file(file.path(
      staging_directory, "replicate_results.csv"
    ))
    cutoff_summary_sha256 <- roce_sha256_file(file.path(
      staging_directory, "cutoff_summary.csv"
    ))
    writeLines(c(
      "grouped_cutoff_diagnostic=passed",
      paste0("expected_replications=", expected_replications),
      paste0(
        "simulation_ids=", min(manifest$sim_id), "--", max(manifest$sim_id)
      ),
      paste0("rho_values=", paste(expected_rhos, collapse = ";")),
      paste0("cutoffs=", paste(expected_cutoffs, collapse = ";")),
      paste0("primary_cutoff=", primary_cutoff),
      paste0("package_fingerprint=", provenance$fingerprint),
      paste0("workflow_fingerprint=", workflow_fingerprints[[1L]]),
      paste0("manifest_fingerprint=", expected_manifest_sha256),
      paste0("summary_workflow_fingerprint=", summary_workflow_fingerprint),
      paste0("replicate_results_sha256=", replicate_results_sha256),
      paste0("cutoff_summary_sha256=", cutoff_summary_sha256),
      "primary_cutoff_identity=exact",
      paste0(
        "production_primary_cutoff_identity=",
        if (identical(
          production_identity_mode, "require-production-identity"
        )) {
          "exact"
        } else {
          "not_run_disjoint_pilot"
        }
      ),
      paste0(
        "production_identity_fields=",
        if (length(production_identity_fields) > 0L) {
          paste(production_identity_fields, collapse = ";")
        } else {
          "none"
        }
      ),
      "nuisance_refits_for_cutoff_grid=0"
    ), file.path(staging_directory, "audit_passed.txt"))
  },
  caller = "grouped cutoff diagnostic summary"
)
message("[pass] wrote ", output_directory)
print(summary, row.names = FALSE)
