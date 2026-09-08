#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop(
    paste0(
      "usage: summarize_rhc_cutoff_sensitivity.R RHC_RESULT.rds ",
      "[OUTPUT_ROOT] [CUTOFFS]"
    ),
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_root <- if (length(args) >= 2L) args[[2L]] else dirname(input_path)
cutoff_text <- if (length(args) >= 3L) args[[3L]] else "1,1.5,2,2.5,3"
cutoffs <- suppressWarnings(as.numeric(
  trimws(strsplit(cutoff_text, ",", fixed = TRUE)[[1L]])
))
if (length(cutoffs) < 1L || any(!is.finite(cutoffs)) || any(cutoffs <= 0) ||
    anyDuplicated(cutoffs)) {
  stop("CUTOFFS must be unique, positive, comma-separated numbers.",
       call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
derivation_workflow_fingerprint <- roce_sha256_environment(
  "ROCE_RHC_SENSITIVITY_WORKFLOW_FINGERPRINT"
)

result <- readRDS(input_path)
required_result_fields <- c("tate_fit", "metadata", "weights")
missing_result_fields <- setdiff(required_result_fields, names(result))
if (length(missing_result_fields) > 0L) {
  stop(
    "RHC result is missing field(s): ",
    paste(missing_result_fields, collapse = ", "),
    call. = FALSE
  )
}
metadata_scalar <- function(name) {
  value <- result$metadata[[name]]
  if (is.null(value) || length(value) != 1L) NA_character_ else as.character(value)
}
provenance <- roce_result_provenance(
  data.frame(
    package_library = metadata_scalar("package_library"),
    package_fingerprint = metadata_scalar("package_fingerprint"),
    stringsAsFactors = FALSE
  ),
  project_library
)
if (!isTRUE(provenance$passed)) {
  stop(
    "RHC result failed tested-package provenance validation: ",
    provenance$detail,
    call. = FALSE
  )
}
primary_workflow_fingerprint <- tolower(metadata_scalar(
  "rhc_workflow_fingerprint"
))
primary_data_fingerprint <- tolower(metadata_scalar("rhc_data_fingerprint"))
if (!grepl("^[0-9a-f]{64}$", primary_workflow_fingerprint) ||
    !grepl("^[0-9a-f]{64}$", primary_data_fingerprint)) {
  stop("RHC result has invalid workflow or data provenance.", call. = FALSE)
}
arm_results <- result$tate_fit$arm_results
if (!is.list(arm_results) || !all(c("mu1", "mu0") %in% names(arm_results))) {
  stop("RHC result does not contain both arm-specific cross-fit results.",
       call. = FALSE)
}

site_sizes <- result$metadata$site_n
if (is.null(names(site_sizes)) || !identical(names(site_sizes)[[1L]], "t") ||
    any(!is.finite(site_sizes)) || any(site_sizes <= 0)) {
  stop("RHC metadata must contain positive named site sizes beginning with t.",
       call. = FALSE)
}
# The aggregation helper uses data_split only for ordered site names and sample
# sizes. This minimal reconstruction keeps the cutoff sensitivity independent
# of individual-level RHC records and reuses the exact fitted nuisances and
# held-out pseudo-values.
data_split_sizes <- lapply(as.integer(site_sizes), function(n) list(n = n))
names(data_split_sizes) <- names(site_sizes)

fits <- lapply(cutoffs, function(cutoff) {
  calculate_tate_crossfit_aggregation(
    data_split = data_split_sizes,
    mu1_result = arm_results$mu1,
    mu0_result = arm_results$mu0,
    lambda_selection = 1 / cutoff,
    lambda_rule = "min",
    verbose = FALSE
  )
})
names(fits) <- format(cutoffs, scientific = FALSE, trim = TRUE)

reference_cutoff <- as.numeric(result$metadata$aggregation_cutoff)
reference_index <- which(abs(cutoffs - reference_cutoff) <= 1e-12)
if (length(reference_index) != 1L) {
  stop("the cutoff grid must contain the primary RHC cutoff exactly once.",
       call. = FALSE)
} else {
  reference_fit <- fits[[reference_index]]
  reference_checks <- c(
    estimate = abs(reference_fit$estimate - result$tate_fit$estimate),
    se = abs(reference_fit$se - result$tate_fit$se),
    weights = max(abs(reference_fit$weights - result$weights))
  )
  if (any(!is.finite(reference_checks)) || any(reference_checks > 1e-10)) {
    stop(
      sprintf(
        paste0(
          "reference cutoff did not reproduce the primary fit; discrepancies: ",
          "estimate=%.3g, se=%.3g, weights=%.3g."
        ),
        reference_checks[["estimate"]], reference_checks[["se"]],
        reference_checks[["weights"]]
      ),
      call. = FALSE
    )
  }
}

summary <- do.call(rbind, lapply(seq_along(fits), function(index) {
  fit <- fits[[index]]
  data.frame(
    package_library = provenance$library,
    package_fingerprint = provenance$fingerprint,
    primary_rhc_workflow_fingerprint = primary_workflow_fingerprint,
    primary_rhc_data_fingerprint = primary_data_fingerprint,
    derivation_workflow_fingerprint = derivation_workflow_fingerprint,
    primary_cutoff = reference_cutoff,
    primary_identity_checked =
      abs(cutoffs[[index]] - reference_cutoff) <= 1e-12,
    primary_identity_exact = if (
        abs(cutoffs[[index]] - reference_cutoff) <= 1e-12) TRUE else NA,
    cutoff = cutoffs[[index]],
    lambda = 1 / cutoffs[[index]],
    estimate = fit$estimate,
    se = fit$se,
    ci_lower = fit$ci_lower,
    ci_upper = fit$ci_upper,
    target_anchor_weight = 1 - sum(fit$weights),
    max_abs_source_weight = max(abs(fit$weights)),
    mean_penalty_activation_fraction = mean(
      fit$fold_penalty_coefficients > 0
    ),
    stringsAsFactors = FALSE
  )
}))
weights <- do.call(rbind, lapply(seq_along(fits), function(index) {
  fit <- fits[[index]]
  data.frame(
    package_library = provenance$library,
    package_fingerprint = provenance$fingerprint,
    primary_rhc_workflow_fingerprint = primary_workflow_fingerprint,
    primary_rhc_data_fingerprint = primary_data_fingerprint,
    derivation_workflow_fingerprint = derivation_workflow_fingerprint,
    primary_cutoff = reference_cutoff,
    primary_identity_checked =
      abs(cutoffs[[index]] - reference_cutoff) <= 1e-12,
    primary_identity_exact = if (
        abs(cutoffs[[index]] - reference_cutoff) <= 1e-12) TRUE else NA,
    cutoff = cutoffs[[index]],
    lambda = 1 / cutoffs[[index]],
    source = names(fit$weights),
    weight = as.numeric(fit$weights),
    mean_wald_statistic = colMeans(fit$fold_wald_statistics),
    max_wald_statistic = apply(fit$fold_wald_statistics, 2L, max),
    penalty_activation_fraction = colMeans(
      fit$fold_penalty_coefficients > 0
    ),
    stringsAsFactors = FALSE
  )
}))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
summary_path <- file.path(output_root, "rhc_direct_tate_cutoff_sensitivity.csv")
weights_path <- file.path(
  output_root, "rhc_direct_tate_cutoff_sensitivity_weights.csv"
)
all_paths <- c(summary_path, weights_path)
existing_paths <- all_paths[file.exists(all_paths)]
if (length(existing_paths) > 0L) {
  stop(
    "refusing to overwrite existing cutoff-sensitivity output(s): ",
    paste(existing_paths, collapse = ", "),
    call. = FALSE
  )
}

roce_commit_csv_bundle(
  values = list(summary, weights),
  output_paths = c(summary_path, weights_path)
)

print(summary, row.names = FALSE, digits = 5)
message("[done] wrote ", summary_path, " and ", weights_path)
