#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
output_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/rhc_supported_primary"
}
result_path <- file.path(output_root, "rhc_direct_tate.rds")
if (!file.exists(result_path)) {
  stop("RHC result not found: ", result_path, call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "atomic_output.R"))

result <- readRDS(result_path)
fit <- result$tate_fit
methods <- result$methods
pairwise <- result$pairwise
tolerance <- 1e-10
expected_primary_cutoff <- 1 / get(
  "AGG_WALD_LAMBDA", envir = asNamespace("RoCE")
)

checks <- list()
add_check <- function(name, passed, detail) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
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
add_check("tested_package_provenance", provenance$passed, provenance$detail)
expected_rhc_workflow_fingerprint <- roce_sha256_environment(
  "ROCE_EXPECT_RHC_WORKFLOW_FINGERPRINT"
)
expected_rhc_data_fingerprint <- roce_sha256_environment(
  "ROCE_EXPECT_RHC_DATA_FINGERPRINT"
)
add_check(
  "rhc_workflow_provenance",
  identical(
    tolower(metadata_scalar("rhc_workflow_fingerprint")),
    expected_rhc_workflow_fingerprint
  ),
  metadata_scalar("rhc_workflow_fingerprint")
)
add_check(
  "rhc_data_provenance",
  identical(
    tolower(metadata_scalar("rhc_data_fingerprint")),
    expected_rhc_data_fingerprint
  ),
  metadata_scalar("rhc_data_fingerprint")
)
scheduler_fields <- c(
  "slurm_job_id", "slurm_array_job_id", "slurm_array_task_id",
  "slurm_cluster_name", "slurm_job_partition", "slurm_node_list",
  "compute_hostname"
)
scheduler_values <- vapply(scheduler_fields, metadata_scalar, character(1L))
add_check(
  "scheduler_provenance",
  all(!is.na(scheduler_values) & nzchar(scheduler_values)),
  paste(paste0(scheduler_fields, "=", scheduler_values), collapse = ";")
)
add_check(
  "target_anchor_nuisance_grid",
  identical(as.integer(result$metadata$target_anchor_nlambda), 100L),
  metadata_scalar("target_anchor_nlambda")
)

expected_nuisance_lambda_rule <- trimws(Sys.getenv(
  "ROCE_EXPECT_NUISANCE_LAMBDA_RULE", "min"
))
if (!expected_nuisance_lambda_rule %in% c("min", "1se")) {
  stop(
    "ROCE_EXPECT_NUISANCE_LAMBDA_RULE must be 'min' or '1se'.",
    call. = FALSE
  )
}
expected_nlambda <- suppressWarnings(as.integer(Sys.getenv(
  "ROCE_EXPECT_NLAMBDA_INIT", "100"
)))
if (length(expected_nlambda) != 1L || is.na(expected_nlambda) ||
    expected_nlambda < 2L) {
  stop("ROCE_EXPECT_NLAMBDA_INIT must be one integer >= 2.", call. = FALSE)
}
add_check(
  "rhc_analysis_configuration",
  identical(as.integer(result$metadata$n_folds), 10L) &&
    identical(as.integer(result$metadata$nlambda_init), expected_nlambda) &&
    identical(as.character(result$metadata$nuisance_lambda_rule),
              expected_nuisance_lambda_rule) &&
    isTRUE(all.equal(as.numeric(result$metadata$M_tau), 5,
                     tolerance = 0)) &&
    isTRUE(all.equal(as.numeric(result$metadata$M_tau_inference), 5,
                     tolerance = 0)) &&
    isTRUE(all.equal(
      as.numeric(result$metadata$aggregation_lambda),
      1 / expected_primary_cutoff,
      tolerance = 1e-12
    )) &&
    isTRUE(all.equal(
      as.numeric(result$metadata$aggregation_cutoff),
      expected_primary_cutoff,
      tolerance = 1e-12
    )),
  sprintf(
    paste0(
      "folds=%s;nlambda=%s;rule=%s;M_fit=%s;M_inf=%s;",
      "lambda=%s;cutoff=%s"
    ),
    metadata_scalar("n_folds"), metadata_scalar("nlambda_init"),
    metadata_scalar("nuisance_lambda_rule"), metadata_scalar("M_tau"),
    metadata_scalar("M_tau_inference"),
    metadata_scalar("aggregation_lambda"),
    metadata_scalar("aggregation_cutoff")
  )
)

required_method_columns <- c(
  "method", "estimate", "se", "ci_lower", "ci_upper"
)
add_check(
  "method_schema",
  all(required_method_columns %in% names(methods)),
  paste(setdiff(required_method_columns, names(methods)), collapse = ",")
)
if (all(required_method_columns %in% names(methods))) {
  finite_methods <- all(is.finite(unlist(
    methods[c("estimate", "se", "ci_lower", "ci_upper")]
  )))
  add_check("method_values_finite", finite_methods, nrow(methods))
  add_check("method_standard_errors_positive", all(methods$se > 0),
            paste(signif(methods$se, 6), collapse = ","))
  expected_lower <- methods$estimate - stats::qnorm(0.975) * methods$se
  expected_upper <- methods$estimate + stats::qnorm(0.975) * methods$se
  ci_error <- max(abs(c(
    methods$ci_lower - expected_lower,
    methods$ci_upper - expected_upper
  )))
  add_check("method_ci_arithmetic", ci_error <= tolerance,
            sprintf("max_abs_error=%.3g", ci_error))
  method_labels <- as.character(methods$method)
  face_index <- match("RoCE", method_labels)
  target_index <- match("Target-only", method_labels)
  method_indices_valid <- !anyNA(c(face_index, target_index))
  method_fit_error <- if (method_indices_valid) {
    max(abs(c(
      methods$estimate[face_index] - fit$estimate,
      methods$se[face_index] - fit$se,
      methods$estimate[target_index] - fit$target_only$estimate,
      methods$se[target_index] - fit$target_only$se
    )))
  } else {
    Inf
  }
  add_check(
    "method_table_matches_core_fit",
    method_indices_valid && is.finite(method_fit_error) &&
      method_fit_error <= tolerance,
    sprintf("max_abs_error=%.3g", method_fit_error)
  )
}

site_sizes <- as.integer(result$metadata$site_n)
site_labels <- c(
  unname(result$metadata$target_site),
  unname(result$metadata$source_sites)
)
site_components <- result$metadata$site_components
component_names_expected <- c(
  "target", paste0("s", seq_len(result$metadata$K))
)
components_valid <- is.character(site_components) &&
  identical(names(site_components), component_names_expected) &&
  length(site_components) == length(site_labels) &&
  all(nzchar(site_components))
raw_components <- if (components_valid) {
  unlist(strsplit(unname(site_components), ";", fixed = TRUE), use.names = FALSE)
} else {
  character(0)
}
excluded_components <- if (
    is.character(result$metadata$excluded_site_levels) &&
    length(result$metadata$excluded_site_levels) == 1L &&
    !identical(result$metadata$excluded_site_levels, "none")) {
  strsplit(
    result$metadata$excluded_site_levels, ";", fixed = TRUE
  )[[1L]]
} else {
  character(0)
}
add_check(
  "site_component_metadata",
  components_valid && !anyDuplicated(raw_components) &&
    !any(excluded_components %in% raw_components),
  if (components_valid) {
    paste(paste(site_labels, site_components, sep = "<-"), collapse = ",")
  } else {
    "missing_or_malformed"
  }
)
expected_exclusion_text <- trimws(Sys.getenv(
  "ROCE_EXPECT_RHC_EXCLUDED_SITES", ""
))
expected_exclusions <- if (nzchar(expected_exclusion_text)) {
  unique(trimws(strsplit(
    expected_exclusion_text, ",", fixed = TRUE
  )[[1L]]))
} else {
  character(0)
}
add_check(
  "requested_support_screen",
  setequal(excluded_components, expected_exclusions),
  sprintf(
    "expected={%s};observed={%s}",
    paste(expected_exclusions, collapse = ","),
    paste(excluded_components, collapse = ",")
  )
)
direct_phi_valid <- is.numeric(fit$all_phi_tau) &&
  length(fit$all_phi_tau) == sum(site_sizes) &&
  all(is.finite(fit$all_phi_tau))
add_check(
  "direct_tate_pseudovalue_schema",
  direct_phi_valid,
  sprintf(
    "type=%s;length=%d;expected=%d;nonfinite=%d",
    typeof(fit$all_phi_tau), length(fit$all_phi_tau), sum(site_sizes),
    if (is.numeric(fit$all_phi_tau)) sum(!is.finite(fit$all_phi_tau)) else NA_integer_
  )
)
expected_variance <- if (direct_phi_valid) {
  RoCE:::.multisite_pseudovalue_variance(fit$all_phi_tau, site_sizes)
} else {
  NA_real_
}
variance_error <- abs(fit$variance - expected_variance)
estimate_error <- if (direct_phi_valid) {
  abs(fit$estimate - mean(fit$all_phi_tau))
} else {
  NA_real_
}
add_check(
  "direct_tate_pseudovalue_estimate",
  estimate_error <= tolerance,
  sprintf("abs_error=%.3g", estimate_error)
)
add_check(
  "direct_tate_site_centered_variance",
  variance_error <= tolerance,
  sprintf("abs_error=%.3g", variance_error)
)
add_check(
  "direct_tate_se_variance_identity",
  abs(fit$se^2 - fit$variance) <= tolerance,
  sprintf("abs_error=%.3g", abs(fit$se^2 - fit$variance))
)

reconstruct_target_phi <- function(fold_info) {
  if (is.null(fold_info) || length(fold_info) == 0L) return(NULL)
  indices <- unlist(lapply(fold_info, `[[`, "target_idx"), use.names = FALSE)
  if (!is.numeric(indices) || length(indices) == 0L ||
      any(!is.finite(indices)) || any(indices < 1L) ||
      any(indices != floor(indices))) {
    return(NULL)
  }
  target_phi <- rep(NA_real_, max(indices))
  assignment_count <- integer(length(target_phi))
  for (info in fold_info) {
    fold_indices <- as.integer(info$target_idx)
    fold_values <- as.numeric(info$varphi_ot) + info$fold_target_estimate
    if (length(fold_indices) != length(fold_values)) return(NULL)
    target_phi[fold_indices] <- fold_values
    assignment_count[fold_indices] <- assignment_count[fold_indices] + 1L
  }
  if (any(assignment_count != 1L)) return(NULL)
  target_phi
}

target_phi_stored <- fit$intermediates$target_only_phi
target_phi_from_arms <- if (
    is.numeric(fit$arm_results$mu1$intermediates$target_only_phi) &&
    is.numeric(fit$arm_results$mu0$intermediates$target_only_phi) &&
    length(fit$arm_results$mu1$intermediates$target_only_phi) ==
      length(fit$arm_results$mu0$intermediates$target_only_phi)) {
  fit$arm_results$mu1$intermediates$target_only_phi -
    fit$arm_results$mu0$intermediates$target_only_phi
} else {
  NULL
}
target_phi_reconstructed <- reconstruct_target_phi(fit$intermediates$fold_info)
target_phi_source <- if (is.numeric(target_phi_stored) &&
                         length(target_phi_stored) > 0L) {
  "stored"
} else if (is.numeric(target_phi_reconstructed) &&
           length(target_phi_reconstructed) > 0L) {
  "direct_fold_info"
} else if (is.numeric(target_phi_from_arms) &&
           length(target_phi_from_arms) > 0L) {
  "arm_difference"
} else {
  "unavailable"
}
target_phi <- if (identical(target_phi_source, "stored")) {
  target_phi_stored
} else if (identical(target_phi_source, "direct_fold_info")) {
  target_phi_reconstructed
} else if (identical(target_phi_source, "arm_difference")) {
  target_phi_from_arms
} else {
  NULL
}
target_phi_valid <- is.numeric(target_phi) &&
  length(target_phi) == site_sizes[[1L]] && all(is.finite(target_phi))
add_check(
  "target_only_pseudovalue_schema",
  target_phi_valid,
  sprintf(
    "source=%s;length=%d;expected=%d;nonfinite=%d",
    target_phi_source, length(target_phi), site_sizes[[1L]],
    if (is.numeric(target_phi)) sum(!is.finite(target_phi)) else NA_integer_
  )
)
target_variance <- if (target_phi_valid) {
  RoCE:::.multisite_pseudovalue_variance(target_phi, length(target_phi))
} else {
  NA_real_
}
add_check(
  "target_only_pseudovalue_estimate",
  target_phi_valid &&
    abs(fit$target_only$estimate - mean(target_phi)) <= tolerance,
  sprintf(
    "abs_error=%.3g",
    abs(fit$target_only$estimate - mean(target_phi))
  )
)
add_check(
  "target_only_site_centered_variance",
  target_phi_valid &&
    abs(fit$target_only$variance - target_variance) <= tolerance,
  sprintf(
    "abs_error=%.3g",
    abs(fit$target_only$variance - target_variance)
  )
)

weights <- as.numeric(fit$weights)
anchor <- 1 - sum(weights)
add_check(
  "common_source_weights_finite",
  length(weights) == result$metadata$K && all(is.finite(weights)),
  paste(signif(weights, 6), collapse = ",")
)
add_check(
  "target_anchor_finite",
  is.finite(anchor),
  sprintf("anchor=%.6f", anchor)
)
fold_weights <- as.matrix(fit$fold_weights)
valid_fold_weight_shape <- identical(
  dim(fold_weights), c(as.integer(fit$n_folds), as.integer(result$metadata$K))
)
add_check(
  "fold_weight_shape_and_finiteness",
  valid_fold_weight_shape && all(is.finite(fold_weights)),
  paste(dim(fold_weights), collapse = "x")
)
if (valid_fold_weight_shape) {
  average_weight_error <- max(abs(weights - colMeans(fold_weights)))
  add_check(
    "reported_weights_equal_fold_average",
    average_weight_error <= tolerance,
    sprintf("max_abs_error=%.3g", average_weight_error)
  )
}
wald <- as.matrix(fit$fold_wald_statistics)
penalty <- as.matrix(fit$fold_penalty_coefficients)
lambda <- as.numeric(fit$fold_lambdas)
valid_fold_diagnostic_shape <- identical(dim(wald), dim(fold_weights)) &&
  identical(dim(penalty), dim(fold_weights)) &&
  length(lambda) == nrow(fold_weights)
add_check(
  "fold_wald_penalty_shape_and_finiteness",
  valid_fold_diagnostic_shape && all(is.finite(c(wald, penalty, lambda))),
  sprintf(
    "wald=%s;penalty=%s;lambda=%d",
    paste(dim(wald), collapse = "x"),
    paste(dim(penalty), collapse = "x"),
    length(lambda)
  )
)
if (valid_fold_diagnostic_shape) {
  expected_penalty <- pmax(lambda * wald - 1, 0)
  penalty_error <- max(abs(penalty - expected_penalty))
  add_check(
    "fold_truncated_wald_penalty_identity",
    penalty_error <= tolerance,
    sprintf("max_abs_error=%.3g", penalty_error)
  )
}

required_pairwise_columns <- c(
  "source", "estimate", "weight", "mean_wald_statistic",
  "max_wald_statistic", "penalty_activation_fraction"
)
pairwise_ok <- all(required_pairwise_columns %in% names(pairwise))
add_check(
  "source_diagnostic_schema",
  pairwise_ok,
  paste(setdiff(required_pairwise_columns, names(pairwise)), collapse = ",")
)
if (pairwise_ok) {
  pairwise_numeric <- unlist(pairwise[setdiff(
    required_pairwise_columns, "source"
  )])
  add_check(
    "source_diagnostics_finite",
    all(is.finite(pairwise_numeric)),
    nrow(pairwise)
  )
  valid_ranges <- all(pairwise$mean_wald_statistic >= 0) &&
    all(pairwise$max_wald_statistic >= pairwise$mean_wald_statistic) &&
    all(pairwise$penalty_activation_fraction >= 0) &&
    all(pairwise$penalty_activation_fraction <= 1)
  add_check("source_diagnostic_ranges", valid_ranges,
            "Wald>=0; activation in [0,1]")
  pairwise_weight_error <- max(abs(pairwise$weight - weights))
  add_check(
    "source_table_matches_reported_weights",
    pairwise_weight_error <= tolerance,
    sprintf("max_abs_error=%.3g", pairwise_weight_error)
  )
}

comparison_rows <- lapply(names(result$comparison_results), function(method) {
  comparison <- result$comparison_results[[method]]
  component <- comparison$components
  dr_diagnostics <- RoCE:::.comparison_variance_diagnostics(comparison)
  covariance_identity_error <- abs(
    component$cross_arm_covariance -
      (component$mu1_influence_se^2 + component$mu0_influence_se^2 -
         component$se_analytic^2) / 2
  )
  data.frame(
    method = method,
    covariance_identity_error = covariance_identity_error,
    correlation = component$cross_arm_correlation,
    dr_weight_n = dr_diagnostics$dr_weight_n,
    dr_weight_n_clipped = dr_diagnostics$dr_weight_n_clipped,
    dr_weight_fraction_clipped =
      dr_diagnostics$dr_weight_fraction_clipped,
    dr_weight_max_site_fraction_clipped =
      dr_diagnostics$dr_weight_max_site_fraction_clipped,
    dr_weight_min_before_clipping =
      dr_diagnostics$dr_weight_min_before_clipping,
    dr_weight_max_before_clipping =
      dr_diagnostics$dr_weight_max_before_clipping,
    stringsAsFactors = FALSE
  )
})
comparison_audit <- do.call(rbind, comparison_rows)
if (is.null(comparison_audit)) {
  comparison_audit <- data.frame(
    method = character(0),
    covariance_identity_error = numeric(0),
    correlation = numeric(0),
    dr_weight_n = integer(0),
    dr_weight_n_clipped = integer(0),
    dr_weight_fraction_clipped = numeric(0),
    dr_weight_max_site_fraction_clipped = numeric(0),
    dr_weight_min_before_clipping = numeric(0),
    dr_weight_max_before_clipping = numeric(0)
  )
}
add_check(
  "comparison_cross_arm_covariance_identity",
  nrow(comparison_audit) == 0L ||
    all(comparison_audit$covariance_identity_error <= tolerance),
  if (nrow(comparison_audit) == 0L) "no comparisons" else sprintf(
    "max_abs_error=%.3g", max(comparison_audit$covariance_identity_error)
  )
)
add_check(
  "comparison_cross_arm_correlation_range",
  nrow(comparison_audit) == 0L ||
    all(is.finite(comparison_audit$correlation)) &&
    all(abs(comparison_audit$correlation) <= 1 + 1e-8),
  if (nrow(comparison_audit) == 0L) "no comparisons" else paste(
    signif(comparison_audit$correlation, 6), collapse = ","
  )
)
dr_comparisons <- comparison_audit$method %in% c(
  "federated_dr", "pooled_dr"
)
dr_rows <- comparison_audit[dr_comparisons, , drop = FALSE]
dr_diagnostic_values <- unlist(dr_rows[c(
  "dr_weight_n", "dr_weight_n_clipped",
  "dr_weight_fraction_clipped", "dr_weight_max_site_fraction_clipped",
  "dr_weight_min_before_clipping", "dr_weight_max_before_clipping"
)])
dr_diagnostics_valid <- nrow(dr_rows) == 2L &&
  all(is.finite(dr_diagnostic_values)) &&
  all(dr_rows$dr_weight_n >= 1L) &&
  all(dr_rows$dr_weight_n == floor(dr_rows$dr_weight_n)) &&
  all(dr_rows$dr_weight_n_clipped >= 0L) &&
  all(dr_rows$dr_weight_n_clipped <= dr_rows$dr_weight_n) &&
  all(dr_rows$dr_weight_n_clipped == floor(dr_rows$dr_weight_n_clipped)) &&
  all(dr_rows$dr_weight_fraction_clipped >= 0) &&
  all(dr_rows$dr_weight_fraction_clipped <= 1) &&
  all(dr_rows$dr_weight_max_site_fraction_clipped >=
        dr_rows$dr_weight_fraction_clipped) &&
  all(dr_rows$dr_weight_max_site_fraction_clipped <= 1) &&
  all(dr_rows$dr_weight_min_before_clipping > 0) &&
  all(dr_rows$dr_weight_max_before_clipping >=
        dr_rows$dr_weight_min_before_clipping) &&
  all(abs(
    dr_rows$dr_weight_fraction_clipped -
      dr_rows$dr_weight_n_clipped / dr_rows$dr_weight_n
  ) <= tolerance)
add_check(
  "comparison_density_ratio_clipping_diagnostics",
  dr_diagnostics_valid,
  if (nrow(dr_rows) == 0L) "no DR comparisons" else paste(
    sprintf(
      "%s:%d/%d",
      dr_rows$method, dr_rows$dr_weight_n_clipped, dr_rows$dr_weight_n
    ),
    collapse = ";"
  )
)

expected_n_bootstrap <- suppressWarnings(as.integer(Sys.getenv(
  "ROCE_EXPECT_N_BOOTSTRAP", "5000"
)))
if (length(expected_n_bootstrap) != 1L || is.na(expected_n_bootstrap) ||
    expected_n_bootstrap < 2L) {
  stop("ROCE_EXPECT_N_BOOTSTRAP must be one integer >= 2.", call. = FALSE)
}
bootstrap_checks <- vapply(
  result$comparison_results,
  function(comparison) {
    component <- comparison$components
    identical(component$variance_method, "bootstrap") &&
      is.finite(component$se_bootstrap) && component$se_bootstrap > 0 &&
      identical(as.integer(component$n_bootstrap), expected_n_bootstrap) &&
      abs(comparison$se - component$se_bootstrap) <= tolerance
  },
  logical(1L)
)
metadata_bootstrap_ok <-
  identical(result$metadata$variance_method, "bootstrap") &&
  identical(as.integer(result$metadata$n_bootstrap), expected_n_bootstrap)
add_check(
  "comparison_bootstrap_configuration",
  length(bootstrap_checks) > 0L && all(bootstrap_checks) &&
    metadata_bootstrap_ok,
  sprintf(
    "expected=%d;methods=%d;valid=%d;metadata=%s",
    expected_n_bootstrap, length(bootstrap_checks), sum(bootstrap_checks),
    metadata_bootstrap_ok
  )
)

nuisance <- fit$nuisance_fit_diagnostics
nonconvergence_columns <- grep(
  "_nonconverged$", names(nuisance), value = TRUE
)
nonconverged_total <- if (length(nonconvergence_columns) == 0L) {
  NA_real_
} else {
  sum(unlist(nuisance[nonconvergence_columns]), na.rm = TRUE)
}
add_check(
  "nuisance_convergence",
  is.finite(nonconverged_total) && nonconverged_total == 0,
  sprintf("nonconverged_fits=%s", nonconverged_total)
)
cv_exclusion_columns <- grep(
  "_cv_invalid_(fold_fits|lambdas)$", names(nuisance), value = TRUE
)
cv_exclusion_values <- if (length(cv_exclusion_columns) > 0L) {
  suppressWarnings(as.numeric(unlist(nuisance[cv_exclusion_columns])))
} else {
  numeric(0)
}
add_check(
  "nuisance_cv_candidate_exclusion_diagnostics",
  length(cv_exclusion_values) > 0L &&
    all(is.finite(cv_exclusion_values)) && all(cv_exclusion_values >= 0),
  sprintf(
    "columns=%d;total_exclusions=%s",
    length(cv_exclusion_columns),
    if (length(cv_exclusion_values) > 0L) {
      format(sum(cv_exclusion_values), scientific = FALSE)
    } else {
      "unavailable"
    }
  )
)

density_ratio_fits <- list()
outcome_fits <- list()
for (arm_name in intersect(c("mu1", "mu0"), names(fit$arm_results))) {
  arm_result <- fit$arm_results[[arm_name]]
  for (outer_index in seq_along(arm_result$fold_results)) {
    source_results <- arm_result$fold_results[[outer_index]]$source_results
    for (source_name in names(source_results)) {
      source_result <- source_results[[source_name]]
      prefix <- paste(arm_name, outer_index, source_name, sep = "/")
      for (inner_name in names(source_result$per_k2_gamma)) {
        density_ratio_fits[[paste(prefix, "initial", inner_name, sep = "/")]] <-
          source_result$per_k2_gamma[[inner_name]]
      }
      density_ratio_fits[[paste(prefix, "calibrated", sep = "/")]] <-
        source_result$gamma_s
      outcome_fits[[paste(prefix, "calibrated", sep = "/")]] <-
        source_result$alpha_ts
    }
  }
}
fit_max_abs <- function(fits) vapply(
  fits,
  function(value) {
    value <- as.numeric(value)
    if (length(value) > 0L && all(is.finite(value))) max(abs(value)) else Inf
  },
  numeric(1L)
)
fit_line_search_failures <- function(fits) vapply(
  fits,
  function(value) {
    diagnostic <- suppressWarnings(as.integer(
      attr(value, "line_search_failures")
    ))
    if (length(diagnostic) == 1L && !is.na(diagnostic)) diagnostic else NA_integer_
  },
  integer(1L)
)
fit_update_ratios <- function(fits) vapply(
  fits,
  function(value) {
    diagnostic <- suppressWarnings(as.numeric(
      attr(value, "update_to_threshold_ratio")
    ))
    if (length(diagnostic) == 1L && is.finite(diagnostic) &&
        diagnostic >= 0) diagnostic else NA_real_
  },
  numeric(1L)
)
parameter_bound <- get("PARAM_MAX", envir = asNamespace("RoCE"))
add_parameter_fit_checks <- function(label, fits) {
  max_abs <- fit_max_abs(fits)
  line_search_failures <- fit_line_search_failures(fits)
  update_ratios <- fit_update_ratios(fits)
  add_check(
    paste0(label, "_fit_schema_and_finiteness"),
    length(fits) > 0L && all(is.finite(max_abs)) &&
      all(!is.na(line_search_failures)) && all(is.finite(update_ratios)),
    sprintf("fits=%d", length(fits))
  )
  add_check(
    paste0(label, "_coefficients_away_from_hard_bound"),
    length(max_abs) > 0L && max(max_abs) < 0.9 * parameter_bound,
    sprintf(
      "max_abs=%.6g;hard_bound=%.6g",
      max(max_abs), parameter_bound
    )
  )
  add_check(
    paste0(label, "_line_search"),
    length(line_search_failures) > 0L &&
      all(line_search_failures == 0L),
    sprintf("failures=%s", sum(line_search_failures, na.rm = TRUE))
  )
  invisible(list(
    max_abs = max_abs,
    line_search_failures = line_search_failures,
    update_ratios = update_ratios
  ))
}

add_parameter_fit_checks("density_ratio", density_ratio_fits)
add_parameter_fit_checks("outcome", outcome_fits)
support_floors <- vapply(
  density_ratio_fits,
  function(value) {
    diagnostic <- suppressWarnings(as.numeric(
      attr(value, "support_penalty_floor")
    ))
    if (length(diagnostic) == 1L && is.finite(diagnostic)) diagnostic else NA_real_
  },
  numeric(1L)
)
support_floor_applied <- vapply(
  density_ratio_fits,
  function(value) isTRUE(attr(value, "support_penalty_floor_applied")),
  logical(1L)
)
add_check(
  "density_ratio_support_floor_metadata",
  length(support_floors) > 0L && all(is.finite(support_floors)),
  sprintf(
    "applied_fits=%d;max_floor=%.6g",
    sum(support_floor_applied), max(support_floors, na.rm = TRUE)
  )
)

clip <- fit$clip_diagnostics$total
clip_values <- unlist(clip[c(
  "n_obs", "logit_truncated", "logit_truncation_fraction",
  "weight_min_clipped", "weight_max_clipped", "ratio_min_clipped",
  "ratio_max_clipped", "max_abs_logit", "max_raw_weight"
)])
add_check(
  "clipping_diagnostics_finite",
  length(clip_values) > 0L && all(is.finite(clip_values)),
  paste(signif(clip_values, 6), collapse = ",")
)
clip_fraction_error <- abs(
  clip$logit_truncation_fraction - clip$logit_truncated / clip$n_obs
)
add_check(
  "clipping_fraction_identity",
  is.finite(clip_fraction_error) && clip$n_obs > 0L &&
    clip_fraction_error <= tolerance,
  sprintf("abs_error=%.3g", clip_fraction_error)
)

audit <- do.call(rbind, checks)
failed <- audit[!audit$passed, , drop = FALSE]
job_label <- gsub(
  "[^A-Za-z0-9_-]", "_", Sys.getenv("SLURM_JOB_ID", "local")
)
audit_label <- sprintf(
  "%s_job%s", format(Sys.time(), "%Y%m%d_%H%M%S"), job_label
)
audit_directory <- file.path(output_root, "audits", audit_label)
roce_write_atomic_directory(
  audit_directory,
  writer = function(staging_directory) {
    write.csv(
      audit, file.path(staging_directory, "rhc_direct_tate_audit.csv"),
      row.names = FALSE
    )
    write.csv(
      comparison_audit,
      file.path(staging_directory, "rhc_direct_tate_comparison_audit.csv"),
      row.names = FALSE
    )
    gate_name <- if (nrow(failed) == 0L) {
      "rhc_direct_tate_audit_passed.txt"
    } else {
      "rhc_direct_tate_audit_failed.txt"
    }
    writeLines(c(
      paste0(
        "rhc_direct_tate_audit=",
        if (nrow(failed) == 0L) "passed" else "failed"
      ),
      paste0("failed_checks=", nrow(failed)),
      paste0("package_fingerprint=", provenance$fingerprint),
      paste0(
        "rhc_workflow_fingerprint=", expected_rhc_workflow_fingerprint
      ),
      paste0("rhc_data_fingerprint=", expected_rhc_data_fingerprint),
      paste0("rhc_result_md5=", unname(tools::md5sum(result_path))),
      paste0(
        "audit_driver_md5=",
        unname(tools::md5sum(file.path(
          "scripts", "slurm", "audit_rhc_direct_tate.R"
        )))
      )
    ), file.path(staging_directory, gate_name))
  },
  caller = "RHC TATE audit"
)
print(audit, row.names = FALSE)
if (nrow(failed) > 0L) {
  stop(
    "RHC audit failed: ",
    paste(failed$check, collapse = ", "), "; inspect ", audit_directory,
    call. = FALSE
  )
}
message("[done] RHC TATE audit passed: ", audit_directory)
}

main()
