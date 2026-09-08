#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/smoke_final_direct_tate_raw"
}
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
if (!requireNamespace("RoCE", quietly = TRUE)) {
  stop("the tested RoCE package is unavailable for the smoke audit.",
       call. = FALSE)
}
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
expected_nlambda <- suppressWarnings(as.integer(
  Sys.getenv("ROCE_EXPECT_NLAMBDA", "100")
))
if (length(expected_nlambda) != 1L || is.na(expected_nlambda) ||
    expected_nlambda < 2L) {
  stop("ROCE_EXPECT_NLAMBDA must be one integer of at least 2.",
       call. = FALSE)
}
result_path <- file.path(output_root, "task_000001.csv")
if (!file.exists(result_path)) {
  stop("smoke result not found: ", result_path, call. = FALSE)
}

result <- read.csv(result_path, stringsAsFactors = FALSE)
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

expected_methods <- c(
  "one_round_crossfit", "target_only", "sample_size", "inverse_variance",
  "federated_dr", "pooled_dr", "one_round_crossfit_ate_armwise",
  "one_round_crossfit_ate", "target_only_ate", "sample_size_ate",
  "inverse_variance_ate", "federated_dr_ate", "pooled_dr_ate"
)
required_columns <- c(
  "sim_id", "method", "estimate", "se", "bias", "coverage", "ci_width",
  "estimand_scope", "truth", "ci_lower", "ci_upper", "config", "p", "K",
  "outcome_family",
  "rho", "cutoff", "n_site", "n_folds", "aggregation_lambda",
  "aggregation_cutoff", "primary_cutoff", "nlambda_init",
  "nuisance_lambda_rule",
  "n_bootstrap", "M_tau", "M_tau_inference", "target_anchor_nlambda",
  "min_site_arm_outcome_cell_n", "min_target_arm_outcome_cell_n",
  "n_site_arm_outcome_cells_below_8",
  "dr_weight_n", "dr_weight_n_clipped",
  "dr_weight_fraction_clipped", "dr_weight_max_site_fraction_clipped",
  "dr_weight_min_before_clipping", "dr_weight_max_before_clipping",
  roce_resource_metadata_columns(),
  roce_scheduler_provenance_columns(),
  "task_id", "experiment", "package_library", "package_fingerprint",
  "workflow_fingerprint", "manifest_fingerprint",
  "task_elapsed_seconds"
)
missing_columns <- setdiff(required_columns, names(result))
add_check(
  "required_schema",
  length(missing_columns) == 0L,
  if (length(missing_columns) == 0L) "complete" else {
    paste(missing_columns, collapse = ",")
  }
)
add_check(
  "exact_method_set",
  nrow(result) == length(expected_methods) &&
    identical(sort(result$method), sort(expected_methods)) &&
    !anyDuplicated(result$method),
  paste(result$method, collapse = ",")
)

if (length(missing_columns) == 0L) {
  expected_cv_threads <- roce_read_positive_integer_env(
    "ROCE_EXPECT_NUISANCE_CV_THREADS", 5L, 5L
  )
  resource_metadata_ok <- isTRUE(tryCatch(
    roce_validate_result_resource_metadata(
      result, expected_threads = expected_cv_threads
    ),
    error = function(error) FALSE
  ))
  add_check(
    "slurm_resource_metadata",
    resource_metadata_ok,
    sprintf(
      "allocated=%s;threads=%s;workers=%s;parallel=%s;full=%s",
      paste(unique(result$allocated_cores), collapse = ","),
      paste(unique(result$nuisance_cv_threads), collapse = ","),
      paste(unique(result$source_workers_per_arm), collapse = ","),
      paste(unique(result$parallel_treatment_arms), collapse = ","),
      paste(unique(result$fully_parallel_cores), collapse = ",")
    )
  )
  scheduler_provenance_ok <- isTRUE(tryCatch(
    roce_validate_scheduler_provenance(result, require_slurm = TRUE),
    error = function(error) FALSE
  ))
  add_check(
    "scheduler_provenance",
    scheduler_provenance_ok,
    paste(unique(result$slurm_job_id), collapse = ",")
  )
  provenance <- roce_result_provenance(result, project_library)
  add_check(
    "tested_package_provenance",
    provenance$passed,
    provenance$detail
  )
  expected_workflow <- tolower(trimws(Sys.getenv(
    "ROCE_WORKFLOW_FINGERPRINT", ""
  )))
  observed_workflow <- unique(tolower(trimws(
    as.character(result$workflow_fingerprint)
  )))
  observed_manifest <- unique(tolower(trimws(
    as.character(result$manifest_fingerprint)
  )))
  add_check(
    "simulation_workflow_provenance",
    length(observed_workflow) == 1L &&
      grepl("^[0-9a-f]{64}$", observed_workflow) &&
      identical(observed_workflow, expected_workflow),
    paste(observed_workflow, collapse = ",")
  )
  add_check(
    "manifest_provenance",
    length(observed_manifest) == 1L &&
      grepl("^[0-9a-f]{64}$", observed_manifest),
    paste(observed_manifest, collapse = ",")
  )
  numeric_columns <- c(
    "estimate", "se", "bias", "ci_width", "truth", "ci_lower", "ci_upper",
    "task_elapsed_seconds"
  )
  add_check(
    "finite_core_values",
    all(is.finite(unlist(result[numeric_columns]))),
    paste("rows=", nrow(result), sep = "")
  )
  add_check(
    "positive_standard_errors_and_widths",
    all(result$se > 0) && all(result$ci_width > 0),
    paste(signif(range(result$se), 6), collapse = ",")
  )
  bias_error <- max(abs(result$bias - (result$estimate - result$truth)))
  ci_error <- max(abs(c(
    result$ci_lower - (result$estimate - stats::qnorm(0.975) * result$se),
    result$ci_upper - (result$estimate + stats::qnorm(0.975) * result$se),
    result$ci_width - (result$ci_upper - result$ci_lower)
  )))
  coverage_text <- toupper(trimws(as.character(result$coverage)))
  coverage_valid <- !is.na(result$coverage) &
    coverage_text %in% c("TRUE", "T", "1", "FALSE", "F", "0")
  coverage <- if (is.logical(result$coverage)) {
    as.logical(result$coverage)
  } else {
    coverage_text %in% c("TRUE", "T", "1")
  }
  expected_coverage <- result$truth >= result$ci_lower &
    result$truth <= result$ci_upper
  add_check("bias_identity", bias_error <= tolerance,
            sprintf("max_abs_error=%.3g", bias_error))
  add_check("ci_arithmetic", ci_error <= tolerance,
            sprintf("max_abs_error=%.3g", ci_error))
  add_check(
    "coverage_ci_identity",
    all(coverage_valid) && identical(coverage, expected_coverage),
    paste(as.integer(coverage), collapse = ",")
  )

  expected_scope <- ifelse(
    RoCE:::.is_tate_method(result$method), "tate", "treated_mean"
  )
  add_check(
    "estimand_scope_matches_method",
    identical(as.character(result$estimand_scope), expected_scope),
    paste(result$estimand_scope, collapse = ",")
  )
  truth_counts <- vapply(
    split(result$truth, result$estimand_scope),
    function(values) length(unique(values)),
    integer(1L)
  )
  add_check(
    "one_truth_per_estimand",
    length(truth_counts) == 2L && all(truth_counts == 1L),
    paste(names(truth_counts), truth_counts, sep = "=", collapse = ",")
  )

  metadata_ok <- all(result$sim_id == 1L) && all(result$task_id == 1L) &&
    all(result$experiment == "single_task_smoke") &&
    all(result$config == "C3") && all(result$p == 100L) &&
    all(result$outcome_family == "binomial") &&
    all(result$K == 4L) && all(result$rho == 0) &&
    all(abs(result$cutoff - expected_primary_cutoff) <= tolerance) &&
    all(abs(result$primary_cutoff - expected_primary_cutoff) <= tolerance) &&
    all(result$n_site == 1000L) && all(result$n_folds == 5L) &&
    all(abs(
      result$aggregation_lambda - 1 / expected_primary_cutoff
    ) <= tolerance) &&
    all(abs(
      result$aggregation_cutoff - expected_primary_cutoff
    ) <= tolerance) &&
    all(result$nlambda_init == expected_nlambda) &&
    all(result$nuisance_lambda_rule == "min") &&
    all(result$target_anchor_nlambda == 100L) &&
    all(result$n_bootstrap == 5000L)
  add_check(
    "p100_c3_smoke_metadata",
    metadata_ok,
    sprintf(
      "sim=%s;task=%s;C=%s;p=%s;K=%s;rho=%s;nlambda=%s",
      result$sim_id[[1L]], result$task_id[[1L]], result$config[[1L]],
      result$p[[1L]], result$K[[1L]], result$rho[[1L]],
      result$nlambda_init[[1L]]
    )
  )

  comparison_rows <- grepl(
    "^(sample_size|inverse_variance|federated_dr|pooled_dr)(_ate)?$",
    result$method
  )
  bootstrap_columns <- c(
    "comparison_variance_method", "comparison_se_analytic",
    "comparison_se_bootstrap", "comparison_n_bootstrap"
  )
  bootstrap_schema_ok <- all(bootstrap_columns %in% names(result))
  bootstrap_ok <- bootstrap_schema_ok && any(comparison_rows) &&
    all(result$comparison_variance_method[comparison_rows] == "bootstrap") &&
    all(is.finite(result$comparison_se_analytic[comparison_rows])) &&
    all(is.finite(result$comparison_se_bootstrap[comparison_rows])) &&
    all(result$comparison_se_analytic[comparison_rows] > 0) &&
    all(result$comparison_se_bootstrap[comparison_rows] > 0) &&
    all(result$comparison_n_bootstrap[comparison_rows] == 5000L) &&
    max(abs(
      result$se[comparison_rows] -
        result$comparison_se_bootstrap[comparison_rows]
    )) <= tolerance
  add_check(
    "comparison_bootstrap_metadata",
    bootstrap_ok,
    if (bootstrap_schema_ok) {
      sprintf(
        "rows=%d;n_bootstrap=%s;max_se_error=%.3g",
        sum(comparison_rows),
        paste(unique(result$comparison_n_bootstrap[comparison_rows]), collapse = ","),
        max(abs(
          result$se[comparison_rows] -
            result$comparison_se_bootstrap[comparison_rows]
        ))
      )
    } else {
      paste(setdiff(bootstrap_columns, names(result)), collapse = ",")
    }
  )

  truncation_ok <- all(result$M_tau == 5) &&
    all(result$M_tau_inference == 5)
  add_check(
    "truncation_metadata",
    truncation_ok,
    sprintf(
      "M_tau=%s;M_tau_inference=%s",
      result$M_tau[[1L]], result$M_tau_inference[[1L]]
    )
  )
}

direct <- result[result$method == "one_round_crossfit_ate", , drop = FALSE]
required_direct_diagnostics <- unique(c(
  RoCE:::.direct_tate_aggregation_diagnostic_columns(),
  "face_initial_dr_nonconverged",
  "face_calibrated_dr_nonconverged", "face_calibrated_outcome_nonconverged",
  "face_initial_dr_cv_invalid_fold_fits",
  "face_initial_dr_cv_invalid_lambdas",
  "face_calibrated_dr_cv_invalid_fold_fits",
  "face_calibrated_dr_cv_invalid_lambdas",
  "face_calibrated_outcome_cv_invalid_fold_fits",
  "face_calibrated_outcome_cv_invalid_lambdas"
))
missing_direct <- setdiff(required_direct_diagnostics, names(result))
direct_ok <- nrow(direct) == 1L && length(missing_direct) == 0L
direct_values <- if (direct_ok) {
  unlist(direct[required_direct_diagnostics])
} else {
  numeric(0)
}
add_check(
  "direct_tate_diagnostic_schema_and_finiteness",
  direct_ok && all(is.finite(direct_values)),
  if (length(missing_direct) == 0L) "complete" else {
    paste(missing_direct, collapse = ",")
  }
)
if (direct_ok && all(is.finite(direct_values))) {
  add_check(
    "direct_tate_weight_and_wald_ranges",
    abs(direct$target_anchor_weight) <= 10 &&
      direct$mean_abs_source_weight <= direct$max_abs_source_weight &&
      direct$max_abs_source_weight <= 10 &&
      direct$mean_wald_statistic >= 0 &&
      direct$max_wald_statistic >= direct$mean_wald_statistic &&
      direct$penalized_source_fold_fraction >= 0 &&
      direct$penalized_source_fold_fraction <= 1,
    sprintf(
      "anchor=%.4g;max_weight=%.4g;max_wald=%.4g;penalty_fraction=%.4g",
      direct$target_anchor_weight, direct$max_abs_source_weight,
      direct$max_wald_statistic, direct$penalized_source_fold_fraction
    )
  )
  nonconverged <- sum(unlist(direct[c(
    "face_initial_dr_nonconverged", "face_calibrated_dr_nonconverged",
    "face_calibrated_outcome_nonconverged"
  )]))
  add_check(
    "direct_tate_nuisance_convergence",
    nonconverged == 0,
    sprintf("nonconverged_fits=%g", nonconverged)
  )
  cv_exclusion_columns <- c(
    "face_initial_dr_cv_invalid_fold_fits",
    "face_initial_dr_cv_invalid_lambdas",
    "face_calibrated_dr_cv_invalid_fold_fits",
    "face_calibrated_dr_cv_invalid_lambdas",
    "face_calibrated_outcome_cv_invalid_fold_fits",
    "face_calibrated_outcome_cv_invalid_lambdas"
  )
  cv_exclusions <- unlist(direct[cv_exclusion_columns])
  add_check(
    "direct_tate_cv_candidate_exclusion_diagnostics",
    all(is.finite(cv_exclusions)) && all(cv_exclusions >= 0),
    sprintf(
      "invalid_fold_fits=%g;invalid_lambdas=%g",
      sum(cv_exclusions[c(1, 3, 5)]),
      sum(cv_exclusions[c(2, 4, 6)])
    )
  )
  stability_columns <- c(
    "face_initial_dr_line_search_failures",
    "face_calibrated_dr_line_search_failures",
    "face_initial_dr_support_floor_applied",
    "face_calibrated_dr_support_floor_applied",
    "face_max_initial_dr_abs_coefficient",
    "face_max_calibrated_dr_abs_coefficient",
    "face_max_initial_dr_support_floor",
    "face_max_calibrated_dr_support_floor"
  )
  stability_schema_ok <- all(stability_columns %in% names(direct))
  add_check(
    "direct_tate_nuisance_stability_schema",
    stability_schema_ok,
    paste(setdiff(stability_columns, names(direct)), collapse = ",")
  )
  if (stability_schema_ok) {
    line_search_failures <- sum(unlist(direct[c(
      "face_initial_dr_line_search_failures",
      "face_calibrated_dr_line_search_failures"
    )]))
    max_density_ratio_coefficient <- max(unlist(direct[c(
      "face_max_initial_dr_abs_coefficient",
      "face_max_calibrated_dr_abs_coefficient"
    )]))
    support_values <- unlist(direct[c(
      "face_initial_dr_support_floor_applied",
      "face_calibrated_dr_support_floor_applied",
      "face_max_initial_dr_support_floor",
      "face_max_calibrated_dr_support_floor"
    )])
    add_check(
      "direct_tate_density_ratio_stability",
      all(is.finite(c(
        line_search_failures, max_density_ratio_coefficient, support_values
      ))) && line_search_failures == 0 &&
        max_density_ratio_coefficient < 0.9 * RoCE:::PARAM_MAX &&
        all(support_values >= 0),
      sprintf(
        "line_search_failures=%g;max_abs_coefficient=%.6g;support_activations=%g",
        line_search_failures, max_density_ratio_coefficient,
        direct$face_initial_dr_support_floor_applied +
          direct$face_calibrated_dr_support_floor_applied
      )
    )
  }
}

optional_diagnostics <- intersect(c(
  "max_weight_optimizer_iterations", "max_weight_psd_ridge",
  "weight_psd_ridge_fold_fraction", "inference_logit_truncated",
  "inference_logit_truncation_fraction", "inference_max_abs_logit",
  "inference_safety_clip_count"
), names(result))
if (length(optional_diagnostics) > 0L && nrow(direct) == 1L) {
  optional_values <- unlist(direct[optional_diagnostics])
  optional_ok <- all(is.finite(optional_values))
  if ("max_weight_optimizer_iterations" %in% optional_diagnostics) {
    optional_ok <- optional_ok && direct$max_weight_optimizer_iterations >= 1
  }
  if ("max_weight_psd_ridge" %in% optional_diagnostics) {
    optional_ok <- optional_ok && direct$max_weight_psd_ridge >= 0
  }
  if ("weight_psd_ridge_fold_fraction" %in% optional_diagnostics) {
    optional_ok <- optional_ok && direct$weight_psd_ridge_fold_fraction >= 0 &&
      direct$weight_psd_ridge_fold_fraction <= 1
  }
  if ("inference_logit_truncation_fraction" %in% optional_diagnostics) {
    optional_ok <- optional_ok &&
      direct$inference_logit_truncation_fraction >= 0 &&
      direct$inference_logit_truncation_fraction <= 1
  }
  if ("inference_safety_clip_count" %in% optional_diagnostics) {
    optional_ok <- optional_ok && direct$inference_safety_clip_count == 0
  }
  add_check(
    "optional_optimizer_and_clipping_diagnostics",
    optional_ok,
    paste(optional_diagnostics, signif(optional_values, 5),
          sep = "=", collapse = ",")
  )
}

sensitivity_root <- Sys.getenv(
  "ROCE_SENSITIVITY_OUTPUT_ROOT",
  paste0(sub("[/]+$", "", output_root), "_reused_sensitivity")
)
sensitivity_path <- file.path(sensitivity_root, "task_000001.csv")
add_check(
  "reused_sensitivity_sidecar_exists",
  file.exists(sensitivity_path),
  sensitivity_path
)
if (file.exists(sensitivity_path)) {
  sensitivity <- read.csv(sensitivity_path, stringsAsFactors = FALSE)
  sensitivity_required <- c(
    "task_id", "sim_id", "method", "estimate", "se", "truth", "bias",
    "coverage", "ci_lower", "ci_upper", "ci_width", "config", "p", "K",
    "rho", "cutoff", "aggregation_lambda", "aggregation_cutoff",
    "primary_cutoff",
    "M_tau", "M_tau_inference", "sensitivity_kind", "experiment",
    "primary_experiment", "dgp_type", "outcome_family",
    "heterogeneity_type", "estimand_type",
    "nuisance_lambda_rule", "min_site_arm_outcome_cell_n",
    "min_target_arm_outcome_cell_n", "n_site_arm_outcome_cells_below_8",
    "reused_primary_task", "nuisance_refit_count", "package_library",
    "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
    roce_scheduler_provenance_columns()
  )
  missing_sensitivity <- setdiff(sensitivity_required, names(sensitivity))
  add_check(
    "reused_sensitivity_schema",
    length(missing_sensitivity) == 0L,
    if (length(missing_sensitivity) == 0L) "complete" else {
      paste(missing_sensitivity, collapse = ",")
    }
  )
  if (length(missing_sensitivity) == 0L) {
    sensitivity_grid <- unique(sensitivity[c("cutoff", "M_tau_inference")])
    observed_grid <- sensitivity_grid[order(
      sensitivity_grid$M_tau_inference, sensitivity_grid$cutoff
    ), , drop = FALSE]
    expected_grid <- rbind(
      data.frame(cutoff = c(1, 1.5, 2, 2.5, 3), M_tau_inference = 5),
      data.frame(
        cutoff = expected_primary_cutoff,
        M_tau_inference = c(4, 6, Inf)
      )
    )
    expected_grid <- expected_grid[order(
      expected_grid$M_tau_inference, expected_grid$cutoff
    ), , drop = FALSE]
    rownames(observed_grid) <- NULL
    rownames(expected_grid) <- NULL
    method_counts <- table(sensitivity$method)
    add_check(
      "reused_sensitivity_exact_grid",
      nrow(sensitivity) == 16L &&
        identical(observed_grid, expected_grid) &&
        identical(
          as.integer(method_counts[c(
            "one_round_crossfit_ate", "target_only_ate"
          )]),
          c(8L, 8L)
        ) &&
        !anyDuplicated(sensitivity[c(
          "method", "cutoff", "M_tau_inference"
        )]),
      sprintf("rows=%d;settings=%d", nrow(sensitivity), nrow(observed_grid))
    )
    sensitivity_provenance <- roce_result_provenance(
      sensitivity, project_library
    )
    add_check(
      "reused_sensitivity_provenance",
      sensitivity_provenance$passed &&
        isTRUE(tryCatch(
          roce_validate_scheduler_provenance(
            sensitivity, require_slurm = TRUE
          ),
          error = function(error) FALSE
        )) &&
        identical(
          unique(tolower(sensitivity$workflow_fingerprint)),
          unique(tolower(result$workflow_fingerprint))
        ) &&
        identical(
          unique(tolower(sensitivity$manifest_fingerprint)),
          unique(tolower(result$manifest_fingerprint))
        ),
      sensitivity_provenance$detail
    )
    arithmetic_error <- max(abs(c(
      sensitivity$bias - (sensitivity$estimate - sensitivity$truth),
      sensitivity$ci_lower -
        (sensitivity$estimate - stats::qnorm(0.975) * sensitivity$se),
      sensitivity$ci_upper -
        (sensitivity$estimate + stats::qnorm(0.975) * sensitivity$se),
      sensitivity$ci_width -
        (sensitivity$ci_upper - sensitivity$ci_lower)
    )))
    add_check(
      "reused_sensitivity_arithmetic_and_reuse",
      all(is.finite(unlist(sensitivity[c(
        "estimate", "se", "truth", "bias", "ci_lower", "ci_upper",
        "ci_width"
      )]))) &&
      all(sensitivity$se > 0) && all(sensitivity$M_tau == 5) &&
        all(sensitivity$experiment == "c3_reused_sensitivity") &&
        all(sensitivity$primary_experiment == "single_task_smoke") &&
        all(sensitivity$dgp_type == "face") &&
        all(sensitivity$outcome_family == "binomial") &&
        all(sensitivity$heterogeneity_type == "none") &&
        all(sensitivity$estimand_type == "superpopulation") &&
        all(sensitivity$nuisance_lambda_rule == "min") &&
        all(abs(
          sensitivity$primary_cutoff - expected_primary_cutoff
        ) <= tolerance) &&
        all(sensitivity$aggregation_cutoff == sensitivity$cutoff) &&
        all(abs(
          sensitivity$aggregation_lambda - 1 / sensitivity$cutoff
        ) <= tolerance) &&
        all(is.finite(sensitivity$min_site_arm_outcome_cell_n)) &&
        all(is.finite(sensitivity$min_target_arm_outcome_cell_n)) &&
        all(is.finite(sensitivity$n_site_arm_outcome_cells_below_8)) &&
        all(sensitivity$reused_primary_task) &&
        all(sensitivity$nuisance_refit_count == 0L) &&
        arithmetic_error <= tolerance,
      sprintf("max_arithmetic_error=%.3g", arithmetic_error)
    )
    reference <- sensitivity[
      sensitivity$sensitivity_kind == "reference_identity",
      , drop = FALSE
    ]
    primary_reference <- result[
      result$method %in% c("one_round_crossfit_ate", "target_only_ate"),
      , drop = FALSE
    ]
    reference <- reference[order(reference$method), , drop = FALSE]
    primary_reference <- primary_reference[
      order(primary_reference$method), , drop = FALSE
    ]
    identity_columns <- c(
      "method", "estimate", "se", "truth", "bias", "coverage",
      "ci_lower", "ci_upper", "ci_width"
    )
    add_check(
      "reused_sensitivity_reference_identity",
      nrow(reference) == 2L && nrow(primary_reference) == 2L &&
        isTRUE(all.equal(
          reference[identity_columns], primary_reference[identity_columns],
          tolerance = tolerance, check.attributes = FALSE
        )),
      sprintf("reference_rows=%d", nrow(reference))
    )
    target_sensitivity <- sensitivity[
      sensitivity$method == "target_only_ate", , drop = FALSE
    ]
    add_check(
      "reused_sensitivity_target_invariance",
      length(unique(signif(target_sensitivity$estimate, 14))) == 1L &&
        length(unique(signif(target_sensitivity$se, 14))) == 1L,
      sprintf(
        "estimate_range=%.3g;se_range=%.3g",
        diff(range(target_sensitivity$estimate)),
        diff(range(target_sensitivity$se))
      )
    )
  }
}

audit <- do.call(rbind, checks)
audit_path <- file.path(output_root, "direct_tate_smoke_audit.csv")
write.csv(audit, audit_path, row.names = FALSE)
print(audit, row.names = FALSE)

failed <- audit[!audit$passed, , drop = FALSE]
if (nrow(failed) > 0L) {
  stop(
    "TATE smoke audit failed: ",
    paste(failed$check, collapse = ", "),
    call. = FALSE
  )
}

# Commit a compact production gate only after every smoke check passes.  The
# production submitter validates these values against the exact immutable
# package, task workflow, and nuisance-grid size it is about to use.
single_gate_value <- function(values, name) {
  values <- unique(values)
  if (length(values) != 1L || is.na(values) || !nzchar(as.character(values))) {
    stop("smoke gate requires one nonmissing ", name, ".", call. = FALSE)
  }
  as.character(values)
}
gate_values <- c(
  "implementation_audit=passed",
  paste0(
    "package_fingerprint=",
    single_gate_value(result$package_fingerprint, "package fingerprint")
  ),
  paste0(
    "workflow_fingerprint=",
    single_gate_value(result$workflow_fingerprint, "workflow fingerprint")
  ),
  paste0(
    "audit_driver_md5=",
    unname(tools::md5sum(file.path(
      "scripts", "slurm", "audit_direct_tate_smoke.R"
    )))
  ),
  paste0(
    "submit_driver_md5=",
    unname(tools::md5sum(file.path(
      "scripts", "slurm", "submit_main_direct_tate.sh"
    )))
  ),
  paste0("smoke_result_md5=", unname(tools::md5sum(result_path))),
  paste0(
    "smoke_sensitivity_result_md5=",
    unname(tools::md5sum(sensitivity_path))
  ),
  paste0(
    "nlambda_init=",
    single_gate_value(as.integer(result$nlambda_init), "nuisance-grid size")
  ),
  paste0("primary_cutoff=", expected_primary_cutoff),
  "p=100",
  "config=C3",
  "K=4",
  "rho=0"
)
if (any(lengths(lapply(strsplit(gate_values, "=", fixed = TRUE), identity)) !=
        2L)) {
  stop("internal error while constructing the smoke pass gate.", call. = FALSE)
}
gate_path <- file.path(output_root, "direct_tate_smoke_audit_passed.txt")
gate_temporary <- tempfile(
  pattern = ".direct_tate_smoke_audit_passed_",
  tmpdir = output_root, fileext = ".tmp"
)
writeLines(gate_values, gate_temporary)
if (!file.rename(gate_temporary, gate_path)) {
  unlink(gate_temporary)
  stop("failed to atomically write the smoke pass gate.", call. = FALSE)
}
message(
  "[done] TATE smoke implementation audit passed. ",
  "One replicate is not used to judge Monte Carlo coverage or RMSE."
)
