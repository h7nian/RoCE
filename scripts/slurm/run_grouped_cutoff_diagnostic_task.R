#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
source(file.path("scripts", "slurm", "grouped_cutoff_diagnostic_helpers.R"))
if (length(args) != 2L) {
  stop("usage: run_grouped_cutoff_diagnostic_task.R MANIFEST.csv ARRAY_ROW",
       call. = FALSE)
}

manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
row_id <- suppressWarnings(as.integer(args[[2L]]))
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (is.na(row_id) || row_id < 1L || row_id > nrow(manifest)) {
  stop("array row is outside the grouped-cutoff manifest.", call. = FALSE)
}
task <- manifest[row_id, , drop = FALSE]
required_columns <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho_values",
  "cutoffs", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
  "M_tau", "M_tau_inference"
)
missing_columns <- setdiff(required_columns, names(task))
if (length(missing_columns) > 0L) {
  stop("manifest is missing: ", paste(missing_columns, collapse = ", "),
       call. = FALSE)
}
if (as.integer(task$p) != 100L ||
    !as.integer(task$K) %in% c(2L, 4L, 8L) ||
    !as.character(task$config) %in% c("C1", "C2", "C3")) {
  stop("grouped cutoff diagnostics require p=100 and an approved C/K setting.",
       call. = FALSE)
}
if (as.integer(task$n_site) != 1000L ||
    as.integer(task$n_folds) != 5L ||
    as.integer(task$nlambda_init) != 100L ||
    as.integer(task$n_bootstrap) != 5000L ||
    as.numeric(task$M_tau) != 5 ||
    as.numeric(task$M_tau_inference) != 5) {
  stop(
    paste0(
      "grouped cutoff diagnostics must use the locked production settings: ",
      "n_site=1000, n_folds=5, nlambda_init=100, n_bootstrap=5000, ",
      "M_tau=M_tau_inference=5."
    ),
    call. = FALSE
  )
}

rho_values <- roce_parse_semicolon_numeric_grid(
  task$rho_values, "rho_values", require_zero_first = TRUE
)
cutoffs <- roce_parse_semicolon_numeric_grid(
  task$cutoffs, "cutoffs", require_positive = TRUE
)
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))

primary_cutoff <- 1 / RoCE:::AGG_WALD_LAMBDA
if (!any(abs(cutoffs - primary_cutoff) <= 1e-12)) {
  stop("cutoffs must contain the tested package primary cutoff.",
       call. = FALSE)
}
if ("primary_cutoff" %in% names(task) &&
    abs(as.numeric(task$primary_cutoff) - primary_cutoff) > 1e-12) {
  stop("manifest primary cutoff differs from the tested package.",
       call. = FALSE)
}

package_provenance <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment("ROCE_WORKFLOW_FINGERPRINT")
manifest_fingerprint <- roce_sha256_environment("ROCE_MANIFEST_FINGERPRINT")
output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT",
  file.path(dirname(manifest_path), "raw")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_root, sprintf("task_%06d.csv", as.integer(task$task_id))
)
if (file.exists(output_path)) {
  message("[skip] output already exists: ", output_path)
  quit(save = "no", status = 0L)
}
lock_path <- file.path(
  output_root, "locks", sprintf("task_%06d.lock", as.integer(task$task_id))
)
release_lock <- roce_claim_task_lock(lock_path, c(
  paste0("task_id=", task$task_id),
  paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
  paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
))
on.exit(release_lock(), add = TRUE)

allocated_cores <- roce_read_positive_integer_env("SLURM_CPUS_PER_TASK", 1L)
nuisance_cv_threads <- roce_read_positive_integer_env(
  "ROCE_NUISANCE_CV_THREADS", 1L, as.integer(task$n_folds)
)
resource_plan <- roce_slurm_resource_plan(
  source_count = as.integer(task$K),
  n_folds = as.integer(task$n_folds),
  allocated_cores = allocated_cores,
  nuisance_cv_threads = nuisance_cv_threads
)
positive_rho_workers <- roce_read_positive_integer_env(
  "ROCE_POSITIVE_RHO_WORKERS",
  min(length(rho_values) - 1L,
      max(1L, allocated_cores %/% nuisance_cv_threads)),
  length(rho_values) - 1L
)

message(sprintf(
  paste0(
    "[grouped-cutoff %d] sim=%d C=%s p=100 K=%d rhos=%s cutoffs=%s ",
    "allocated=%d cv_threads=%d source_workers=%d positive_rho_workers=%d"
  ),
  task$task_id, task$sim_id, task$config, task$K,
  paste(rho_values, collapse = ","), paste(cutoffs, collapse = ","),
  allocated_cores, nuisance_cv_threads, resource_plan$source_workers,
  positive_rho_workers
))

started_at <- proc.time()[["elapsed"]]
simulation_args <- list(
  sim_id = as.integer(task$sim_id),
  n_total = as.integer(task$n_site) * (as.integer(task$K) + 1L),
  K = as.integer(task$K),
  p = 100L,
  config = as.character(task$config),
  methods = "one_round_crossfit",
  verbose = FALSE,
  n_cores_internal = resource_plan$source_workers,
  nlambda_init = as.integer(task$nlambda_init),
  estimand_type = "superpopulation",
  outcome_type = "binary",
  n_folds = as.integer(task$n_folds),
  aggregation_lambda = 1 / primary_cutoff,
  n_bootstrap = as.integer(task$n_bootstrap),
  M_tau = as.numeric(task$M_tau),
  M_tau_inference = as.numeric(task$M_tau_inference),
  estimate_ate = TRUE,
  parallel_treatment_arms = resource_plan$parallel_treatment_arms,
  dgp_type = "face"
)
grouped <- RoCE:::.run_face_rho_group(
  simulation_args = simulation_args,
  rho_values = rho_values,
  changed_sources = "s1",
  artifact_rhos = rho_values,
  positive_rho_workers = positive_rho_workers
)

make_result_row <- function(fit, method, truth, rho, cutoff) {
  estimate <- as.numeric(fit$estimate)
  se <- as.numeric(fit$se)
  z <- stats::qnorm(0.975)
  n_total <- as.integer(task$n_site) * (as.integer(task$K) + 1L)
  result <- RoCE:::.make_simulation_result_row(
    sim_id = as.integer(task$sim_id),
    method = method,
    estimate = estimate,
    se = se,
    truth = truth,
    n_total = n_total,
    K = as.integer(task$K),
    p = 100L,
    config = as.character(task$config),
    heterogeneity_type = if (rho > 0) "one_deviated_source" else "none",
    estimand_type = "superpopulation"
  )
  result$task_id <- as.integer(task$task_id)
  result$experiment <- as.character(task$experiment)
  result$rho <- rho
  result$cutoff <- cutoff
  result$aggregation_lambda <- 1 / cutoff
  result$aggregation_cutoff <- cutoff
  result$estimand_scope <- "tate"
  result$truth <- truth
  result$ci_lower <- estimate - z * se
  result$ci_upper <- estimate + z * se
  result$n_site <- as.integer(task$n_site)
  result$n_folds <- as.integer(task$n_folds)
  result$nlambda_init <- as.integer(task$nlambda_init)
  result$nuisance_lambda_rule <- "min"
  result$n_bootstrap <- as.integer(task$n_bootstrap)
  result$M_tau <- as.numeric(task$M_tau)
  result$M_tau_inference <- as.numeric(task$M_tau_inference)
  result$dgp_type <- "face"
  result$outcome_family <- "binomial"
  result$primary_cutoff <- primary_cutoff
  result
}

rows <- list()
row_index <- 0L
for (rho_index in seq_along(rho_values)) {
  rho <- rho_values[[rho_index]]
  rho_key <- format(rho, scientific = FALSE, trim = TRUE)
  artifacts <- grouped$artifacts[[rho_key]]
  fitted_tate <- artifacts$direct_tate_results$one_round_crossfit
  tate_truth <- as.numeric(artifacts$tate_truth)
  if (is.null(fitted_tate) || is.null(artifacts$data_split) ||
      length(tate_truth) != 1L || !is.finite(tate_truth)) {
    stop("grouped rho run did not retain a complete artifact for rho=", rho,
         call. = FALSE)
  }
  grid <- data.frame(
    cutoff = cutoffs,
    M_tau_inference = rep(as.numeric(task$M_tau_inference), length(cutoffs))
  )
  evaluated <- reaggregate_tate_sensitivity_grid(
    data_split = artifacts$data_split,
    fitted_tate = fitted_tate,
    sensitivity_grid = grid,
    verbose = FALSE
  )
  if (evaluated$n_nuisance_refits != 0L ||
      evaluated$n_inference_refreshes != 0L) {
    stop("cutoff-only reaggregation unexpectedly refreshed a nuisance or IF fit.",
         call. = FALSE)
  }

  primary <- grouped$results[[rho_index]]
  primary_direct <- primary[primary$method == "one_round_crossfit_ate", , drop = FALSE]
  if (nrow(primary_direct) != 1L) {
    stop("rho-group result is missing its unique TATE row.",
         call. = FALSE)
  }
  nuisance_diagnostics <- c(
    RoCE:::.summarize_tate_crossfit_timing(
      fitted_tate$arm_results$mu1, fitted_tate$arm_results$mu0
    ),
    RoCE:::.summarize_tate_nuisance_diagnostics(
      fitted_tate$arm_results$mu1, fitted_tate$arm_results$mu0
    )
  )
  cell_diagnostics <- RoCE:::.binary_cell_diagnostics(
    artifacts$data_split, "binomial"
  )

  for (cutoff_index in seq_along(cutoffs)) {
    cutoff <- cutoffs[[cutoff_index]]
    fit <- evaluated$results[[cutoff_index]]
    face <- make_result_row(
      fit, "one_round_crossfit_ate", tate_truth, rho, cutoff
    )
    aggregation_diagnostics <-
      RoCE:::.summarize_tate_aggregation_diagnostics(fit)
    all_diagnostics <- c(aggregation_diagnostics, nuisance_diagnostics)
    if (anyDuplicated(names(all_diagnostics))) {
      stop("aggregation and nuisance diagnostic names must be unique.",
           call. = FALSE)
    }
    for (name in names(all_diagnostics)) {
      face[[name]] <- all_diagnostics[[name]]
    }
    for (name in names(cell_diagnostics)) {
      face[[name]] <- cell_diagnostics[[name]]
    }

    fold_weights <- as.matrix(fit$fold_weights)
    fold_wald <- as.matrix(fit$fold_wald_statistics)
    fold_penalty <- as.matrix(fit$fold_penalty_coefficients)
    if (!identical(dim(fold_weights), dim(fold_wald)) ||
        !identical(dim(fold_weights), dim(fold_penalty)) ||
        ncol(fold_weights) != as.integer(task$K) ||
        any(!is.finite(fold_weights)) || any(!is.finite(fold_wald)) ||
        any(!is.finite(fold_penalty)) || any(fold_penalty < 0)) {
      stop("source-level fold diagnostics have invalid dimensions or values.",
           call. = FALSE)
    }
    face$source1_mean_weight <- mean(fold_weights[, 1L])
    face$source1_mean_abs_weight <- mean(abs(fold_weights[, 1L]))
    face$source1_zero_weight_fold_fraction <-
      mean(abs(fold_weights[, 1L]) <= sqrt(.Machine$double.eps))
    face$source1_mean_wald_statistic <- mean(fold_wald[, 1L])
    face$source1_max_wald_statistic <- max(fold_wald[, 1L])
    face$source1_penalty_activation_fraction <- mean(fold_penalty[, 1L] > 0)
    face$informative_sources_mean_abs_weight <- if (ncol(fold_weights) > 1L) {
      mean(abs(fold_weights[, -1L, drop = FALSE]))
    } else {
      NA_real_
    }
    if (is.null(fit$source_estimates) ||
        !"s1" %in% names(fit$source_estimates)) {
      stop("TATE fit is missing the named s1 source estimate.",
           call. = FALSE)
    }
    source1_estimate <- as.numeric(fit$source_estimates[["s1"]])
    if (length(source1_estimate) != 1L || !is.finite(source1_estimate)) {
      stop("the s1 source estimate must be one finite value.", call. = FALSE)
    }
    face$source1_estimate <- source1_estimate
    face$source1_minus_target_estimate <-
      source1_estimate - as.numeric(fit$target_only$estimate)
    phase1b <- fit$intermediates$phase1b
    if (is.null(phase1b$avg_source_est) ||
        is.null(phase1b$avg_target_est) ||
        ncol(as.matrix(phase1b$avg_source_est)) != as.integer(task$K) ||
        nrow(as.matrix(phase1b$avg_source_est)) !=
          length(phase1b$avg_target_est) ||
        any(!is.finite(phase1b$avg_source_est)) ||
        any(!is.finite(phase1b$avg_target_est))) {
      stop("training-fold target/source discrepancies are incomplete.",
           call. = FALSE)
    }
    face$source1_mean_training_discrepancy <- mean(
      phase1b$avg_source_est[, 1L] - phase1b$avg_target_est
    )

    is_primary_cutoff <- abs(cutoff - primary_cutoff) <= 1e-12
    face$primary_cutoff_identity_checked <- is_primary_cutoff
    face$primary_cutoff_identity_exact <- if (is_primary_cutoff) {
      identical(as.numeric(fit$estimate), as.numeric(primary_direct$estimate)) &&
        identical(as.numeric(fit$se), as.numeric(primary_direct$se))
    } else {
      NA
    }
    if (is_primary_cutoff && !isTRUE(face$primary_cutoff_identity_exact)) {
      stop("primary-cutoff reaggregation failed the exact identity check.",
           call. = FALSE)
    }

    target <- make_result_row(
      fit$target_only, "target_only_ate", tate_truth, rho, cutoff
    )
    diagnostic_names <- setdiff(names(face), names(target))
    for (name in diagnostic_names) target[[name]] <- NA
    for (name in names(cell_diagnostics)) {
      target[[name]] <- cell_diagnostics[[name]]
    }
    target$target_anchor_weight <- 1
    target$mean_abs_source_weight <- 0
    target$max_abs_source_weight <- 0
    target$primary_cutoff_identity_checked <- is_primary_cutoff
    target$primary_cutoff_identity_exact <- if (is_primary_cutoff) TRUE else NA

    row_index <- row_index + 1L
    rows[[row_index]] <- face
    row_index <- row_index + 1L
    rows[[row_index]] <- target
  }
}

result <- RoCE:::.bind_sim_result_list(rows)
result$allocated_cores <- resource_plan$allocated_cores
result$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
result$source_workers_per_arm <- resource_plan$source_workers
result$parallel_treatment_arms <- resource_plan$parallel_treatment_arms
result$fully_parallel_cores <- resource_plan$fully_parallel_cores
result$rho_group_positive_rho_workers <- grouped$reuse$positive_rho_workers
result$rho_group_positive_rho_backend <- grouped$reuse$positive_rho_backend
result$n_nuisance_refits_for_cutoff_grid <- 0L
result$task_elapsed_seconds <- proc.time()[["elapsed"]] - started_at
result$package_library <- package_provenance$library
result$package_fingerprint <- package_provenance$fingerprint
result$workflow_fingerprint <- workflow_fingerprint
result$manifest_fingerprint <- manifest_fingerprint
scheduler_provenance <- roce_scheduler_provenance()
for (field in names(scheduler_provenance)) {
  result[[field]] <- scheduler_provenance[[field]]
}
roce_commit_csv_bundle(list(result), output_path)
message("[done] wrote ", output_path)
}

main()
