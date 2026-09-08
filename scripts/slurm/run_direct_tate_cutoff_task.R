#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
if (length(args) != 2L) {
  stop("usage: run_direct_tate_cutoff_task.R MANIFEST.csv ARRAY_TASK_ID")
}

manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
row_id <- as.integer(args[[2L]])
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (is.na(row_id) || row_id < 1L || row_id > nrow(manifest)) {
  stop(sprintf(
    "array task ID %s is outside manifest rows 1:%d.",
    args[[2L]], nrow(manifest)
  ))
}
task <- manifest[row_id, , drop = FALSE]
task_started_at <- proc.time()[["elapsed"]]
task_dimension <- suppressWarnings(as.integer(task$p))
if (length(task_dimension) != 1L || is.na(task_dimension) ||
    task_dimension != 100L) {
  stop("the cutoff workflow is restricted to p=100; found p=", task$p)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
package_provenance <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment(
  "ROCE_WORKFLOW_FINGERPRINT"
)
manifest_fingerprint <- roce_sha256_environment(
  "ROCE_MANIFEST_FINGERPRINT"
)

output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT",
  file.path(dirname(manifest_path), "raw_cutoff_diagnostic")
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
  output_root, "cutoff_task_locks",
  sprintf("task_%06d.lock", as.integer(task$task_id))
)
release_lock <- roce_claim_task_lock(lock_path, c(
  paste0("task_id=", as.integer(task$task_id)),
  paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
  paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
))
on.exit(release_lock(), add = TRUE)

cutoffs <- as.numeric(strsplit(task$cutoffs, ",", fixed = TRUE)[[1L]])
if (length(cutoffs) == 0L || any(!is.finite(cutoffs)) || any(cutoffs <= 0)) {
  stop("manifest cutoffs must be a comma-separated list of positive numbers.")
}
primary_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "ROCE_PRIMARY_CUTOFF", as.character(1 / RoCE:::AGG_WALD_LAMBDA)
)))
manifest_primary_cutoff <- if ("primary_cutoff" %in% names(task)) {
  suppressWarnings(as.numeric(task$primary_cutoff))
} else {
  primary_cutoff
}
if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
    primary_cutoff <= 0 ||
    length(manifest_primary_cutoff) != 1L ||
    !is.finite(manifest_primary_cutoff) ||
    abs(primary_cutoff - manifest_primary_cutoff) > 1e-12 ||
    !any(abs(cutoffs - primary_cutoff) <= 1e-12)) {
  stop(
    paste0(
      "ROCE_PRIMARY_CUTOFF and the manifest primary cutoff must agree ",
      "on one positive cutoff in the grid."
    ),
       call. = FALSE)
}
M_tau <- if ("M_tau" %in% names(task)) as.numeric(task$M_tau) else 5
M_tau_inference <- if ("M_tau_inference" %in% names(task)) {
  as.numeric(task$M_tau_inference)
} else {
  5
}
nlambda_init <- if ("nlambda_init" %in% names(task)) {
  suppressWarnings(as.integer(task$nlambda_init))
} else {
  NA_integer_
}
if (!is.finite(M_tau) || M_tau <= 0 ||
    is.na(M_tau_inference) || M_tau_inference <= 0) {
  stop("manifest truncation radii must be positive (inference may be Inf).")
}
if (length(nlambda_init) != 1L || is.na(nlambda_init) ||
    nlambda_init < 2L) {
  stop("manifest nlambda_init must be one integer >= 2.")
}

allocated_cores <- roce_read_positive_integer_env(
  "SLURM_CPUS_PER_TASK", 1L
)
nuisance_cv_threads <- roce_read_positive_integer_env(
  "ROCE_NUISANCE_CV_THREADS", 1L, as.integer(task$n_folds)
)
resource_plan <- roce_slurm_resource_plan(
  source_count = as.integer(task$K),
  n_folds = as.integer(task$n_folds),
  allocated_cores = allocated_cores,
  nuisance_cv_threads = nuisance_cv_threads
)

message(sprintf(
  paste0(
    "[task %d] experiment=%s sim=%d C=%s p=%d K=%d rho=%g ",
    paste0(
      "cutoffs=%s M_fit=%g M_inf=%g nlambda=%d n_site=%d folds=%d ",
      "allocated=%d cv_threads=%d source_workers_per_arm=%d parallel_arms=%s"
    )
  ),
  task$task_id, task$experiment, task$sim_id, task$config,
  task$p, task$K, task$rho, paste(cutoffs, collapse = ","),
  M_tau, M_tau_inference, nlambda_init, task$n_site, task$n_folds,
  resource_plan$allocated_cores, resource_plan$nuisance_cv_threads,
  resource_plan$source_workers, resource_plan$parallel_treatment_arms
))

parallel_treatment_arms <- resource_plan$parallel_treatment_arms
source_cores <- resource_plan$source_workers

set.seed(as.integer(task$sim_id))
generated <- generate_simulation_data(
  n_total = as.integer(task$n_site) * (as.integer(task$K) + 1L),
  K = as.integer(task$K),
  p = as.integer(task$p),
  config = task$config,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  dgp_type = "face",
  ate_deviation = as.numeric(task$rho),
  n_deviated_sites = if (as.numeric(task$rho) > 0) 1L else 0L,
  warn_ignored = FALSE
)
data_split <- split_data_by_site(generated)
cell_diagnostics <- RoCE:::.binary_cell_diagnostics(
  data_split, "binomial"
)

# The arm-specific nuisance fits do not depend on the aggregation cutoff.
# Fit them once at the selected primary cutoff, retain their fold intermediates, and
# recompute only the TATE aggregation for the remaining cutoffs.
reference_fit <- run_tate_crossfit(
  data_split = data_split,
  n_folds = as.integer(task$n_folds),
  communication_mode = "one_round",
  lambda_selection = 1 / primary_cutoff,
  verbose = FALSE,
  n_cores = source_cores,
  nlambda_init = nlambda_init,
  family = "binomial",
  M_tau = M_tau,
  M_tau_inference = M_tau_inference,
  parallel_arms = parallel_treatment_arms
)
reference_diagnostics <- c(
  RoCE:::.summarize_tate_crossfit_timing(
    reference_fit$arm_results$mu1,
    reference_fit$arm_results$mu0
  ),
  RoCE:::.summarize_tate_nuisance_diagnostics(
    reference_fit$arm_results$mu1,
    reference_fit$arm_results$mu0
  )
)

truth <- generated$mu1_true - generated$mu0_true
n_total <- sum(vapply(data_split, function(site) site$n, integer(1L)))
make_row <- function(method, fit, cutoff) {
  estimate <- fit$estimate
  se <- fit$se
  z_alpha_05 <- stats::qnorm(0.975)
  ci_lower <- estimate - z_alpha_05 * se
  ci_upper <- estimate + z_alpha_05 * se
  data.frame(
    sim_id = as.integer(task$sim_id),
    method = method,
    estimate = estimate,
    se = se,
    truth = truth,
    bias = estimate - truth,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    coverage = truth >= ci_lower & truth <= ci_upper,
    ci_width = ci_upper - ci_lower,
    n_total = n_total,
    K = as.integer(task$K),
    p = as.integer(task$p),
    config = task$config,
    heterogeneity_type = RoCE:::.face_heterogeneity_type(
      ate_deviation = as.numeric(task$rho),
      n_deviated_sites = as.integer(as.numeric(task$rho) > 0)
    ),
    estimand_type = "superpopulation",
    estimand_scope = "tate",
    dgp_type = "face",
    outcome_family = "binomial",
    aggregation_lambda = 1 / cutoff,
    aggregation_cutoff = cutoff,
    nlambda_init = nlambda_init,
    # No comparison estimator is fitted in this reuse-only diagnostic. Keep
    # the common result schema explicit so the setting-level QC pipeline can
    # group and audit cutoff results without inventing bootstrap metadata.
    n_bootstrap = NA_integer_,
    M_tau = M_tau,
    M_tau_inference = M_tau_inference,
    primary_cutoff = primary_cutoff,
    primary_cutoff_identity_checked =
      abs(cutoff - primary_cutoff) <= 1e-12,
    primary_cutoff_identity_exact =
      abs(cutoff - primary_cutoff) <= 1e-12,
    task_id = as.integer(task$task_id),
    experiment = task$experiment,
    rho = as.numeric(task$rho),
    cutoff = cutoff,
    n_site = as.integer(task$n_site),
    n_folds = as.integer(task$n_folds),
    stringsAsFactors = FALSE
  )
}

rows <- vector("list", 2L * length(cutoffs))
for (index in seq_along(cutoffs)) {
  cutoff <- cutoffs[[index]]
  lambda <- 1 / cutoff
  fit <- if (abs(cutoff - primary_cutoff) <= 1e-12) {
    reference_fit
  } else {
    calculate_tate_crossfit_aggregation(
      data_split = data_split,
      mu1_result = reference_fit$arm_results$mu1,
      mu0_result = reference_fit$arm_results$mu0,
      lambda_selection = lambda,
      verbose = FALSE
    )
  }

  face_row <- make_row("one_round_crossfit_ate", fit, cutoff)
  aggregation_diagnostics <-
    RoCE:::.summarize_tate_aggregation_diagnostics(fit)
  for (diagnostic_name in names(aggregation_diagnostics)) {
    face_row[[diagnostic_name]] <- aggregation_diagnostics[[diagnostic_name]]
  }

  target_row <- make_row(
    "target_only_ate",
    list(
      estimate = fit$target_only$estimate,
      se = fit$target_only$se
    ),
    cutoff
  )
  for (diagnostic_name in names(aggregation_diagnostics)) {
    target_row[[diagnostic_name]] <- NA_real_
  }
  target_row$target_anchor_weight <- 1
  target_row$mean_abs_source_weight <- 0
  target_row$max_abs_source_weight <- 0
  for (diagnostic_name in names(reference_diagnostics)) {
    face_row[[diagnostic_name]] <- reference_diagnostics[[diagnostic_name]]
    target_row[[diagnostic_name]] <- NA_real_
  }
  for (diagnostic_name in names(cell_diagnostics)) {
    face_row[[diagnostic_name]] <- cell_diagnostics[[diagnostic_name]]
    target_row[[diagnostic_name]] <- cell_diagnostics[[diagnostic_name]]
  }

  rows[[2L * index - 1L]] <- face_row
  rows[[2L * index]] <- target_row
}
result <- do.call(rbind, rows)
result$target_anchor_nlambda <- RoCE:::LAMBDA_GRID_SIZE_STANDARD
result$allocated_cores <- resource_plan$allocated_cores
result$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
result$source_workers_per_arm <- resource_plan$source_workers
result$parallel_treatment_arms <- resource_plan$parallel_treatment_arms
result$fully_parallel_cores <- resource_plan$fully_parallel_cores
result$task_elapsed_seconds <- proc.time()[["elapsed"]] - task_started_at
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
