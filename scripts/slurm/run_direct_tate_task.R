#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
if (length(args) != 2L) {
  stop("usage: run_direct_tate_task.R MANIFEST.csv ARRAY_TASK_ID")
}

manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
task_id <- as.integer(args[[2L]])
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (is.na(task_id) || task_id < 1L || task_id > nrow(manifest)) {
  stop(sprintf("array task ID %s is outside manifest rows 1:%d.",
               args[[2L]], nrow(manifest)))
}
task <- manifest[task_id, , drop = FALSE]
task_started_at <- proc.time()[["elapsed"]]
task_dimension <- suppressWarnings(as.integer(task$p))
if (length(task_dimension) != 1L || is.na(task_dimension) ||
    task_dimension != 100L) {
  stop(
    "the MSI TATE production workflow is restricted to p=100; found p=",
    task$p,
    ". Use a separate diagnostic manifest for any other dimension."
  )
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "direct_tate_sensitivity_rows.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
package_provenance <- roce_runtime_package_provenance(
  project_library
)
project_library <- package_provenance$library
package_fingerprint <- package_provenance$fingerprint
workflow_fingerprint <- roce_sha256_environment(
  "ROCE_WORKFLOW_FINGERPRINT"
)
manifest_fingerprint <- roce_sha256_environment(
  "ROCE_MANIFEST_FINGERPRINT"
)

output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT",
  file.path(dirname(manifest_path), "raw")
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(
  output_root,
  sprintf("task_%06d.csv", as.integer(task$task_id))
)
sensitivity_expected <- roce_is_reused_sensitivity_task(task)
sensitivity_output_root <- Sys.getenv(
  "ROCE_SENSITIVITY_OUTPUT_ROOT",
  file.path(dirname(output_root), "reused_sensitivity_raw")
)
sensitivity_output_path <- file.path(
  sensitivity_output_root,
  sprintf("task_%06d.csv", as.integer(task$task_id))
)
if (file.exists(output_path)) {
  if (sensitivity_expected && !file.exists(sensitivity_output_path)) {
    stop(
      "primary output exists but its required sensitivity sidecar is missing: ",
      sensitivity_output_path,
      call. = FALSE
    )
  }
  message("[skip] output already exists: ", output_path)
  quit(save = "no", status = 0L)
}
if (sensitivity_expected && file.exists(sensitivity_output_path)) {
  stop(
    "orphan sensitivity sidecar exists without its primary commit file: ",
    sensitivity_output_path,
    call. = FALSE
  )
}
lock_path <- file.path(
  output_root, "task_locks",
  sprintf("task_%06d.lock", as.integer(task$task_id))
)
release_lock <- roce_claim_task_lock(lock_path, c(
  paste0("task_id=", as.integer(task$task_id)),
  paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
  paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
))
on.exit(release_lock(), add = TRUE)

methods <- roce_parse_method_list(task$methods)
lambda <- 1 / as.numeric(task$cutoff)
M_tau <- if ("M_tau" %in% names(task)) as.numeric(task$M_tau) else 5
M_tau_inference <- if ("M_tau_inference" %in% names(task)) {
  as.numeric(task$M_tau_inference)
} else {
  5
}
n_bootstrap <- if ("n_bootstrap" %in% names(task)) {
  suppressWarnings(as.integer(task$n_bootstrap))
} else {
  5000L
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
if (length(n_bootstrap) != 1L || is.na(n_bootstrap) || n_bootstrap < 2L) {
  stop("manifest n_bootstrap must be one integer >= 2.")
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
parallel_treatment_arms <- resource_plan$parallel_treatment_arms
source_cores <- resource_plan$source_workers
task_verbose <- roce_parse_boolean(
  Sys.getenv("ROCE_TASK_VERBOSE", "false"), "ROCE_TASK_VERBOSE"
)

message(sprintf(
  paste0(
    "[task %d] experiment=%s sim=%d C=%s p=%d K=%d rho=%g ",
    paste0(
      "cutoff=%g M_fit=%g M_inf=%g nlambda=%d boot=%d n_site=%d folds=%d ",
      paste0(
        "allocated=%d cv_threads=%d source_workers_per_arm=%d ",
        "parallel_arms=%s verbose=%s methods=%s"
      )
    )
  ),
  task$task_id, task$experiment, task$sim_id, task$config,
  task$p, task$K, task$rho, task$cutoff, M_tau, M_tau_inference,
  nlambda_init, n_bootstrap,
  task$n_site, task$n_folds, resource_plan$allocated_cores,
  resource_plan$nuisance_cv_threads, resource_plan$source_workers,
  resource_plan$parallel_treatment_arms, task_verbose,
  paste(methods, collapse = "+")
))

simulation_args <- list(
  sim_id = as.integer(task$sim_id),
  n_total = as.integer(task$n_site) * (as.integer(task$K) + 1L),
  K = as.integer(task$K),
  p = as.integer(task$p),
  config = task$config,
  methods = methods,
  verbose = task_verbose,
  n_cores_internal = source_cores,
  nlambda_init = nlambda_init,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  n_folds = as.integer(task$n_folds),
  aggregation_lambda = lambda,
  n_bootstrap = n_bootstrap,
  M_tau = M_tau,
  M_tau_inference = M_tau_inference,
  estimate_ate = TRUE,
  parallel_treatment_arms = parallel_treatment_arms,
  dgp_type = "face",
  ate_deviation = as.numeric(task$rho),
  n_deviated_sites = if (as.numeric(task$rho) > 0) 1L else 0L
)
if (sensitivity_expected) {
  if (!exists(
    "reaggregate_tate_sensitivity_grid",
    envir = asNamespace("RoCE"),
    inherits = FALSE
  ) || !"return_fitted_tate" %in%
      names(formals(RoCE::run_single_simulation))) {
    stop(
      paste0(
        "the installed package does not support provenance-locked reused ",
        "sensitivities; install the audited final library before production."
      ),
      call. = FALSE
    )
  }
  simulation_args$return_fitted_tate <- TRUE
}
result <- do.call(run_single_simulation, simulation_args)
simulation_artifacts <- attr(result, "roce_simulation_artifacts")
attr(result, "roce_simulation_artifacts") <- NULL

sensitivity_result <- NULL
if (sensitivity_expected) {
  sensitivity_result <- roce_build_reused_sensitivity_rows(
    task = task,
    artifacts = simulation_artifacts,
    sensitivity_grid = roce_default_sensitivity_grid(),
    nlambda_init = nlambda_init,
    n_bootstrap = n_bootstrap,
    outcome_family = "binomial"
  )
  roce_validate_sensitivity_identity(sensitivity_result, result)
}

task_elapsed_seconds <- proc.time()[["elapsed"]] - task_started_at
result <- roce_annotate_direct_tate_rows(
  result = result,
  task = task,
  resource_plan = resource_plan,
  task_elapsed_seconds = task_elapsed_seconds,
  package_provenance = package_provenance,
  workflow_fingerprint = workflow_fingerprint,
  manifest_fingerprint = manifest_fingerprint,
  metadata_policy = "from_task"
)

if (!is.null(sensitivity_result)) {
  sensitivity_result <- roce_annotate_direct_tate_rows(
    result = sensitivity_result,
    task = task,
    resource_plan = resource_plan,
    task_elapsed_seconds = task_elapsed_seconds,
    package_provenance = package_provenance,
    workflow_fingerprint = workflow_fingerprint,
    manifest_fingerprint = manifest_fingerprint,
    metadata_policy = "preserve"
  )
}

if (is.null(sensitivity_result)) {
  roce_commit_csv_bundle(list(result), output_path)
} else {
  # Commit the primary file last so existing completion scanners never treat
  # a task as complete before its required sensitivity sidecar is durable.
  roce_commit_csv_bundle(
    list(sensitivity_result, result),
    c(sensitivity_output_path, output_path)
  )
}
message("[done] wrote ", output_path)
}

main()
