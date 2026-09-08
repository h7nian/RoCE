#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
if (length(args) != 3L) {
  stop(
    paste(
      "usage: run_direct_tate_rho_group_task.R",
      "GROUP_MANIFEST.csv PRIMARY_MANIFEST.csv ARRAY_TASK_ID"
    ),
    call. = FALSE
  )
}

group_manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
primary_manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
group_task_id <- as.integer(args[[3L]])
group_manifest <- read.csv(group_manifest_path, stringsAsFactors = FALSE)
primary_manifest <- read.csv(primary_manifest_path, stringsAsFactors = FALSE)
if (is.na(group_task_id) || group_task_id < 1L ||
    group_task_id > nrow(group_manifest)) {
  stop("array task ID is outside the rho-group manifest.", call. = FALSE)
}
group <- group_manifest[group_task_id, , drop = FALSE]
if (as.integer(group$p) != 100L ||
    as.integer(group$group_task_id) != group_task_id) {
  stop("rho-group production is restricted to correctly indexed p=100 rows.",
       call. = FALSE)
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

package_provenance <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment(
  "ROCE_WORKFLOW_FINGERPRINT"
)
manifest_fingerprint <- roce_sha256_environment(
  "ROCE_MANIFEST_FINGERPRINT"
)
group_manifest_fingerprint <- roce_sha256_environment(
  "ROCE_GROUP_MANIFEST_FINGERPRINT"
)

rho_values <- as.numeric(strsplit(group$rho_values, ";", fixed = TRUE)[[1L]])
primary_task_ids <- as.integer(
  strsplit(group$primary_task_ids, ";", fixed = TRUE)[[1L]]
)
expected_rhos <- c(0, 0.5, 1, 1.5, 2, 2.5)
if (!identical(rho_values, expected_rhos) ||
    length(primary_task_ids) != length(expected_rhos) ||
    anyNA(primary_task_ids) || anyDuplicated(primary_task_ids)) {
  stop("rho-group row does not contain the exact six-rho task mapping.",
       call. = FALSE)
}
primary_rows <- primary_manifest[
  match(primary_task_ids, as.integer(primary_manifest$task_id)), , drop = FALSE
]
if (nrow(primary_rows) != length(expected_rhos) ||
    anyNA(primary_rows$task_id) ||
    !identical(as.integer(primary_rows$task_id), primary_task_ids) ||
    !identical(as.numeric(primary_rows$rho), expected_rhos)) {
  stop("rho-group task IDs do not map exactly to the primary manifest.",
       call. = FALSE)
}
stable_columns <- c(
  "experiment", "sim_id", "config", "p", "K", "cutoff", "n_site",
  "n_folds", "nlambda_init", "n_bootstrap", "M_tau",
  "M_tau_inference", "methods", "deviation_mechanism"
)
for (column in stable_columns) {
  if (length(unique(primary_rows[[column]])) != 1L ||
      !identical(as.character(primary_rows[[column]][[1L]]),
                 as.character(group[[column]][[1L]]))) {
    stop("rho-group invariant mismatch for field '", column, "'.",
         call. = FALSE)
  }
}

output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT", file.path(dirname(primary_manifest_path), "raw")
)
sensitivity_output_root <- Sys.getenv(
  "ROCE_SENSITIVITY_OUTPUT_ROOT",
  file.path(dirname(output_root), "reused_sensitivity_raw")
)
commit_root <- file.path(output_root, "rho_group_commits")
primary_output_paths <- file.path(
  output_root, sprintf("task_%06d.csv", primary_task_ids)
)
sensitivity_indices <- which(vapply(
  seq_len(nrow(primary_rows)),
  function(index) roce_is_reused_sensitivity_task(
    primary_rows[index, , drop = FALSE]
  ),
  logical(1L)
))
sensitivity_output_paths <- file.path(
  sensitivity_output_root,
  sprintf("task_%06d.csv", primary_task_ids[sensitivity_indices])
)
commit_path <- file.path(
  commit_root, sprintf("group_%06d_committed.txt", group_task_id)
)
expected_paths <- c(primary_output_paths, sensitivity_output_paths)
if (file.exists(commit_path)) {
  if (!all(file.exists(expected_paths))) {
    stop("rho-group commit exists but one or more outputs are missing.",
         call. = FALSE)
  }
  message("[skip] committed rho group already exists: ", commit_path)
  quit(save = "no", status = 0L)
}
if (any(file.exists(expected_paths))) {
  stop(
    "partial rho-group outputs exist without a commit sentinel; refusing overwrite.",
    call. = FALSE
  )
}
lock_path <- file.path(
  output_root, "rho_group_locks",
  sprintf("group_%06d.lock", group_task_id)
)
release_lock <- roce_claim_task_lock(lock_path, c(
  paste0("group_task_id=", group_task_id),
  paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
  paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
))
on.exit(release_lock(), add = TRUE)

nlambda_init <- as.integer(group$nlambda_init)
n_bootstrap <- as.integer(group$n_bootstrap)
if (nlambda_init < 2L || n_bootstrap < 2L) {
  stop("rho-group nuisance grid and bootstrap count must both be >= 2.",
       call. = FALSE)
}
allocated_cores <- roce_read_positive_integer_env(
  "SLURM_CPUS_PER_TASK", 1L
)
nuisance_cv_threads <- roce_read_positive_integer_env(
  "ROCE_NUISANCE_CV_THREADS", 1L, as.integer(group$n_folds)
)
resource_plan <- roce_slurm_resource_plan(
  source_count = as.integer(group$K),
  n_folds = as.integer(group$n_folds),
  allocated_cores = allocated_cores,
  nuisance_cv_threads = nuisance_cv_threads
)
default_positive_rho_workers <- min(
  length(expected_rhos) - 1L,
  max(1L, resource_plan$allocated_cores %/%
        resource_plan$nuisance_cv_threads)
)
positive_rho_workers <- roce_read_positive_integer_env(
  "ROCE_POSITIVE_RHO_WORKERS", default_positive_rho_workers,
  length(expected_rhos) - 1L
)
task_verbose <- roce_parse_boolean(
  Sys.getenv("ROCE_TASK_VERBOSE", "false"), "ROCE_TASK_VERBOSE"
)
include_hard_threshold_diagnostic <- roce_parse_boolean(
  Sys.getenv("ROCE_INCLUDE_HARD_THRESHOLD_DIAGNOSTIC", "false"),
  "ROCE_INCLUDE_HARD_THRESHOLD_DIAGNOSTIC"
)
include_quadratic_bias_rule <- roce_parse_boolean(
  Sys.getenv("ROCE_INCLUDE_QUADRATIC_BIAS_RULE", "true"),
  "ROCE_INCLUDE_QUADRATIC_BIAS_RULE"
)
methods <- roce_parse_method_list(group$methods)
simulation_args <- list(
  sim_id = as.integer(group$sim_id),
  n_total = as.integer(group$n_site) * (as.integer(group$K) + 1L),
  K = as.integer(group$K),
  p = as.integer(group$p),
  config = as.character(group$config),
  methods = methods,
  verbose = task_verbose,
  n_cores_internal = resource_plan$source_workers,
  nlambda_init = nlambda_init,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  n_folds = as.integer(group$n_folds),
  aggregation_lambda = 1 / as.numeric(group$cutoff),
  n_bootstrap = n_bootstrap,
  M_tau = as.numeric(group$M_tau),
  M_tau_inference = as.numeric(group$M_tau_inference),
  estimate_ate = TRUE,
  parallel_treatment_arms = resource_plan$parallel_treatment_arms,
  include_hard_threshold_diagnostic =
    include_hard_threshold_diagnostic,
  include_quadratic_bias_rule = include_quadratic_bias_rule,
  dgp_type = "face",
  deviation_mechanism = roce_task_deviation_mechanism(group)
)

message(sprintf(
  paste0(
    "[rho-group %d] sim=%d C=%s p=%d K=%d rhos=%s nlambda=%d ",
    "boot=%d allocated=%d cv_threads=%d source_workers=%d ",
    "positive_rho_workers=%d"
  ),
  group_task_id, group$sim_id, group$config, group$p, group$K,
  paste(expected_rhos, collapse = ","), nlambda_init, n_bootstrap,
  resource_plan$allocated_cores, resource_plan$nuisance_cv_threads,
  resource_plan$source_workers, positive_rho_workers
))
group_started_at <- proc.time()[["elapsed"]]
group_result <- RoCE:::.run_face_rho_group(
  simulation_args = simulation_args,
  rho_values = expected_rhos,
  changed_sources = "s1",
  artifact_rhos = primary_rows$rho[sensitivity_indices],
  positive_rho_workers = positive_rho_workers
)
group_elapsed <- proc.time()[["elapsed"]] - group_started_at

annotated_results <- vector("list", length(expected_rhos))
annotated_sensitivities <- vector("list", length(sensitivity_indices))
for (index in seq_along(expected_rhos)) {
  task <- primary_rows[index, , drop = FALSE]
  annotated <- roce_annotate_direct_tate_rows(
    result = group_result$results[[index]],
    task = task,
    resource_plan = resource_plan,
    task_elapsed_seconds = group_elapsed,
    package_provenance = package_provenance,
    workflow_fingerprint = workflow_fingerprint,
    manifest_fingerprint = manifest_fingerprint,
    metadata_policy = "from_task"
  )
  annotated$rho_group_task_id <- group_task_id
  annotated$rho_group_size <- length(expected_rhos)
  annotated$rho_group_manifest_fingerprint <- group_manifest_fingerprint
  annotated$rho_group_positive_rho_workers <-
    group_result$reuse$positive_rho_workers
  annotated$rho_group_positive_rho_backend <-
    group_result$reuse$positive_rho_backend
  annotated$rho_group_source_workers_per_fit <- if (task$rho == 0) {
    resource_plan$source_workers
  } else {
    group_result$reuse$positive_rho_source_workers_per_fit
  }
  annotated_results[[index]] <- annotated
}
for (sidecar_index in seq_along(sensitivity_indices)) {
  index <- sensitivity_indices[[sidecar_index]]
  task <- primary_rows[index, , drop = FALSE]
  rho_key <- format(task$rho, scientific = FALSE, trim = TRUE)
  artifacts <- group_result$artifacts[[rho_key]]
  if (is.null(artifacts)) {
    stop("missing fitted artifact for required sensitivity rho=", task$rho,
         call. = FALSE)
  }
  sensitivity <- roce_build_reused_sensitivity_rows(
    task = task,
    artifacts = artifacts,
    sensitivity_grid = roce_default_sensitivity_grid(),
    nlambda_init = nlambda_init,
    n_bootstrap = n_bootstrap,
    outcome_family = "binomial"
  )
  roce_validate_sensitivity_identity(
    sensitivity,
    group_result$results[[index]]
  )
  sensitivity <- roce_annotate_direct_tate_rows(
    result = sensitivity,
    task = task,
    resource_plan = resource_plan,
    task_elapsed_seconds = group_elapsed,
    package_provenance = package_provenance,
    workflow_fingerprint = workflow_fingerprint,
    manifest_fingerprint = manifest_fingerprint,
    metadata_policy = "preserve"
  )
  sensitivity$rho_group_task_id <- group_task_id
  sensitivity$rho_group_size <- length(expected_rhos)
  sensitivity$rho_group_manifest_fingerprint <- group_manifest_fingerprint
  sensitivity$rho_group_positive_rho_workers <-
    group_result$reuse$positive_rho_workers
  sensitivity$rho_group_positive_rho_backend <-
    group_result$reuse$positive_rho_backend
  sensitivity$rho_group_source_workers_per_fit <- if (task$rho == 0) {
    resource_plan$source_workers
  } else {
    group_result$reuse$positive_rho_source_workers_per_fit
  }
  annotated_sensitivities[[sidecar_index]] <- sensitivity
}

commit_lines <- c(
  "rho_group_commit=complete",
  paste0("group_task_id=", group_task_id),
  paste0("primary_task_ids=", paste(primary_task_ids, collapse = ";")),
  paste0("package_fingerprint=", package_provenance$fingerprint),
  paste0("workflow_fingerprint=", workflow_fingerprint),
  paste0("manifest_fingerprint=", manifest_fingerprint),
  paste0("group_manifest_fingerprint=", group_manifest_fingerprint)
)
# Sidecars precede primaries, and the sentinel is the final commit marker.
roce_commit_csv_bundle(
  values = c(annotated_sensitivities, annotated_results),
  output_paths = c(sensitivity_output_paths, primary_output_paths),
  sentinel_path = commit_path,
  sentinel_lines = commit_lines
)
message("[done] committed rho group: ", commit_path)
}

main()
