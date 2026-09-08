#!/usr/bin/env Rscript

.weight_calibration_files_fingerprint <- function(files, sha256_file) {
  if (length(files) < 1L || anyNA(files) || any(!file.exists(files)) ||
      any(dir.exists(files))) {
    stop("calibration workflow fingerprint inputs are incomplete.",
         call. = FALSE)
  }
  hashes <- vapply(files, sha256_file, character(1L))
  digest_input <- tempfile("weight_calibration_hashes_")
  on.exit(unlink(digest_input, force = TRUE), add = TRUE)
  writeLines(hashes, digest_input, useBytes = TRUE)
  sha256_file(digest_input)
}

.weight_calibration_installed_package_fingerprint <- function(
    package_library, sha256_file) {
  package_root <- file.path(package_library, "RoCE")
  .weight_calibration_files_fingerprint(c(
    file.path(package_root, "DESCRIPTION"),
    file.path(package_root, "libs", "RoCE.so"),
    file.path(package_root, "R", "RoCE.rdb"),
    file.path(package_root, "R", "RoCE.rdx")
  ), sha256_file)
}

.weight_calibration_integer <- function(value, name, lower, upper) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (length(numeric_value) != 1L || is.na(numeric_value) ||
      !is.finite(numeric_value) || numeric_value != floor(numeric_value) ||
      numeric_value < lower || numeric_value > upper) {
    stop(name, " must be one integer in [", lower, ", ", upper, "].",
         call. = FALSE)
  }
  as.integer(numeric_value)
}

.weight_calibration_scalar_equal <- function(x, expected, name) {
  if (length(x) != 1L || is.na(x) ||
      !identical(as.character(x), as.character(expected))) {
    stop("locked calibration mismatch for '", name, "'.", call. = FALSE)
  }
  invisible(TRUE)
}

.weight_calibration_validate_task <- function(
    group_manifest, primary_manifest, group_task_id) {
  group_required <- c(
    "group_task_id", "experiment", "sim_id", "config", "p", "K",
    "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
    "M_tau", "M_tau_inference", "methods", "rho_values",
    "primary_task_ids"
  )
  primary_required <- c(
    "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
    "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
    "M_tau", "M_tau_inference", "methods"
  )
  missing_group <- setdiff(group_required, names(group_manifest))
  missing_primary <- setdiff(primary_required, names(primary_manifest))
  if (length(missing_group) > 0L || length(missing_primary) > 0L) {
    stop(
      "calibration manifests are missing required columns: ",
      paste(c(missing_group, missing_primary), collapse = ", "),
      call. = FALSE
    )
  }
  manifest_group_ids <- suppressWarnings(as.numeric(
    group_manifest$group_task_id
  ))
  if (anyNA(manifest_group_ids) || any(!is.finite(manifest_group_ids)) ||
      any(manifest_group_ids < 1) ||
      any(manifest_group_ids != floor(manifest_group_ids)) ||
      anyDuplicated(manifest_group_ids)) {
    stop("group manifest task IDs must be unique positive integers.",
         call. = FALSE)
  }
  if (group_task_id > nrow(group_manifest)) {
    stop("GROUP_TASK_ID is outside the group manifest.", call. = FALSE)
  }
  group <- group_manifest[group_task_id, , drop = FALSE]
  observed_group_id <- suppressWarnings(as.numeric(group$group_task_id))
  if (length(observed_group_id) != 1L || is.na(observed_group_id) ||
      !is.finite(observed_group_id) ||
      observed_group_id != floor(observed_group_id) ||
      observed_group_id != group_task_id) {
    stop("GROUP_TASK_ID does not equal the indexed group_manifest row.",
         call. = FALSE)
  }
  rho_values <- suppressWarnings(as.numeric(strsplit(
    as.character(group$rho_values), ";", fixed = TRUE
  )[[1L]]))
  primary_task_ids_numeric <- suppressWarnings(as.numeric(strsplit(
    as.character(group$primary_task_ids), ";", fixed = TRUE
  )[[1L]]))
  expected_rhos <- c(0, 0.5, 1, 1.5, 2, 2.5)
  if (!identical(rho_values, expected_rhos) ||
      length(primary_task_ids_numeric) != 6L ||
      anyNA(primary_task_ids_numeric) ||
      any(!is.finite(primary_task_ids_numeric)) ||
      any(primary_task_ids_numeric < 1) ||
      any(primary_task_ids_numeric != floor(primary_task_ids_numeric)) ||
      anyDuplicated(primary_task_ids_numeric)) {
    stop("group manifest must contain the exact six-rho task mapping.",
         call. = FALSE)
  }
  primary_task_ids <- as.integer(primary_task_ids_numeric)
  manifest_task_ids <- suppressWarnings(as.numeric(primary_manifest$task_id))
  if (anyNA(manifest_task_ids) || any(!is.finite(manifest_task_ids)) ||
      any(manifest_task_ids < 1) ||
      any(manifest_task_ids != floor(manifest_task_ids)) ||
      anyDuplicated(manifest_task_ids)) {
    stop("primary manifest task IDs must be unique positive integers.",
         call. = FALSE)
  }
  primary_index <- match(
    primary_task_ids, manifest_task_ids
  )
  if (anyNA(primary_index)) {
    stop("group task IDs are absent from the primary manifest.",
         call. = FALSE)
  }
  primary_rows <- primary_manifest[primary_index, , drop = FALSE]
  if (!identical(as.integer(primary_rows$task_id), primary_task_ids) ||
      !identical(as.numeric(primary_rows$rho), expected_rhos)) {
    stop("group task IDs do not map exactly to the primary manifest rhos.",
         call. = FALSE)
  }
  stable_columns <- setdiff(group_required, c(
    "group_task_id", "rho_values", "primary_task_ids"
  ))
  for (column in stable_columns) {
    group_value <- as.character(group[[column]][[1L]])
    primary_values <- as.character(primary_rows[[column]])
    if (length(unique(primary_values)) != 1L ||
        !identical(primary_values[[1L]], group_value)) {
      stop("group/primary manifest mismatch for '", column, "'.",
           call. = FALSE)
    }
  }

  locked <- list(
    experiment = "negative_transfer",
    config = "C1", p = 100L, K = 2L, cutoff = 1,
    n_site = 1000L, n_folds = 5L, nlambda_init = 100L,
    n_bootstrap = 5000L, M_tau = 5, M_tau_inference = 5
  )
  for (name in names(locked)) {
    .weight_calibration_scalar_equal(group[[name]], locked[[name]], name)
  }
  sim_id <- .weight_calibration_integer(group$sim_id, "sim_id", 1L, 10L)
  methods <- trimws(strsplit(
    as.character(group$methods), ",", fixed = TRUE
  )[[1L]])
  methods <- methods[nzchar(methods)]
  expected_methods <- c(
    "one_round_crossfit", "target_only", "sample_size",
    "inverse_variance", "federated_dr", "pooled_dr"
  )
  if (anyDuplicated(methods) || !setequal(methods, expected_methods)) {
    stop("calibration requires the exact locked comparison-method set.",
         call. = FALSE)
  }
  list(
    group = group,
    primary_rows = primary_rows,
    sim_id = sim_id,
    methods = methods,
    rho_values = expected_rhos,
    primary_task_ids = primary_task_ids
  )
}

.weight_calibration_validate_results <- function(
    results, task, B, expected_package, expected_workflow,
    expected_manifest, expected_group_manifest) {
  if (!is.data.frame(results) || nrow(results) < 1L) {
    stop("calibration produced no result rows.", call. = FALSE)
  }
  core <- c(
    "sim_id", "method", "estimate", "se", "bias", "coverage", "truth",
    "ci_lower", "ci_upper", "rho", "n_weight_bootstrap",
    "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
    "rho_group_manifest_fingerprint"
  )
  missing <- setdiff(core, names(results))
  if (length(missing) > 0L) {
    stop("calibration result schema is missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_core <- c(
    "sim_id", "estimate", "se", "bias", "truth", "ci_lower", "ci_upper",
    "rho", "n_weight_bootstrap"
  )
  if (any(vapply(results[numeric_core], function(x) {
    !is.numeric(x) || any(!is.finite(x))
  }, logical(1L))) || any(results$se <= 0) ||
      any(results$ci_lower > results$ci_upper) ||
      any(abs(results$bias - (results$estimate - results$truth)) > 1e-10)) {
    stop("calibration point-estimate schema contains invalid values.",
         call. = FALSE)
  }
  if (!identical(sort(unique(as.numeric(results$rho))), task$rho_values) ||
      any(results$sim_id != task$sim_id) ||
      any(results$n_weight_bootstrap != B)) {
    stop("calibration result identifiers do not match the requested task.",
         call. = FALSE)
  }
  for (rho in task$rho_values) {
    methods <- as.character(results$method[results$rho == rho])
    expected_methods <- c(
      task$methods,
      paste0(task$methods, "_ate"),
      "one_round_crossfit_ate_armwise",
      "one_round_crossfit_ate_hard_threshold"
    )
    if (anyDuplicated(methods) || !setequal(methods, expected_methods) ||
        length(methods) != length(expected_methods)) {
      stop("each rho must contain the exact locked result-method schema.",
           call. = FALSE)
    }
  }
  bootstrap_numeric <- c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
    "weight_bootstrap_failures"
  )
  bootstrap_other <- c(
    "weight_bootstrap_multiplier", "weight_bootstrap_relearn_weights",
    "weight_bootstrap_screening_rule"
  )
  missing_bootstrap <- setdiff(
    c(bootstrap_numeric, bootstrap_other), names(results)
  )
  if (length(missing_bootstrap) > 0L) {
    stop("calibration bootstrap schema is missing: ",
         paste(missing_bootstrap, collapse = ", "), call. = FALSE)
  }
  soft <- results$method == "one_round_crossfit_ate"
  nonsoft <- !soft
  if (sum(soft) != 6L || any(vapply(results[soft, bootstrap_numeric],
                                    function(x) any(!is.finite(x)),
                                    logical(1L))) ||
      any(results$weight_relearn_n_bootstrap[soft] != B) ||
      any(results$weight_bootstrap_failures[soft] != 0L) ||
      !is.logical(results$weight_bootstrap_relearn_weights) ||
      any(!results$weight_bootstrap_relearn_weights[soft]) ||
      any(results$weight_bootstrap_screening_rule[soft] != "soft_penalty") ||
      any(!is.na(unlist(
        results[nonsoft, c(bootstrap_numeric, bootstrap_other), drop = FALSE],
        use.names = FALSE
      )))) {
    stop("weight-relearning bootstrap rows failed the strict schema audit.",
         call. = FALSE)
  }
  expected <- c(
    package_fingerprint = expected_package,
    workflow_fingerprint = expected_workflow,
    manifest_fingerprint = expected_manifest,
    rho_group_manifest_fingerprint = expected_group_manifest
  )
  for (name in names(expected)) {
    values <- unique(tolower(trimws(as.character(results[[name]]))))
    if (length(values) != 1L || !identical(values, expected[[name]])) {
      stop("calibration result provenance mismatch for '", name, "'.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

run_weight_bootstrap_calibration_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) {
    stop(paste(
      "usage: run_weight_bootstrap_calibration.R",
      "GROUP_MANIFEST PRIMARY_MANIFEST GROUP_TASK_ID OUTPUT_DIR B"
    ), call. = FALSE)
  }
  project_root <- normalizePath(
    Sys.getenv("ROCE_PROJECT_ROOT", getwd()), mustWork = TRUE
  )
  group_manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
  primary_manifest_path <- normalizePath(args[[2L]], mustWork = TRUE)
  group_task_id <- .weight_calibration_integer(
    args[[3L]], "GROUP_TASK_ID", 1L, 10L
  )
  B <- .weight_calibration_integer(args[[5L]], "B", 200L, 5000L)
  output_directory <- file.path(
    normalizePath(dirname(args[[4L]]), mustWork = FALSE),
    basename(args[[4L]])
  )
  if (!nzchar(basename(output_directory)) ||
      basename(output_directory) %in% c(".", "..")) {
    stop("OUTPUT_DIR must be one concrete directory path.", call. = FALSE)
  }

  project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
  if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  source(file.path(project_root, "scripts", "slurm", "atomic_output.R"))
  source(file.path(
    project_root, "scripts", "slurm", "direct_tate_task_helpers.R"
  ))
  source(file.path(project_root, "scripts", "slurm", "result_provenance.R"))
  source(file.path(project_root, "scripts", "slurm", "resource_topology.R"))
  source(file.path(
    project_root, "scripts", "slurm", "simulation_qc_policy.R"
  ))

  workflow_files <- file.path(project_root, c(
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.sh",
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R",
    "scripts/slurm/atomic_output.R",
    "scripts/slurm/direct_tate_task_helpers.R",
    "scripts/slurm/result_provenance.R",
    "scripts/slurm/resource_topology.R",
    "scripts/slurm/simulation_qc_policy.R",
    "scripts/slurm/package_library_utils.sh"
  ))
  package_provenance <- roce_runtime_package_provenance(project_library)
  installed_package_fingerprint <-
    .weight_calibration_installed_package_fingerprint(
      package_provenance$library, roce_sha256_file
    )
  if (!identical(installed_package_fingerprint,
                 package_provenance$fingerprint)) {
    stop("installed package bytes do not match ROCE_PACKAGE_FINGERPRINT.",
         call. = FALSE)
  }
  expected_workflow <- roce_sha256_environment("ROCE_WORKFLOW_FINGERPRINT")
  expected_manifest <- roce_sha256_environment("ROCE_MANIFEST_FINGERPRINT")
  expected_group_manifest <- roce_sha256_environment(
    "ROCE_GROUP_MANIFEST_FINGERPRINT"
  )
  observed_workflow <- .weight_calibration_files_fingerprint(
    workflow_files, roce_sha256_file
  )
  observed_manifest <- roce_sha256_file(primary_manifest_path)
  observed_group_manifest <- roce_sha256_file(group_manifest_path)
  if (!identical(observed_workflow, expected_workflow) ||
      !identical(observed_manifest, expected_manifest) ||
      !identical(observed_group_manifest, expected_group_manifest)) {
    stop("calibration input/workflow fingerprint mismatch.", call. = FALSE)
  }

  task <- .weight_calibration_validate_task(
    read.csv(group_manifest_path, stringsAsFactors = FALSE),
    read.csv(primary_manifest_path, stringsAsFactors = FALSE),
    group_task_id
  )
  allocated_cores <- roce_read_positive_integer_env(
    "SLURM_CPUS_PER_TASK", 1L
  )
  nuisance_cv_threads <- roce_read_positive_integer_env(
    "ROCE_NUISANCE_CV_THREADS", 5L, 5L
  )
  if (nuisance_cv_threads != 5L) {
    stop("calibration requires exactly five nuisance-CV threads.",
         call. = FALSE)
  }
  resource_plan <- roce_slurm_resource_plan(
    source_count = 2L, n_folds = 5L,
    allocated_cores = allocated_cores,
    nuisance_cv_threads = nuisance_cv_threads
  )
  cores_per_positive_rho <- nuisance_cv_threads *
    if (resource_plan$parallel_treatment_arms) 2L else 1L
  positive_rho_workers <- min(
    5L, max(1L, allocated_cores %/% cores_per_positive_rho)
  )
  if (dir.exists(output_directory) || file.exists(output_directory)) {
    stop("calibration output already exists: ", output_directory,
         call. = FALSE)
  }
  lock_path <- paste0(output_directory, ".lock")
  release_lock <- roce_claim_task_lock(lock_path, c(
    paste0("group_task_id=", group_task_id),
    paste0("sim_id=", task$sim_id),
    paste0("B=", B),
    paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
    paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  ))
  on.exit(release_lock(), add = TRUE)

  simulation_args <- list(
    sim_id = task$sim_id,
    n_total = 3000L,
    K = 2L,
    p = 100L,
    config = "C1",
    methods = task$methods,
    verbose = FALSE,
    n_cores_internal = resource_plan$source_workers,
    nlambda_init = 100L,
    nuisance_lambda_rule = "min",
    estimand_type = "superpopulation",
    outcome_type = "binary",
    n_folds = 5L,
    aggregation_lambda = 1,
    n_bootstrap = 5000L,
    n_weight_bootstrap = B,
    M_tau = 5,
    M_tau_inference = 5,
    estimate_ate = TRUE,
    parallel_treatment_arms = resource_plan$parallel_treatment_arms,
    include_hard_threshold_diagnostic = TRUE,
    dgp_type = "face"
  )
  message(sprintf(
    paste0(
      "[weight-calibration] sim=%d config=C1 p=100 K=2 rhos=%s B=%d ",
      "allocated=%d cv_threads=%d source_workers=%d ",
      "parallel_arms=%s positive_rho_workers=%d package=%s fingerprint=%s"
    ),
    task$sim_id, paste(task$rho_values, collapse = ","), B,
    resource_plan$allocated_cores, resource_plan$nuisance_cv_threads,
    resource_plan$source_workers,
    as.character(resource_plan$parallel_treatment_arms),
    positive_rho_workers, package_provenance$library,
    package_provenance$fingerprint
  ))
  started_at <- proc.time()[["elapsed"]]
  grouped <- RoCE:::.run_face_rho_group(
    simulation_args = simulation_args,
    rho_values = task$rho_values,
    changed_sources = "s1",
    artifact_rhos = task$rho_values,
    positive_rho_workers = positive_rho_workers
  )
  elapsed <- proc.time()[["elapsed"]] - started_at
  if (!identical(names(grouped$artifacts), vapply(
    task$rho_values, format, character(1L), scientific = FALSE, trim = TRUE
  )) || any(vapply(grouped$artifacts, function(x) {
    is.null(x$data_split) ||
      is.null(x$direct_tate_results$one_round_crossfit)
  }, logical(1L)))) {
    stop("calibration did not retain complete artifacts for all six rhos.",
         call. = FALSE)
  }
  for (rho in task$rho_values) {
    rho_key <- format(rho, scientific = FALSE, trim = TRUE)
    fit <- grouped$artifacts[[rho_key]]$direct_tate_results$one_round_crossfit
    RoCE:::.validate_weight_bootstrap_result(fit)
  }

  rows <- RoCE:::.bind_sim_result_list(grouped$results)
  rows$rho_group_task_id <- group_task_id
  rows$task_id <- task$primary_task_ids[match(rows$rho, task$rho_values)]
  rows$rho_group_size <- 6L
  rows$rho_group_manifest_fingerprint <- expected_group_manifest
  rows$allocated_cores <- resource_plan$allocated_cores
  rows$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
  rows$source_workers_per_arm <- resource_plan$source_workers
  rows$parallel_treatment_arms <- resource_plan$parallel_treatment_arms
  rows$fully_parallel_cores <- resource_plan$fully_parallel_cores
  rows$rho_group_positive_rho_workers <- positive_rho_workers
  rows$rho_group_positive_rho_backend <- grouped$reuse$positive_rho_backend
  rows$task_elapsed_seconds <- elapsed
  rows$package_library <- package_provenance$library
  rows$package_fingerprint <- package_provenance$fingerprint
  rows$workflow_fingerprint <- expected_workflow
  rows$manifest_fingerprint <- expected_manifest
  scheduler <- roce_scheduler_provenance()
  for (name in names(scheduler)) rows[[name]] <- scheduler[[name]]

  for (rho in task$rho_values) {
    rho_key <- format(rho, scientific = FALSE, trim = TRUE)
    fit <- grouped$artifacts[[rho_key]]$direct_tate_results$one_round_crossfit
    soft <- rows[
      rows$rho == rho & rows$method == "one_round_crossfit_ate",
      , drop = FALSE
    ]
    if (nrow(soft) != 1L ||
        !isTRUE(all.equal(
          as.numeric(soft$estimate), as.numeric(fit$estimate),
          tolerance = 1e-12, check.attributes = FALSE
        )) ||
        !isTRUE(all.equal(
          as.numeric(soft$se), as.numeric(fit$se),
          tolerance = 1e-12, check.attributes = FALSE
        ))) {
      stop("saved artifact/result identity failed for rho=", rho, ".",
           call. = FALSE)
    }
  }

  .weight_calibration_validate_results(
    rows, task, B, package_provenance$fingerprint, expected_workflow,
    expected_manifest, expected_group_manifest
  )
  qc <- RoCE:::diagnose_simulation_results(
    rows, expected_replications = 1L
  )
  classification <- roce_classify_simulation_qc(qc)
  qc$implementation_failed <- classification$implementation_failed
  qc$statistical_review_required <-
    classification$statistical_review_required
  if (any(qc$implementation_failed) ||
      any(grepl(
        "weight_relearn_bootstrap_failures_detected",
        qc$diagnostic_status, fixed = TRUE
      ))) {
    stop("calibration diagnostic QC failed closed.", call. = FALSE)
  }

  writer <- function(staging_directory) {
    writer_package <- roce_runtime_package_provenance(project_library)
    writer_package_fingerprint <-
      .weight_calibration_installed_package_fingerprint(
        writer_package$library, roce_sha256_file
      )
    writer_workflow <- .weight_calibration_files_fingerprint(
      workflow_files, roce_sha256_file
    )
    writer_manifest <- roce_sha256_file(primary_manifest_path)
    writer_group_manifest <- roce_sha256_file(group_manifest_path)
    if (!identical(writer_package_fingerprint,
                   writer_package$fingerprint) ||
        !identical(writer_package$fingerprint,
                   package_provenance$fingerprint) ||
        !identical(writer_workflow, expected_workflow) ||
        !identical(writer_manifest, expected_manifest) ||
        !identical(writer_group_manifest, expected_group_manifest)) {
      stop("writer provenance changed before atomic publication.",
           call. = FALSE)
    }
    results_path <- file.path(staging_directory, "results.csv")
    qc_path <- file.path(staging_directory, "diagnostic_qc.csv")
    artifacts_path <- file.path(staging_directory, "artifacts.rds")
    metadata_path <- file.path(staging_directory, "metadata.txt")
    utils::write.csv(rows, results_path, row.names = FALSE, na = "")
    utils::write.csv(qc, qc_path, row.names = FALSE, na = "")
    saveRDS(list(
      group_result = grouped,
      task = task,
      simulation_args = simulation_args,
      resource_plan = resource_plan,
      positive_rho_workers = positive_rho_workers
    ), artifacts_path)
    writeLines(c(
      "weight_bootstrap_calibration=complete",
      paste0("group_task_id=", group_task_id),
      paste0("sim_id=", task$sim_id),
      paste0("primary_task_ids=", paste(task$primary_task_ids, collapse = ";")),
      paste0("rho_values=", paste(task$rho_values, collapse = ";")),
      paste0("n_weight_bootstrap=", B),
      paste0("package_library=", writer_package$library),
      paste0("package_fingerprint=", writer_package$fingerprint),
      paste0("workflow_fingerprint=", writer_workflow),
      paste0("manifest_fingerprint=", writer_manifest),
      paste0("group_manifest_fingerprint=", writer_group_manifest),
      paste0("elapsed_seconds=", format(elapsed, digits = 17)),
      paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    ), metadata_path, useBytes = TRUE)
    payloads <- c(results_path, qc_path, artifacts_path, metadata_path)
    payload_hashes <- vapply(payloads, roce_sha256_file, character(1L))
    writeLines(
      paste(payload_hashes, basename(payloads), sep = "  "),
      file.path(staging_directory, "sha256.txt"), useBytes = TRUE
    )
  }
  roce_write_atomic_directory(
    output_directory, writer, caller = "weight bootstrap calibration"
  )
  message("[done] atomically wrote calibration: ", output_directory)
  invisible(output_directory)
}

if (sys.nframe() == 0L) {
  run_weight_bootstrap_calibration_main()
}
