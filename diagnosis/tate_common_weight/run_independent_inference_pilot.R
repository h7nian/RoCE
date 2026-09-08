#!/usr/bin/env Rscript

.inference_pilot_rhos <- function() c(0, 0.5, 1, 1.5, 2, 2.5)

.inference_pilot_methods <- function() c(
  "one_round_crossfit", "target_only", "sample_size",
  "inverse_variance", "federated_dr", "pooled_dr"
)

#' Construct the frozen independent-inference pilot manifest.
#'
#' This helper is intentionally standalone: sourcing this file does not run the
#' pilot.  The returned rows are one rho-group task each, not six scientific
#' result rows.
roce_inference_pilot_manifest <- function() {
  task_id <- seq_len(100L)
  data.frame(
    task_id = task_id,
    sim_id = 10000L + task_id,
    experiment = "independent_inference_pilot_v19",
    config = "C1",
    p = 100L,
    K = 2L,
    cutoff = 1,
    n_site = 1000L,
    n_folds = 5L,
    nlambda_init = 100L,
    nuisance_lambda_rule = "min",
    outcome_type = "binary",
    estimand_type = "superpopulation",
    n_bootstrap = 5000L,
    n_weight_bootstrap = 1000L,
    M_tau = 5,
    M_tau_inference = 5,
    include_hard_threshold_diagnostic = TRUE,
    methods = paste(.inference_pilot_methods(), collapse = ","),
    rho_values = paste(.inference_pilot_rhos(), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

.inference_pilot_scalar_equal <- function(x, y) {
  length(x) == length(y) && !anyNA(x) && !anyNA(y) &&
    identical(as.character(x), as.character(y))
}

roce_validate_inference_pilot_manifest <- function(manifest, task_id) {
  expected <- roce_inference_pilot_manifest()
  if (!is.data.frame(manifest) || !identical(names(manifest), names(expected)) ||
      nrow(manifest) != 100L) {
    stop("inference pilot manifest must have the exact frozen 100-row schema.",
         call. = FALSE)
  }
  for (name in names(expected)) {
    value <- manifest[[name]]
    if (!is.atomic(value) || !is.null(dim(value))) {
      stop("inference pilot manifest column '", name,
           "' must be one-dimensional and atomic.", call. = FALSE)
    }
    if (!.inference_pilot_scalar_equal(manifest[[name]], expected[[name]])) {
      stop("inference pilot manifest mismatch for '", name, "'.",
           call. = FALSE)
    }
  }
  if (!is.atomic(task_id) || !is.null(dim(task_id)) || is.object(task_id) ||
      is.logical(task_id) || is.complex(task_id) ||
      !typeof(task_id) %in% c("integer", "double", "character")) {
    stop("TASK_ID must be an unclassed integer, double, or character scalar.",
         call. = FALSE)
  }
  numeric_task <- suppressWarnings(as.numeric(task_id))
  if (length(numeric_task) != 1L || is.na(numeric_task) ||
      !is.finite(numeric_task) || numeric_task != floor(numeric_task) ||
      numeric_task < 1L || numeric_task > 100L) {
    stop("TASK_ID must be one integer in [1, 100].", call. = FALSE)
  }
  expected[as.integer(numeric_task), , drop = FALSE]
}

roce_build_inference_audit <- function(results) {
  required <- c(
    "sim_id", "rho", "method", "estimate", "truth", "se", "ci_lower",
    "ci_upper", "coverage", "se_fixed_weight_bootstrap",
    "se_weight_relearn_bootstrap", "variance_weight_relearn_bootstrap",
    "weight_uncertainty_ratio", "weight_relearn_n_bootstrap",
    "weight_bootstrap_failures"
  )
  missing <- setdiff(required, names(results))
  if (!is.data.frame(results) || length(missing)) {
    stop("inference audit input is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  soft <- results[results$method == "one_round_crossfit_ate", , drop = FALSE]
  expected_rhos <- .inference_pilot_rhos()
  if (!is.numeric(soft$rho) || !is.null(dim(soft$rho)) ||
      any(!is.finite(soft$rho)) || nrow(soft) != 6L ||
      anyDuplicated(soft$rho) ||
      !identical(sort(as.numeric(soft$rho)), expected_rhos)) {
    stop("inference audit requires exactly one soft-TATE row per frozen rho.",
         call. = FALSE)
  }
  soft <- soft[match(expected_rhos, soft$rho), , drop = FALSE]
  if (!is.numeric(soft$sim_id) || !is.null(dim(soft$sim_id)) ||
      any(!is.finite(soft$sim_id)) || any(soft$sim_id != floor(soft$sim_id)) ||
      length(unique(soft$sim_id)) != 1L ||
      soft$sim_id[[1L]] < 10001L || soft$sim_id[[1L]] > 10100L) {
    stop("inference audit rows must share one valid pilot sim_id.",
         call. = FALSE)
  }
  numeric_names <- c(
    "estimate", "truth", "se", "ci_lower", "ci_upper",
    "se_fixed_weight_bootstrap", "se_weight_relearn_bootstrap",
    "variance_weight_relearn_bootstrap", "weight_uncertainty_ratio",
    "weight_relearn_n_bootstrap", "weight_bootstrap_failures"
  )
  if (any(vapply(soft[numeric_names], function(value) {
    !is.numeric(value) || any(!is.finite(value))
  }, logical(1L))) || any(soft$se <= 0) ||
      any(soft$se_fixed_weight_bootstrap <= 0) ||
      any(soft$se_weight_relearn_bootstrap <= 0) ||
      any(soft$weight_relearn_n_bootstrap != 1000L) ||
      any(soft$weight_bootstrap_failures != 0L) ||
      any(abs(soft$variance_weight_relearn_bootstrap -
              soft$se_weight_relearn_bootstrap^2) > 1e-10)) {
    stop("inference audit found invalid bootstrap or analytic scalars.",
         call. = FALSE)
  }
  ratio_expected <- soft$se_weight_relearn_bootstrap /
    soft$se_fixed_weight_bootstrap
  ratio_relative_error <- abs(soft$weight_uncertainty_ratio - ratio_expected) /
    pmax(abs(ratio_expected), .Machine$double.eps)
  if (any(ratio_relative_error > 1e-10)) {
    stop("weight uncertainty ratio identity failed.", call. = FALSE)
  }
  coverage_value <- soft$coverage
  if (is.logical(coverage_value)) {
    if (anyNA(coverage_value)) {
      stop("coverage must contain only logical or numeric 0/1 values.",
           call. = FALSE)
    }
    coverage_logical <- coverage_value
  } else if (is.numeric(coverage_value) && is.null(dim(coverage_value)) &&
             !anyNA(coverage_value) && all(is.finite(coverage_value)) &&
             all(coverage_value %in% c(0, 1))) {
    coverage_logical <- coverage_value == 1
  } else {
    stop("coverage must contain only logical or numeric 0/1 values.",
         call. = FALSE)
  }
  analytic_lower <- soft$estimate - stats::qnorm(0.975) * soft$se
  analytic_upper <- soft$estimate + stats::qnorm(0.975) * soft$se
  analytic_coverage <- soft$truth >= soft$ci_lower & soft$truth <= soft$ci_upper
  if (max(abs(soft$ci_lower - analytic_lower),
          abs(soft$ci_upper - analytic_upper)) > 1e-10 ||
      !identical(coverage_logical, analytic_coverage)) {
    stop("stored analytic CI or coverage identity failed.", call. = FALSE)
  }
  relearn_lower <- soft$estimate -
    stats::qnorm(0.975) * soft$se_weight_relearn_bootstrap
  relearn_upper <- soft$estimate +
    stats::qnorm(0.975) * soft$se_weight_relearn_bootstrap
  data.frame(
    sim_id = soft$sim_id,
    rho = soft$rho,
    method = soft$method,
    estimate = soft$estimate,
    truth = soft$truth,
    analytic_se = soft$se,
    analytic_ci_lower = soft$ci_lower,
    analytic_ci_upper = soft$ci_upper,
    analytic_coverage = analytic_coverage,
    fixed_weight_bootstrap_se = soft$se_fixed_weight_bootstrap,
    weight_relearn_bootstrap_se = soft$se_weight_relearn_bootstrap,
    weight_relearn_ci_lower = relearn_lower,
    weight_relearn_ci_upper = relearn_upper,
    weight_relearn_coverage =
      soft$truth >= relearn_lower & soft$truth <= relearn_upper,
    weight_uncertainty_ratio = soft$weight_uncertainty_ratio,
    n_weight_bootstrap = soft$weight_relearn_n_bootstrap,
    inference_status = "diagnostic_only",
    stringsAsFactors = FALSE
  )
}

.inference_pilot_output_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || basename(path) %in% c("", ".", "..")) {
    stop("OUTPUT_DIR must be one concrete directory path.", call. = FALSE)
  }
  file.path(normalizePath(dirname(path), mustWork = FALSE), basename(path))
}

.inference_pilot_validate_scheduler_task <- function(task_id) {
  scheduler_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
  if (!nzchar(scheduler_id)) return(invisible(TRUE))
  if (!grepl("^[1-9][0-9]*$", scheduler_id) ||
      suppressWarnings(as.numeric(scheduler_id)) != as.numeric(task_id)) {
    stop("SLURM_ARRAY_TASK_ID must be a positive integer equal to TASK_ID.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.inference_pilot_write_failure <- function(
    output, task, message, stage, hashes, sha256_file, atomic_writer) {
  writer <- function(directory) {
    attempt <- data.frame(
      task_id = task$task_id, sim_id = task$sim_id, status = "failed",
      failure_stage = stage, failure_message = as.character(message),
      retry_authorized = FALSE, stringsAsFactors = FALSE
    )
    utils::write.csv(attempt, file.path(directory, "attempt_failure.csv"),
                     row.names = FALSE, na = "")
    writeLines(c(
      "independent_inference_pilot=failed",
      paste0("task_id=", task$task_id), paste0("sim_id=", task$sim_id),
      paste0("failure_stage=", stage), "retry_authorized=FALSE",
      paste0("package_fingerprint=", hashes$package),
      paste0("workflow_fingerprint=", hashes$workflow),
      paste0("manifest_fingerprint=", hashes$manifest)
    ), file.path(directory, "metadata.txt"), useBytes = TRUE)
    files <- c("attempt_failure.csv", "metadata.txt")
    digests <- vapply(file.path(directory, files), sha256_file, character(1L))
    writeLines(paste(digests, files, sep = "  "),
               file.path(directory, "sha256.txt"), useBytes = TRUE)
  }
  atomic_writer(output, writer, caller = "independent inference pilot failure")
}

roce_publish_inference_pilot <- function(
    output, task, hashes, writer, sha256_file, atomic_writer) {
  publication_error <- tryCatch({
    atomic_writer(output, writer, caller = "independent inference pilot")
    NULL
  }, error = function(error) error)
  if (is.null(publication_error)) return(invisible(output))
  if (!file.exists(output) && !dir.exists(output)) {
    # Failure publication is best effort.  The original publication condition
    # remains authoritative and is always rethrown below.
    try(.inference_pilot_write_failure(
      output, task, conditionMessage(publication_error), "publication",
      hashes, sha256_file, atomic_writer
    ), silent = TRUE)
  }
  stop(publication_error)
}

run_independent_inference_pilot_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop("usage: run_independent_inference_pilot.R MANIFEST TASK_ID OUTPUT_DIR",
         call. = FALSE)
  }
  project_root <- normalizePath(Sys.getenv("ROCE_PROJECT_ROOT", getwd()),
                                mustWork = TRUE)
  manifest_path <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- .inference_pilot_output_path(args[[3L]])
  if (file.exists(output) || dir.exists(output)) {
    stop("independent inference pilot output already exists: ", output,
         call. = FALSE)
  }

  source(file.path(project_root, "scripts/slurm/atomic_output.R"))
  source(file.path(project_root, "scripts/slurm/result_provenance.R"))
  source(file.path(project_root, "scripts/slurm/resource_topology.R"))
  source(file.path(project_root, "scripts/slurm/direct_tate_task_helpers.R"))
  source(file.path(project_root, "scripts/slurm/simulation_qc_policy.R"))
  source(file.path(project_root,
                   "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R"))

  task <- roce_validate_inference_pilot_manifest(
    read.csv(manifest_path, stringsAsFactors = FALSE), args[[2L]]
  )
  .inference_pilot_validate_scheduler_task(task$task_id)
  expected_basename <- sprintf("seed_%06d", as.integer(task$sim_id))
  if (!identical(basename(output), expected_basename)) {
    stop("OUTPUT_DIR basename must be ", expected_basename, ".",
         call. = FALSE)
  }
  project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
  if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  package_provenance <- roce_runtime_package_provenance(project_library)
  package_hash <- .weight_calibration_installed_package_fingerprint(
    package_provenance$library, roce_sha256_file
  )
  frozen_package_hash <-
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7"
  workflow_files <- file.path(project_root, c(
    "diagnosis/tate_common_weight/run_independent_inference_pilot.R",
    "diagnosis/tate_common_weight/run_independent_inference_pilot.sh",
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R",
    "scripts/slurm/atomic_output.R", "scripts/slurm/result_provenance.R",
    "scripts/slurm/resource_topology.R",
    "scripts/slurm/direct_tate_task_helpers.R",
    "scripts/slurm/simulation_qc_policy.R",
    "scripts/slurm/package_library_utils.sh"
  ))
  hashes <- list(
    package = package_provenance$fingerprint,
    workflow = .weight_calibration_files_fingerprint(workflow_files,
                                                      roce_sha256_file),
    manifest = roce_sha256_file(manifest_path)
  )
  if (!identical(package_hash, package_provenance$fingerprint) ||
      !identical(package_hash, frozen_package_hash) ||
      !identical(hashes$workflow,
                 roce_sha256_environment("ROCE_WORKFLOW_FINGERPRINT")) ||
      !identical(hashes$manifest,
                 roce_sha256_environment("ROCE_MANIFEST_FINGERPRINT"))) {
    stop("independent inference pilot provenance gate failed.", call. = FALSE)
  }
  allocated <- roce_read_positive_integer_env("SLURM_CPUS_PER_TASK", 1L)
  cv_threads <- roce_read_positive_integer_env(
    "ROCE_NUISANCE_CV_THREADS", 5L, 5L
  )
  if (cv_threads != 5L) stop("pilot requires exactly five CV threads.", call. = FALSE)
  resource_plan <- roce_slurm_resource_plan(
    source_count = 2L, n_folds = 5L, allocated_cores = allocated,
    nuisance_cv_threads = cv_threads
  )
  cores_per_rho <- cv_threads *
    if (resource_plan$parallel_treatment_arms) 2L else 1L
  rho_workers <- min(5L, max(1L, allocated %/% cores_per_rho))
  release_lock <- roce_claim_task_lock(paste0(output, ".lock"), c(
    paste0("task_id=", task$task_id), paste0("sim_id=", task$sim_id),
    paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local"))
  ))
  on.exit(release_lock(), add = TRUE)

  simulation_args <- list(
    sim_id = as.integer(task$sim_id), n_total = 3000L, K = 2L, p = 100L,
    config = "C1", methods = .inference_pilot_methods(), verbose = FALSE,
    n_cores_internal = resource_plan$source_workers, nlambda_init = 100L,
    nuisance_lambda_rule = "min", estimand_type = "superpopulation",
    outcome_type = "binary", n_folds = 5L, aggregation_lambda = 1,
    n_bootstrap = 5000L, n_weight_bootstrap = 1000L,
    M_tau = 5, M_tau_inference = 5, estimate_ate = TRUE,
    parallel_treatment_arms = resource_plan$parallel_treatment_arms,
    include_hard_threshold_diagnostic = TRUE, dgp_type = "face"
  )
  started <- proc.time()[["elapsed"]]
  scientific <- tryCatch({
    grouped <- RoCE:::.run_face_rho_group(
      simulation_args = simulation_args, rho_values = .inference_pilot_rhos(),
      changed_sources = "s1", artifact_rhos = .inference_pilot_rhos(),
      positive_rho_workers = rho_workers
    )
    expected_keys <- vapply(.inference_pilot_rhos(), format, character(1L),
                            scientific = FALSE, trim = TRUE)
    if (!identical(names(grouped$artifacts), expected_keys) ||
        any(vapply(grouped$artifacts, function(entry) {
          is.null(entry$data_split) ||
            is.null(entry$direct_tate_results$one_round_crossfit)
        }, logical(1L)))) stop("pilot did not retain six complete artifacts.")
    rows <- RoCE:::.bind_sim_result_list(grouped$results)
    rows$task_id <- as.integer(task$task_id)
    rows$rho_group_size <- 6L
    rows$rho_group_manifest_fingerprint <- hashes$manifest
    rows$allocated_cores <- resource_plan$allocated_cores
    rows$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
    rows$source_workers_per_arm <- resource_plan$source_workers
    rows$parallel_treatment_arms <- resource_plan$parallel_treatment_arms
    rows$rho_group_positive_rho_workers <- rho_workers
    rows$rho_group_positive_rho_backend <- grouped$reuse$positive_rho_backend
    rows$package_library <- package_provenance$library
    rows$package_fingerprint <- package_provenance$fingerprint
    rows$workflow_fingerprint <- hashes$workflow
    rows$manifest_fingerprint <- hashes$manifest
    scheduler <- roce_scheduler_provenance()
    for (name in names(scheduler)) rows[[name]] <- scheduler[[name]]
    validation_task <- list(
      rho_values = .inference_pilot_rhos(), sim_id = as.integer(task$sim_id),
      methods = .inference_pilot_methods()
    )
    .weight_calibration_validate_results(
      rows, validation_task, 1000L, hashes$package, hashes$workflow,
      hashes$manifest, hashes$manifest
    )
    for (rho in .inference_pilot_rhos()) {
      key <- format(rho, scientific = FALSE, trim = TRUE)
      fit <- grouped$artifacts[[key]]$direct_tate_results$one_round_crossfit
      RoCE:::.validate_weight_bootstrap_result(fit)
      row <- rows[rows$rho == rho &
                    rows$method == "one_round_crossfit_ate", , drop = FALSE]
      if (nrow(row) != 1L || abs(row$estimate - fit$estimate) > 1e-12 ||
          abs(row$se - fit$se) > 1e-12) {
        stop("artifact/result identity failed for rho=", rho, ".")
      }
    }
    inference <- roce_build_inference_audit(rows)
    qc <- RoCE:::diagnose_simulation_results(rows, expected_replications = 1L)
    classification <- roce_classify_simulation_qc(qc)
    qc$implementation_failed <- classification$implementation_failed
    qc$statistical_review_required <- classification$statistical_review_required
    if (any(qc$implementation_failed) ||
        any(grepl("weight_relearn_bootstrap_failures_detected",
                  qc$diagnostic_status, fixed = TRUE))) {
      stop("pilot diagnostic QC failed closed.")
    }
    list(grouped = grouped, rows = rows, inference = inference, qc = qc)
  }, error = function(error) error)
  elapsed <- proc.time()[["elapsed"]] - started
  if (inherits(scientific, "error")) {
    .inference_pilot_write_failure(
      output, task, conditionMessage(scientific), "scientific_or_validation",
      hashes, roce_sha256_file, roce_write_atomic_directory
    )
    stop("independent inference pilot failed; attempt record published: ",
         conditionMessage(scientific), call. = FALSE)
  }

  writer <- function(directory) {
    if (!identical(.weight_calibration_installed_package_fingerprint(
      package_provenance$library, roce_sha256_file), hashes$package) ||
        !identical(.weight_calibration_files_fingerprint(
          workflow_files, roce_sha256_file), hashes$workflow) ||
        !identical(roce_sha256_file(manifest_path), hashes$manifest)) {
      stop("pilot inputs changed before atomic publication.", call. = FALSE)
    }
    utils::write.csv(scientific$rows, file.path(directory, "results.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(scientific$inference,
                     file.path(directory, "inference_audit.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(scientific$qc,
                     file.path(directory, "diagnostic_qc.csv"),
                     row.names = FALSE, na = "")
    saveRDS(list(
      group_result = scientific$grouped, task = task,
      simulation_args = simulation_args, resource_plan = resource_plan,
      positive_rho_workers = rho_workers
    ), file.path(directory, "artifacts.rds"))
    writeLines(c(
      "independent_inference_pilot=complete",
      paste0("task_id=", task$task_id), paste0("sim_id=", task$sim_id),
      paste0("rho_values=", task$rho_values),
      "n_weight_bootstrap=1000", "comparison_bootstrap=5000",
      "primary_estimate_se_ci_unchanged=TRUE",
      "weight_relearn_interval_status=diagnostic_only",
      paste0("package_library=", package_provenance$library),
      paste0("package_fingerprint=", hashes$package),
      paste0("workflow_fingerprint=", hashes$workflow),
      paste0("manifest_fingerprint=", hashes$manifest),
      paste0("elapsed_seconds=", format(elapsed, digits = 17)),
      paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
    ), file.path(directory, "metadata.txt"), useBytes = TRUE)
    files <- c("results.csv", "inference_audit.csv", "diagnostic_qc.csv",
               "artifacts.rds", "metadata.txt")
    digests <- vapply(file.path(directory, files), roce_sha256_file,
                      character(1L))
    writeLines(paste(digests, files, sep = "  "),
               file.path(directory, "sha256.txt"), useBytes = TRUE)
  }
  roce_publish_inference_pilot(
    output, task, hashes, writer, roce_sha256_file,
    roce_write_atomic_directory
  )
  message("[done] atomically wrote independent inference pilot: ", output)
  invisible(output)
}

if (sys.nframe() == 0L) run_independent_inference_pilot_main()
