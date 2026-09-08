#!/usr/bin/env Rscript
# One diagnostic draw per job; draw zero is an exact-data identity refit.
# Run from repository root. No primary results or package files are modified.
# Current protocol uses origin-grouped nuisance CV; the exact former row-CV
# driver is preserved under full_refit_row_cv_v18_protocol_source/.

roce_full_refit_draw_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) {
    stop("Usage: run_full_refit_draw.R BUNDLE_DIR RHO DRAW_ID OUTPUT_DIR")
  }
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("scripts/slurm/direct_tate_task_helpers.R")
  source("diagnosis/tate_common_weight/full_refit_resampling.R")
  source("diagnosis/tate_common_weight/nuisance_cv_partition_audit.R")
  bundle <- normalizePath(args[[1L]], mustWork = TRUE)
  rho <- suppressWarnings(as.numeric(args[[2L]]))
  if (length(rho) != 1L || is.na(rho) || !rho %in% c(0, .5, 1, 1.5, 2, 2.5)) {
    stop("RHO must be one of the six prespecified values.")
  }
  draw_id <- .weight_calibration_integer(args[[3L]], "DRAW_ID", 0L, 1000L)
  output <- args[[4L]]
  if (file.exists(output) || dir.exists(output)) stop("Output already exists.")
  payload_names <- c("results.csv", "diagnostic_qc.csv", "artifacts.rds", "metadata.txt")
  hash_lines <- readLines(file.path(bundle, "sha256.txt"), warn = FALSE)
  expected_hashes <- substr(hash_lines, 1L, 64L)
  payloads <- substring(hash_lines, 67L)
  if (length(hash_lines) != 4L || !setequal(payloads, payload_names) ||
      anyDuplicated(payloads) || any(substr(hash_lines, 65L, 66L) != "  ") ||
      any(!grepl("^[0-9a-f]{64}$", expected_hashes))) {
    stop("Malformed source-bundle checksum manifest.")
  }
  verify_bundle <- function() {
    observed <- vapply(file.path(bundle, payloads), roce_sha256_file, character(1L))
    if (!identical(unname(observed), expected_hashes)) stop("Source bundle checksum mismatch.")
    invisible(TRUE)
  }
  verify_bundle()
  metadata_lines <- readLines(file.path(bundle, "metadata.txt"), warn = FALSE)
  metadata <- stats::setNames(sub("^[^=]*=", "", metadata_lines),
                              sub("=.*$", "", metadata_lines))
  if (anyDuplicated(names(metadata))) stop("Duplicate source metadata keys.")
  locked_metadata <- c(
    package_fingerprint = "27ebfb159e6e7f6db955be810e3ef89724741a6f02f68e83730a8f808e774b45",
    workflow_fingerprint = "10db72dc532e3fb8811c4f6a297f9ad9b8e3fd09cf19604ee7d6228d15392831",
    manifest_fingerprint = "35b8a5d09a59c96d21f56edb2a7fb67ba0dc29e698560766ffca6e479f5d8c43",
    group_manifest_fingerprint = "ad4af6e38e2fb1f6ce69a02a9bfdc51ae8897fe6e15174a377e80448a53a8341"
  )
  if (!identical(metadata[names(locked_metadata)], locked_metadata)) {
    stop("This bounded diagnostic requires the frozen v18 C1/K2 calibration.")
  }
  project_library <- normalizePath(Sys.getenv("ROCE_PROJECT_LIB"), mustWork = TRUE)
  .libPaths(c(project_library, .libPaths()))
  library(RoCE)
  if (!identical(normalizePath(find.package("RoCE")),
                 normalizePath(file.path(project_library, "RoCE")))) {
    stop("Wrong RoCE installation loaded.")
  }
  package_hash <- .weight_calibration_installed_package_fingerprint(
    project_library, roce_sha256_file
  )
  if (!identical(package_hash,
                 "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")) {
    stop("Origin-grouped refitting requires the frozen tested v19 installation.")
  }
  gate <- readLines(file.path(project_library, "audit_tests_passed.txt"))
  if (!"package_tests=passed" %in% gate ||
      !paste0("package_fingerprint=", package_hash) %in% gate) {
    stop("Matching installed-package test gate is absent.")
  }
  workflow_files <- c(
    "diagnosis/tate_common_weight/full_refit_resampling.R",
    "diagnosis/tate_common_weight/nuisance_cv_partition_audit.R",
    "diagnosis/tate_common_weight/run_full_refit_draw.R",
    "diagnosis/tate_common_weight/run_full_refit_draw.sh",
    "diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R",
    "scripts/slurm/result_provenance.R", "scripts/slurm/atomic_output.R",
    "scripts/slurm/direct_tate_task_helpers.R"
  )
  workflow_hash <- .weight_calibration_files_fingerprint(workflow_files, roce_sha256_file)
  source_bundle_hash <- roce_sha256_file(file.path(bundle, "sha256.txt"))
  saved <- readRDS(file.path(bundle, "artifacts.rds"))
  locked_args <- list(dgp_type = "face", config = "C1", K = 2L,
                      n_total = 3000L, p = 100L, n_folds = 5L,
                      nlambda_init = 100L, nuisance_lambda_rule = "min",
                      aggregation_lambda = 1, M_tau = 5, M_tau_inference = 5)
  for (name in names(locked_args)) {
    .weight_calibration_scalar_equal(saved$simulation_args[[name]],
                                    locked_args[[name]], name)
  }
  seed_id <- .weight_calibration_integer(saved$task$sim_id, "sim_id", 1L, 10L)
  rho_key <- format(rho, scientific = FALSE, trim = TRUE)
  entry <- saved$group_result$artifacts[[rho_key]]
  if (is.null(entry)) stop("Source bundle is missing the selected rho artifact.")
  rm(saved)
  gc(verbose = FALSE)
  reference <- entry$direct_tate_results$one_round_crossfit
  source_rows <- read.csv(file.path(bundle, "results.csv"), stringsAsFactors = FALSE)
  source_row <- source_rows[source_rows$rho == rho &
                            source_rows$method == "one_round_crossfit_ate", ]
  if (nrow(source_row) != 1L || source_row$sim_id != seed_id ||
      abs(source_row$estimate - reference$estimate) > 1e-12 ||
      abs(source_row$se - reference$se) > 1e-12) {
    stop("Source artifact and result row disagree.")
  }
  # Same seed/draw mapping for every rho, preserving paired resampling IDs.
  resample_seed <- as.integer(190000L + seed_id * 10000L + draw_id)
  cores <- .weight_calibration_integer(Sys.getenv("ROCE_REFIT_CORES", "1"),
                                       "ROCE_REFIT_CORES", 1L, 2L)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  release_lock <- roce_claim_task_lock(paste0(output, ".lock"), c(
    paste0("sim_id=", seed_id), paste0("rho=", rho), paste0("draw_id=", draw_id),
    paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local"))
  ))
  on.exit(release_lock(), add = TRUE)
  warnings <- character()
  started <- proc.time()[["elapsed"]]
  result <- tryCatch(withCallingHandlers(
    roce_refit_once(entry$data_split, reference, resample_seed,
                   identity = draw_id == 0L, n_cores = cores,
                   parallel_arms = FALSE, nuisance_cv_grouping = "origin"),
    warning = function(w) warnings <<- c(warnings, conditionMessage(w))
  ), error = function(e) list(failure_message = conditionMessage(e),
                              failure_stage = "setup_or_postfit_validation"))
  succeeded <- is.null(result$failure_message)
  summary <- data.frame(
    sim_id = seed_id, rho = rho, draw_id = draw_id,
    resample_seed = resample_seed, identity = draw_id == 0L,
    status = if (succeeded) "completed" else "failed",
    failure_message = if (succeeded) "" else result$failure_message,
    failure_stage = if (succeeded) "" else result$failure_stage,
    estimate_refit_relearned_weights = if (succeeded)
      result$estimate_refit_relearned_weights else NA_real_,
    estimate_refit_original_weights = if (succeeded)
      result$estimate_refit_original_weights else NA_real_,
    estimate_fixed_nuisance_original_weights = if (!is.null(result$matched_fixed_nuisance))
      result$matched_fixed_nuisance$estimate_original_weights else NA_real_,
    estimate_fixed_nuisance_relearned_weights = if (!is.null(result$matched_fixed_nuisance))
      result$matched_fixed_nuisance$estimate_relearned_weights else NA_real_,
    estimate_target_refit = if (succeeded) result$fitted$target_only$estimate else NA_real_,
    estimate_reference = reference$estimate,
    se_analytic_reference = reference$se,
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    nuisance_refit = TRUE, inference_validated = FALSE,
    conditional_on_saved_design_and_outer_partition = TRUE,
    nuisance_cv_partition = "origin_grouped",
    nuisance_cv_groups_duplicate_origins = TRUE,
    nuisance_cv_duplicate_origin_leakage_possible = FALSE,
    nuisance_cv_partition_count = length(result$nuisance_cv_partitions),
    nuisance_cv_partition_audit_passed = succeeded &&
      length(result$nuisance_cv_partitions) > 0L &&
      all(vapply(result$nuisance_cv_partitions, roce_cv_partition_record_valid, logical(1L))),
    resampling_scheme = "site_by_original_outer_fold_multinomial",
    package_fingerprint = package_hash, workflow_fingerprint = workflow_hash,
    source_package_fingerprint = unname(metadata["package_fingerprint"]),
    source_bundle_fingerprint = source_bundle_hash,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", "local"),
    slurm_node = Sys.getenv("SLURMD_NODENAME", "unknown"),
    allocated_cpus = Sys.getenv("SLURM_CPUS_PER_TASK", "unknown"),
    stringsAsFactors = FALSE
  )
  if (succeeded && any(!is.finite(unlist(summary[c(
    "estimate_refit_relearned_weights", "estimate_refit_original_weights",
    "estimate_fixed_nuisance_original_weights", "estimate_fixed_nuisance_relearned_weights",
    "estimate_target_refit"
  )])))) stop("Nonfinite full-refit estimates.")
  verify_bundle()
  if (!identical(package_hash, .weight_calibration_installed_package_fingerprint(
      project_library, roce_sha256_file)) ||
      !identical(workflow_hash, .weight_calibration_files_fingerprint(
        workflow_files, roce_sha256_file))) stop("Inputs changed during refit.")
  roce_write_atomic_directory(output, function(staging) {
    write.csv(summary, file.path(staging, "result.csv"), row.names = FALSE)
    saveRDS(list(result = result, summary = summary, warnings = warnings,
                 source_bundle = bundle), file.path(staging, "draw.rds"))
    files <- c("result.csv", "draw.rds")
    hashes <- vapply(file.path(staging, files), roce_sha256_file, character(1L))
    writeLines(paste(hashes, files, sep = "  "), file.path(staging, "sha256.txt"))
  }, caller = "full nuisance-refit diagnostic")
  print(summary)
  if (!succeeded) stop("Full refit failed; failure record retained, no silent redraw.")
  invisible(summary)
}

if (sys.nframe() == 0L) roce_full_refit_draw_main()
