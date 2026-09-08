#!/usr/bin/env Rscript

.siip_check_sha256 <- function(directory, sha256_file, expected_files = NULL) {
  manifest <- file.path(directory, "sha256.txt")
  if (!file.exists(manifest)) stop("missing checksum manifest: ", directory)
  lines <- readLines(manifest, warn = FALSE)
  if (!length(lines) || any(substr(lines, 65L, 66L) != "  ") ||
      any(!grepl("^[0-9a-f]{64}  [^/]+$", lines)) ||
      anyDuplicated(substring(lines, 67L))) stop("malformed checksum manifest: ", directory)
  files <- substring(lines, 67L)
  if (!is.null(expected_files) && !setequal(files, expected_files)) {
    stop("checksum payload schema mismatch: ", directory)
  }
  if (any(!file.exists(file.path(directory, files))) ||
      any(vapply(file.path(directory, files), sha256_file, character(1L)) !=
          substr(lines, 1L, 64L))) stop("checksum mismatch: ", directory)
  sha256_file(manifest)
}

.siip_read_gate <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (!length(lines) || any(!grepl("^[^=]+=[^=]*$", lines))) {
    stop("malformed seed audit gate: ", path)
  }
  keys <- sub("=.*$", "", lines)
  if (anyDuplicated(keys)) stop("duplicate seed audit gate keys: ", path)
  stats::setNames(sub("^[^=]*=", "", lines), keys)
}

.siip_compare_frames <- function(observed, expected, tolerance = 1e-12) {
  if (!is.data.frame(observed) || !is.data.frame(expected) ||
      !identical(names(observed), names(expected)) ||
      !identical(dim(observed), dim(expected))) return(FALSE)
  all(vapply(names(expected), function(name) {
    x <- observed[[name]]; y <- expected[[name]]
    if (is.numeric(y)) {
      is.numeric(x) && identical(is.na(x), is.na(y)) &&
        all(abs(x[!is.na(y)] - y[!is.na(y)]) <= tolerance)
    } else {
      identical(x, y)
    }
  }, logical(1L)))
}

.siip_binom_interval <- function(covered) {
  covered <- as.logical(covered)
  unname(stats::binom.test(sum(covered), length(covered), conf.level = 0.95)$conf.int)
}

.siip_mcse <- function(values) {
  if (length(values) < 2L) return(NA_real_)
  stats::sd(values) / sqrt(length(values))
}

.siip_sd <- function(values) if (length(values) < 2L) NA_real_ else stats::sd(values)

.siip_ratio <- function(numerator, denominator) {
  if (length(denominator) != 1L || is.na(denominator) || denominator <= 0) NA_real_
  else numerator / denominator
}

roce_summarize_independent_inference <- function(raw_rows, inference_rows) {
  rhos <- .inference_pilot_rhos()
  expected_ids <- sort(unique(raw_rows$sim_id))
  summaries <- lapply(rhos, function(rho) {
    soft <- raw_rows[raw_rows$rho == rho &
                       raw_rows$method == "one_round_crossfit_ate", , drop = FALSE]
    target <- raw_rows[raw_rows$rho == rho &
                         raw_rows$method == "target_only_ate", , drop = FALSE]
    audit <- inference_rows[inference_rows$rho == rho, , drop = FALSE]
    if (!nrow(soft) || nrow(soft) != nrow(target) || nrow(soft) != nrow(audit) ||
        anyDuplicated(soft$sim_id) || anyDuplicated(target$sim_id) ||
        anyDuplicated(audit$sim_id) ||
        !setequal(soft$sim_id, expected_ids) ||
        !is.logical(target$coverage) || anyNA(target$coverage) ||
        !is.logical(audit$analytic_coverage) || anyNA(audit$analytic_coverage) ||
        !is.logical(audit$weight_relearn_coverage) ||
        anyNA(audit$weight_relearn_coverage)) {
      stop("incomplete or invalid per-rho paired estimator rows.")
    }
    ids <- sort(soft$sim_id)
    soft <- soft[match(ids, soft$sim_id), , drop = FALSE]
    target <- target[match(ids, target$sim_id), , drop = FALSE]
    audit <- audit[match(ids, audit$sim_id), , drop = FALSE]
    if (anyNA(target$sim_id) || anyNA(audit$sim_id)) stop("unpaired sim IDs within rho.")
    if (any(target$truth != soft$truth)) stop("paired methods disagree on truth within rho.")
    error <- soft$estimate - soft$truth
    target_error <- target$estimate - target$truth
    squared_difference <- error^2 - target_error^2
    empirical_sd <- .siip_sd(soft$estimate)
    analytic_coverage <- as.logical(audit$analytic_coverage)
    relearn_coverage <- as.logical(audit$weight_relearn_coverage)
    target_coverage <- as.logical(target$coverage)
    analytic_interval <- .siip_binom_interval(analytic_coverage)
    relearn_interval <- .siip_binom_interval(relearn_coverage)
    analytic_difference <- as.numeric(analytic_coverage) - as.numeric(target_coverage)
    relearn_difference <- as.numeric(relearn_coverage) - as.numeric(target_coverage)
    relearn_minus_analytic <-
      as.numeric(relearn_coverage) - as.numeric(analytic_coverage)
    data.frame(
      rho = rho, n_independent_seeds = length(ids),
      mean_error = mean(error), mean_error_mcse = .siip_mcse(error),
      empirical_sd = empirical_sd, rmse = sqrt(mean(error^2)),
      target_rmse = sqrt(mean(target_error^2)),
      paired_mse_difference_soft_minus_target = mean(squared_difference),
      paired_mse_difference_mcse = .siip_mcse(squared_difference),
      mean_analytic_se = mean(audit$analytic_se),
      mean_fixed_weight_bootstrap_se = mean(audit$fixed_weight_bootstrap_se),
      mean_weight_relearn_bootstrap_se = mean(audit$weight_relearn_bootstrap_se),
      analytic_se_to_empirical_sd = .siip_ratio(mean(audit$analytic_se), empirical_sd),
      fixed_weight_se_to_empirical_sd = .siip_ratio(
        mean(audit$fixed_weight_bootstrap_se), empirical_sd),
      weight_relearn_se_to_empirical_sd = .siip_ratio(
        mean(audit$weight_relearn_bootstrap_se), empirical_sd),
      mean_analytic_ci_length = mean(audit$analytic_ci_upper - audit$analytic_ci_lower),
      mean_weight_relearn_ci_length = mean(
        audit$weight_relearn_ci_upper - audit$weight_relearn_ci_lower),
      analytic_coverage = mean(analytic_coverage),
      analytic_coverage_exact_lower = analytic_interval[[1L]],
      analytic_coverage_exact_upper = analytic_interval[[2L]],
      weight_relearn_coverage = mean(relearn_coverage),
      weight_relearn_coverage_exact_lower = relearn_interval[[1L]],
      weight_relearn_coverage_exact_upper = relearn_interval[[2L]],
      analytic_coverage_minus_target = mean(analytic_difference),
      analytic_coverage_minus_target_mcse = .siip_mcse(analytic_difference),
      weight_relearn_coverage_minus_target = mean(relearn_difference),
      weight_relearn_coverage_minus_target_mcse = .siip_mcse(relearn_difference),
      weight_relearn_coverage_minus_analytic = mean(relearn_minus_analytic),
      weight_relearn_coverage_minus_analytic_mcse =
        .siip_mcse(relearn_minus_analytic),
      warning_count = NA_integer_,
      warning_capture_complete = FALSE,
      inference_validated = FALSE,
      statistical_review_required = TRUE,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summaries)
}

summarize_independent_inference_pilot_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: summarize_independent_inference_pilot.R ROOT N OUTPUT_DIR")
  root <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- file.path(normalizePath(dirname(args[[3L]]), mustWork = FALSE), basename(args[[3L]]))
  if (file.exists(output) || dir.exists(output)) stop("summary output already exists.")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
  manifest_path <- file.path(root, "manifest.csv")
  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
  expected <- roce_inference_pilot_manifest()
  roce_validate_inference_pilot_manifest(manifest, 1L)
  n <- suppressWarnings(as.numeric(args[[2L]]))
  checkpoints <- c(1L, 5L, 10L, 25L, 50L, 100L)
  if (length(n) != 1L || !is.finite(n) || n != floor(n) || !n %in% checkpoints) {
    stop("N must be one prespecified checkpoint: 1,5,10,25,50,100.")
  }
  n <- as.integer(n); tasks <- expected[seq_len(n), , drop = FALSE]
  manifest_hash <- roce_sha256_file(manifest_path)
  raw <- audits <- references <- vector("list", n)
  attempt_status <- data.frame(task_id = tasks$task_id, sim_id = tasks$sim_id,
                               status = "pending", detail = "", stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    seed <- tasks$sim_id[[i]]
    bundle <- file.path(root, sprintf("seed_%06d", seed))
    audit_directory <- file.path(root, "audits", sprintf("seed_%06d", seed))
    if (!dir.exists(bundle)) {
      attempt_status$status[[i]] <- "missing_bundle"; next
    }
    if (file.exists(file.path(bundle, "attempt_failure.csv"))) {
      attempt_status$status[[i]] <- "failed_attempt"
      attempt_status$detail[[i]] <- paste(readLines(file.path(bundle, "attempt_failure.csv")), collapse = " | ")
      next
    }
    if (!dir.exists(audit_directory)) {
      attempt_status$status[[i]] <- "missing_seed_audit"; next
    }
    raw_hash <- .siip_check_sha256(bundle, roce_sha256_file, c(
      "results.csv", "inference_audit.csv", "diagnostic_qc.csv",
      "artifacts.rds", "metadata.txt"
    ))
    audit_hash <- .siip_check_sha256(audit_directory, roce_sha256_file, c(
      "audit_passed.txt", "summary.csv", "numerical_checks.csv",
      "nuisance_diagnostics.csv"
    ))
    gate <- .siip_read_gate(file.path(audit_directory, "audit_passed.txt"))
    required_gate <- c("independent_inference_seed_audit", "bundle_directory",
      "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
      "bundle_checksum_manifest_fingerprint", "task_id", "sim_id")
    frozen_package <-
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7"
    frozen_workflow <-
      "9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4"
    if (!all(required_gate %in% names(gate)) ||
        gate[["independent_inference_seed_audit"]] != "passed" ||
        normalizePath(gate[["bundle_directory"]], mustWork = TRUE) != bundle ||
        gate[["bundle_checksum_manifest_fingerprint"]] != raw_hash ||
        gate[["package_fingerprint"]] != frozen_package ||
        gate[["workflow_fingerprint"]] != frozen_workflow ||
        gate[["manifest_fingerprint"]] != manifest_hash ||
        gate[["task_id"]] != as.character(tasks$task_id[[i]]) ||
        gate[["sim_id"]] != as.character(seed)) stop("seed audit gate mismatch: ", seed)
    rows <- read.csv(file.path(bundle, "results.csv"), stringsAsFactors = FALSE)
    inference <- read.csv(file.path(bundle, "inference_audit.csv"), stringsAsFactors = FALSE)
    regenerated <- roce_build_inference_audit(rows)
    if (!.siip_compare_frames(inference, regenerated)) stop("inference audit regeneration mismatch: ", seed)
    for (field in c("package_fingerprint", "workflow_fingerprint", "manifest_fingerprint")) {
      values <- unique(as.character(rows[[field]]))
      if (length(values) != 1L || values != gate[[field]]) stop("raw/audit provenance mismatch: ", seed)
    }
    if (any(rows$sim_id != seed) || any(inference$sim_id != seed)) stop("raw sim ID mismatch: ", seed)
    raw[[i]] <- rows; audits[[i]] <- inference
    references[[i]] <- data.frame(task_id = tasks$task_id[[i]], sim_id = seed,
      bundle_directory = bundle, bundle_checksum_manifest_hash = raw_hash,
      seed_audit_directory = audit_directory,
      seed_audit_checksum_manifest_hash = audit_hash, stringsAsFactors = FALSE)
    attempt_status$status[[i]] <- "complete"
  }
  if (any(attempt_status$status != "complete")) {
    print(attempt_status, row.names = FALSE)
    stop("checkpoint incomplete; no summary gate was published.")
  }
  raw_rows <- do.call(rbind, raw); inference_rows <- do.call(rbind, audits)
  summary <- roce_summarize_independent_inference(raw_rows, inference_rows)
  script_hash <- roce_sha256_file("diagnosis/tate_common_weight/summarize_independent_inference_pilot.R")
  roce_write_atomic_directory(output, function(stage) {
    write.csv(summary, file.path(stage, "rho_summary.csv"), row.names = FALSE, na = "")
    write.csv(raw_rows, file.path(stage, "all_raw_rows.csv"), row.names = FALSE, na = "")
    write.csv(inference_rows, file.path(stage, "all_inference_rows.csv"), row.names = FALSE, na = "")
    write.csv(do.call(rbind, references), file.path(stage, "input_references.csv"), row.names = FALSE)
    write.csv(attempt_status, file.path(stage, "attempt_status.csv"), row.names = FALSE)
    writeLines(c("independent_inference_summary=complete", paste0("n=", n),
      "inference_validated=FALSE", "statistical_review_required=TRUE",
      "warning_capture_complete=FALSE", paste0("manifest_fingerprint=", manifest_hash),
      paste0("summary_script_fingerprint=", script_hash)), file.path(stage, "metadata.txt"))
    files <- c("rho_summary.csv", "all_raw_rows.csv", "all_inference_rows.csv",
      "input_references.csv", "attempt_status.csv", "metadata.txt")
    writeLines(paste(vapply(file.path(stage, files), roce_sha256_file, character(1L)),
                     files, sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "independent inference checkpoint summary")
  invisible(summary)
}

if (sys.nframe() == 0L) summarize_independent_inference_pilot_main()
