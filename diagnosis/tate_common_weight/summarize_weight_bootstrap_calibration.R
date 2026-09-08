#!/usr/bin/env Rscript

.wbc_parse_seed_spec <- function(spec) {
  if (!is.character(spec) || length(spec) != 1L || is.na(spec) ||
      !nzchar(trimws(spec))) stop("SEED_IDS must be nonempty.", call. = FALSE)
  spec <- trimws(spec)
  maximum_seed_count <- 100000L
  if (grepl("^[0-9]+:[0-9]+$", spec)) {
    endpoints <- as.numeric(strsplit(spec, ":", fixed = TRUE)[[1L]])
    if (length(endpoints) != 2L || any(!is.finite(endpoints)) ||
        any(endpoints < 1) || any(endpoints != floor(endpoints)) ||
        any(endpoints > .Machine$integer.max) ||
        endpoints[[1L]] > endpoints[[2L]]) {
      stop("SEED_IDS range must contain increasing positive integers.",
           call. = FALSE)
    }
    requested_count <- endpoints[[2L]] - endpoints[[1L]] + 1
    if (!is.finite(requested_count) ||
        requested_count > maximum_seed_count) {
      stop("SEED_IDS may request at most ", maximum_seed_count,
           " seeds.", call. = FALSE)
    }
    values <- seq(endpoints[[1L]], endpoints[[2L]])
  } else if (grepl("^[0-9]+([,;][0-9]+)*$", spec)) {
    pieces <- strsplit(spec, "[,;]")[[1L]]
    if (length(pieces) > maximum_seed_count) {
      stop("SEED_IDS may request at most ", maximum_seed_count,
           " seeds.", call. = FALSE)
    }
    values <- as.numeric(pieces)
  } else {
    stop("SEED_IDS must use A:B or a comma/semicolon integer list.",
         call. = FALSE)
  }
  if (any(!is.finite(values)) || any(values < 1) ||
      any(values != floor(values)) || any(values > .Machine$integer.max) ||
      anyDuplicated(values)) {
    stop("SEED_IDS must be unique positive integers.", call. = FALSE)
  }
  as.integer(values)
}

.wbc_read_metadata <- function(path) {
  lines <- readLines(path, warn = FALSE)
  positions <- regexpr("=", lines, fixed = TRUE)
  if (length(lines) < 1L || any(positions < 2L))
    stop("malformed bundle metadata: ", path, call. = FALSE)
  keys <- substr(lines, 1L, positions - 1L)
  values <- substr(lines, positions + 1L, nchar(lines))
  if (anyDuplicated(keys) || any(!nzchar(keys)) || any(!nzchar(values)))
    stop("metadata keys/values must be unique and nonempty: ", path,
         call. = FALSE)
  stats::setNames(values, keys)
}

.wbc_strict_logical <- function(x, name) {
  if (is.logical(x) && !anyNA(x)) return(x)
  normalized <- toupper(trimws(as.character(x)))
  if (anyNA(normalized) || any(!normalized %in% c("TRUE", "FALSE")))
    stop(name, " must contain strict logical values.", call. = FALSE)
  normalized == "TRUE"
}

.wbc_verify_sha256 <- function(bundle, sha256_file) {
  manifest_path <- file.path(bundle, "sha256.txt")
  lines <- readLines(manifest_path, warn = FALSE)
  matched <- regexec("^([0-9A-Fa-f]{64})  ([^/]+)$", lines)
  pieces <- regmatches(lines, matched)
  if (length(lines) != 4L || any(lengths(pieces) != 3L))
    stop("invalid sha256.txt format: ", manifest_path, call. = FALSE)
  hashes <- tolower(vapply(pieces, `[[`, character(1L), 2L))
  files <- vapply(pieces, `[[`, character(1L), 3L)
  expected <- c("results.csv", "diagnostic_qc.csv", "artifacts.rds",
                "metadata.txt")
  if (anyDuplicated(files) || !setequal(files, expected))
    stop("sha256.txt does not name the exact bundle payload.", call. = FALSE)
  observed <- vapply(file.path(bundle, files), sha256_file, character(1L))
  if (!identical(unname(observed), unname(hashes)))
    stop("bundle SHA-256 verification failed: ", bundle, call. = FALSE)
  invisible(stats::setNames(hashes, files))
}

.wbc_scalar <- function(x, name) {
  values <- trimws(as.character(x))
  if (length(values) < 1L || anyNA(values) || any(!nzchar(values))) {
    stop(name, " must be complete and nonempty.", call. = FALSE)
  }
  values <- unique(values)
  if (length(values) != 1L) stop(name, " must have one value.", call. = FALSE)
  values[[1L]]
}

.wbc_audit_bundle <- function(bundle, seed, sha256_file) {
  required_files <- c("results.csv", "diagnostic_qc.csv", "artifacts.rds",
                      "metadata.txt", "sha256.txt")
  paths <- file.path(bundle, required_files)
  if (!dir.exists(bundle) || any(!file.exists(paths)) || any(dir.exists(paths)))
    stop("missing or incomplete expected seed bundle: ", bundle,
         call. = FALSE)
  payload_hashes <- .wbc_verify_sha256(bundle, sha256_file)
  metadata <- .wbc_read_metadata(file.path(bundle, "metadata.txt"))
  metadata_required <- c(
    "weight_bootstrap_calibration", "group_task_id", "sim_id",
    "rho_values", "n_weight_bootstrap", "package_library",
    "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
    "group_manifest_fingerprint"
  )
  if (length(setdiff(metadata_required, names(metadata))) > 0L ||
      metadata[["weight_bootstrap_calibration"]] != "complete" ||
      suppressWarnings(as.numeric(metadata[["sim_id"]])) != seed ||
      suppressWarnings(as.numeric(metadata[["group_task_id"]])) != seed) {
    stop("bundle metadata does not match expected seed ", seed, ".",
         call. = FALSE)
  }
  digest_names <- c("package_fingerprint", "workflow_fingerprint",
                    "manifest_fingerprint", "group_manifest_fingerprint")
  if (any(!grepl("^[0-9a-f]{64}$", tolower(metadata[digest_names]))))
    stop("bundle metadata contains an invalid provenance digest.",
         call. = FALSE)

  results <- read.csv(file.path(bundle, "results.csv"),
                      stringsAsFactors = FALSE, check.names = FALSE,
                      na.strings = c("", "NA"))
  qc <- read.csv(file.path(bundle, "diagnostic_qc.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE,
                 na.strings = c("", "NA"))
  expected_rhos <- c(0, 0.5, 1, 1.5, 2, 2.5)
  expected_base <- c("one_round_crossfit", "target_only", "sample_size",
                     "inverse_variance", "federated_dr", "pooled_dr")
  expected_methods <- c(expected_base, paste0(expected_base, "_ate"),
                        "one_round_crossfit_ate_armwise",
                        "one_round_crossfit_ate_hard_threshold")
  core <- c("sim_id", "rho", "method", "estimate", "se", "truth", "bias",
            "coverage", "ci_lower", "ci_upper", "n_weight_bootstrap",
            "package_fingerprint", "workflow_fingerprint",
            "manifest_fingerprint", "rho_group_manifest_fingerprint",
            "se_fixed_weight_bootstrap", "se_weight_relearn_bootstrap",
            "variance_weight_relearn_bootstrap", "weight_uncertainty_ratio",
            "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
            "weight_bootstrap_failures", "weight_bootstrap_multiplier",
            "weight_bootstrap_relearn_weights",
            "weight_bootstrap_screening_rule")
  missing <- setdiff(core, names(results))
  if (length(missing) > 0L) stop("results schema missing: ",
    paste(missing, collapse = ", "), call. = FALSE)
  if (nrow(results) != 84L || any(results$sim_id != seed) ||
      !identical(sort(unique(results$rho)), expected_rhos))
    stop("results have invalid row/seed/rho structure.", call. = FALSE)
  for (rho in expected_rhos) {
    methods <- results$method[results$rho == rho]
    if (anyDuplicated(methods) || !setequal(methods, expected_methods) ||
        length(methods) != length(expected_methods))
      stop("rho method schema mismatch for seed ", seed, ".", call. = FALSE)
  }
  numeric_core <- c("estimate", "se", "truth", "bias", "ci_lower", "ci_upper")
  if (any(vapply(results[numeric_core], function(x)
      !is.numeric(x) || any(!is.finite(x)), logical(1L))) ||
      any(results$se <= 0) || any(results$ci_lower > results$ci_upper) ||
      any(abs(results$bias - (results$estimate - results$truth)) > 1e-10))
    stop("invalid point-estimate values in seed ", seed, ".", call. = FALSE)
  coverage <- .wbc_strict_logical(results$coverage, "coverage")
  coverage_identity <- results$truth >= results$ci_lower &
    results$truth <= results$ci_upper
  if (!identical(coverage, coverage_identity))
    stop("coverage/CI identity failed for seed ", seed, ".", call. = FALSE)
  B <- suppressWarnings(as.numeric(metadata[["n_weight_bootstrap"]]))
  if (!is.finite(B) || B < 2 || B != floor(B) ||
      any(results$n_weight_bootstrap != B))
    stop("invalid bootstrap draw count in seed ", seed, ".", call. = FALSE)
  soft <- results$method == "one_round_crossfit_ate"
  bootstrap_numeric <- c("se_fixed_weight_bootstrap",
    "se_weight_relearn_bootstrap", "variance_weight_relearn_bootstrap",
    "weight_uncertainty_ratio", "weight_relearn_n_bootstrap",
    "weight_bootstrap_seed", "weight_bootstrap_failures")
  if (sum(soft) != 6L || any(vapply(results[soft, bootstrap_numeric],
      function(x) any(!is.finite(x)), logical(1L))) ||
      any(results$se_fixed_weight_bootstrap[soft] <= 0) ||
      any(results$se_weight_relearn_bootstrap[soft] <= 0) ||
      any(results$weight_relearn_n_bootstrap[soft] != B) ||
      any(results$weight_bootstrap_failures[soft] != 0) ||
      any(abs(results$variance_weight_relearn_bootstrap[soft] -
              results$se_weight_relearn_bootstrap[soft]^2) > 1e-10) ||
      any(abs(results$weight_uncertainty_ratio[soft] -
              results$se_weight_relearn_bootstrap[soft] /
              results$se_fixed_weight_bootstrap[soft]) > 1e-10) ||
      any(results$weight_bootstrap_multiplier[soft] !=
            "site_stratified_exponential") ||
      any(results$weight_bootstrap_screening_rule[soft] != "soft_penalty") ||
      any(!.wbc_strict_logical(results$weight_bootstrap_relearn_weights[soft],
                               "weight_bootstrap_relearn_weights")))
    stop("soft bootstrap schema failed for seed ", seed, ".", call. = FALSE)
  nonsoft_columns <- c(bootstrap_numeric, "weight_bootstrap_multiplier",
                       "weight_bootstrap_relearn_weights",
                       "weight_bootstrap_screening_rule")
  if (any(!is.na(unlist(results[!soft, nonsoft_columns, drop = FALSE],
                        use.names = FALSE))))
    stop("non-soft rows unexpectedly contain bootstrap diagnostics.",
         call. = FALSE)

  provenance_map <- c(
    package_fingerprint = "package_fingerprint",
    workflow_fingerprint = "workflow_fingerprint",
    manifest_fingerprint = "manifest_fingerprint",
    rho_group_manifest_fingerprint = "group_manifest_fingerprint"
  )
  for (column in names(provenance_map)) {
    if (.wbc_scalar(tolower(results[[column]]), column) !=
        tolower(metadata[[provenance_map[[column]]]]))
      stop("result/metadata provenance mismatch: ", column, call. = FALSE)
  }
  qc_required <- c("diagnostic_status", "implementation_failed",
                   "statistical_review_required")
  if (length(setdiff(qc_required, names(qc))) > 0L || nrow(qc) != 84L ||
      any(.wbc_strict_logical(qc$implementation_failed,
                              "implementation_failed")) ||
      any(grepl("invalid_weight_relearn_bootstrap_diagnostics|weight_relearn_bootstrap_failures_detected|nonfinite_values|bias_truth_inconsistent|coverage_ci_inconsistent",
                qc$diagnostic_status)))
    stop("stored diagnostic QC contains a schema/implementation failure.",
         call. = FALSE)

  list(
    results = results, qc = qc, metadata = metadata,
    payload_hashes = payload_hashes,
    checksum_manifest_sha256 = sha256_file(file.path(bundle, "sha256.txt")),
    B = as.integer(B)
  )
}

.wbc_mean_mcse <- function(x) {
  c(mean = mean(x), mcse = if (length(x) > 1L) stats::sd(x) /
      sqrt(length(x)) else NA_real_)
}

.wbc_summarize <- function(replicates) {
  split_rho <- split(replicates, replicates$rho)
  rows <- lapply(split_rho, function(x) {
    n <- nrow(x); z <- stats::qnorm(0.975)
    error <- x$soft_estimate - x$truth
    centered_error <- error - mean(error)
    cover_fixed <- abs(error) <= z * x$se_fixed
    cover_relearn <- abs(error) <= z * x$se_relearn
    cover_analytic_bc <- abs(centered_error) <= z * x$se_analytic
    cover_fixed_bc <- abs(centered_error) <= z * x$se_fixed
    cover_relearn_bc <- abs(centered_error) <= z * x$se_relearn
    empirical_sd <- if (n > 1L) stats::sd(x$soft_estimate) else NA_real_
    data.frame(
      rho = x$rho[[1L]], n_replications = n,
      tiny_n_checkpoint = n < 30L,
      mean_se_analytic = mean(x$se_analytic),
      mean_se_fixed_weight_bootstrap = mean(x$se_fixed),
      mean_se_weight_relearn_bootstrap = mean(x$se_relearn),
      mean_fixed_to_analytic_se = mean(x$se_fixed / x$se_analytic),
      mean_relearn_to_analytic_se = mean(x$se_relearn / x$se_analytic),
      mean_relearn_to_fixed_se = mean(x$se_relearn / x$se_fixed),
      empirical_sd_soft = empirical_sd,
      analytic_se_to_empirical_sd = mean(x$se_analytic) / empirical_sd,
      fixed_se_to_empirical_sd = mean(x$se_fixed) / empirical_sd,
      relearn_se_to_empirical_sd = mean(x$se_relearn) / empirical_sd,
      soft_bias = mean(error), soft_rmse = sqrt(mean(error^2)),
      target_rmse = sqrt(mean((x$target_estimate - x$truth)^2)),
      hard_rmse = sqrt(mean((x$hard_estimate - x$truth)^2)),
      coverage_analytic = mean(x$cover_analytic),
      coverage_fixed_weight = mean(cover_fixed),
      coverage_weight_relearn = mean(cover_relearn),
      coverage_analytic_bias_corrected = mean(cover_analytic_bc),
      coverage_fixed_bias_corrected = mean(cover_fixed_bc),
      coverage_relearn_bias_corrected = mean(cover_relearn_bc),
      coverage_relearn_mcse = sqrt(mean(cover_relearn) *
        (1 - mean(cover_relearn)) / n),
      coverage_target_analytic = mean(x$cover_target),
      coverage_hard_analytic = mean(x$cover_hard),
      paired_fixed_minus_analytic_coverage = mean(
        as.numeric(cover_fixed) - as.numeric(x$cover_analytic)
      ),
      paired_relearn_minus_analytic_coverage = mean(
        as.numeric(cover_relearn) - as.numeric(x$cover_analytic)
      ),
      paired_soft_minus_target_coverage = mean(
        as.numeric(x$cover_analytic) - as.numeric(x$cover_target)
      ),
      paired_soft_minus_hard_coverage = mean(
        as.numeric(x$cover_analytic) - as.numeric(x$cover_hard)
      ),
      paired_soft_minus_target_mse = mean(error^2 -
        (x$target_estimate - x$truth)^2),
      paired_soft_minus_h_mse = mean(error^2 -
        (x$hard_estimate - x$truth)^2),
      paired_relearn_minus_analytic_variance = mean(x$se_relearn^2 -
        x$se_analytic^2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_weight_bootstrap_calibration_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop(paste(
    "usage: summarize_weight_bootstrap_calibration.R",
    "BUNDLE_ROOT SEED_IDS OUTPUT_DIR"), call. = FALSE)
  bundle_root <- normalizePath(args[[1L]], mustWork = TRUE)
  seeds <- .wbc_parse_seed_spec(args[[2L]])
  output <- file.path(normalizePath(dirname(args[[3L]]), mustWork = FALSE),
                      basename(args[[3L]]))
  if (dir.exists(output) || file.exists(output))
    stop("summary destination already exists: ", output, call. = FALSE)
  project_root <- normalizePath(Sys.getenv("ROCE_PROJECT_ROOT", getwd()),
                                mustWork = TRUE)
  summary_script_path <- file.path(
    project_root, "diagnosis", "tate_common_weight",
    "summarize_weight_bootstrap_calibration.R"
  )
  provenance_environment <- new.env(parent = baseenv())
  sys.source(file.path(project_root, "scripts/slurm/result_provenance.R"),
             provenance_environment)
  atomic_environment <- new.env(parent = baseenv())
  sys.source(file.path(project_root, "scripts/slurm/atomic_output.R"),
             atomic_environment)
  summary_script_fingerprint <-
    provenance_environment$roce_sha256_file(summary_script_path)
  audited <- lapply(seeds, function(seed) .wbc_audit_bundle(
    file.path(bundle_root, sprintf("seed_%06d", seed)), seed,
    provenance_environment$roce_sha256_file
  ))
  provenance_fields <- c("package_library", "package_fingerprint",
    "workflow_fingerprint", "manifest_fingerprint",
    "group_manifest_fingerprint", "n_weight_bootstrap")
  for (field in provenance_fields) {
    values <- vapply(audited, function(x) x$metadata[[field]], character(1L))
    if (length(unique(values)) != 1L)
      stop("bundles disagree on provenance/config field: ", field,
           call. = FALSE)
  }
  replicate_rows <- do.call(rbind, lapply(seq_along(audited), function(i) {
    d <- audited[[i]]$results
    soft <- d[d$method == "one_round_crossfit_ate", , drop = FALSE]
    target <- d[d$method == "target_only_ate", , drop = FALSE]
    hard <- d[d$method == "one_round_crossfit_ate_hard_threshold", , drop = FALSE]
    target <- target[match(soft$rho, target$rho), , drop = FALSE]
    hard <- hard[match(soft$rho, hard$rho), , drop = FALSE]
    data.frame(
      sim_id = seeds[[i]], rho = soft$rho, truth = soft$truth,
      soft_estimate = soft$estimate, target_estimate = target$estimate,
      hard_estimate = hard$estimate, se_analytic = soft$se,
      se_fixed = soft$se_fixed_weight_bootstrap,
      se_relearn = soft$se_weight_relearn_bootstrap,
      uncertainty_ratio = soft$weight_uncertainty_ratio,
      cover_analytic = .wbc_strict_logical(soft$coverage, "coverage"),
      cover_target = .wbc_strict_logical(target$coverage, "coverage"),
      cover_hard = .wbc_strict_logical(hard$coverage, "coverage"),
      stringsAsFactors = FALSE
    )
  }))
  summary <- .wbc_summarize(replicate_rows)
  audit <- data.frame(
    sim_id = seeds,
    bundle = normalizePath(file.path(bundle_root,
      sprintf("seed_%06d", seeds)), mustWork = TRUE),
    sha256_verified = TRUE, schema_verified = TRUE,
    implementation_qc_passed = TRUE,
    package_fingerprint = vapply(audited, function(x)
      x$metadata[["package_fingerprint"]], character(1L)),
    workflow_fingerprint = vapply(audited, function(x)
      x$metadata[["workflow_fingerprint"]], character(1L)),
    checksum_manifest_sha256 = vapply(audited, function(x)
      x$checksum_manifest_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  writer <- function(stage) {
    utils::write.csv(audit, file.path(stage, "bundle_audit.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(replicate_rows,
      file.path(stage, "replicate_diagnostics.csv"), row.names = FALSE,
      na = "")
    utils::write.csv(summary, file.path(stage, "rho_summary.csv"),
                     row.names = FALSE, na = "")
    writeLines(c(
      "weight_bootstrap_calibration_summary=complete",
      paste0("bundle_root=", bundle_root),
      paste0("seed_ids=", paste(seeds, collapse = ";")),
      paste0("n_replications=", length(seeds)),
      paste0("n_weight_bootstrap=",
             audited[[1L]]$metadata[["n_weight_bootstrap"]]),
      paste0("package_fingerprint=",
             audited[[1L]]$metadata[["package_fingerprint"]]),
      paste0("workflow_fingerprint=",
             audited[[1L]]$metadata[["workflow_fingerprint"]]),
      paste0("summary_script_fingerprint=", summary_script_fingerprint),
      paste0(
        "input_checksum_manifest_sha256=",
        paste(
          paste0(seeds, ":", vapply(audited, function(x)
            x$checksum_manifest_sha256, character(1L))),
          collapse = ";"
        )
      ),
      "interpretation=tiny-n checkpoint; coverage and empirical-SD estimates are descriptive and must not be treated as stable Monte Carlo conclusions"
    ), file.path(stage, "metadata.txt"), useBytes = TRUE)
    payloads <- file.path(stage, c("bundle_audit.csv",
      "replicate_diagnostics.csv", "rho_summary.csv", "metadata.txt"))
    hashes <- vapply(payloads, provenance_environment$roce_sha256_file,
                     character(1L))
    writeLines(paste(hashes, basename(payloads), sep = "  "),
               file.path(stage, "sha256.txt"), useBytes = TRUE)
  }
  atomic_environment$roce_write_atomic_directory(
    output, writer, caller = "weight-bootstrap calibration summary"
  )
  message("[done] wrote calibration summary: ", output)
  invisible(output)
}

if (sys.nframe() == 0L) summarize_weight_bootstrap_calibration_main()
