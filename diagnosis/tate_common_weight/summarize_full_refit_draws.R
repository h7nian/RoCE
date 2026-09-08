#!/usr/bin/env Rscript

source("diagnosis/tate_common_weight/audit_full_refit_intermediates.R")

.frs_parse_draw_ids <- function(spec) {
  if (!is.character(spec) || length(spec) != 1L || is.na(spec) ||
      !nzchar(trimws(spec))) stop("DRAW_IDS must be nonempty.", call. = FALSE)
  spec <- trimws(spec)
  if (grepl("^[0-9]+:[0-9]+$", spec)) {
    endpoints <- as.numeric(strsplit(spec, ":", fixed = TRUE)[[1L]])
    if (any(!is.finite(endpoints)) || any(endpoints != floor(endpoints)) ||
        endpoints[[1L]] > endpoints[[2L]] || endpoints[[1L]] < 0 ||
        endpoints[[2L]] > 1000 ||
        endpoints[[2L]] - endpoints[[1L]] + 1 > 1001) {
      stop("DRAW_IDS range must be increasing integers in [0, 1000].",
           call. = FALSE)
    }
    values <- seq(endpoints[[1L]], endpoints[[2L]])
  } else if (grepl("^[0-9]+([,;][0-9]+)*$", spec)) {
    pieces <- strsplit(spec, "[,;]")[[1L]]
    if (length(pieces) > 1001L)
      stop("DRAW_IDS may contain at most 1001 IDs.", call. = FALSE)
    values <- as.numeric(pieces)
  } else {
    stop("DRAW_IDS must use A:B or a comma/semicolon integer list.",
         call. = FALSE)
  }
  if (any(!is.finite(values)) || any(values < 0) || any(values > 1000) ||
      any(values != floor(values)) || anyDuplicated(values))
    stop("DRAW_IDS must be unique integers in [0, 1000].", call. = FALSE)
  values <- as.integer(values)
  if (!0L %in% values)
    stop("DRAW_IDS must explicitly include identity draw 0.", call. = FALSE)
  values
}

.frs_integer <- function(x, name, lower = 0L) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < lower || value != floor(value) ||
      value > .Machine$integer.max)
    stop(name, " must be one integer >= ", lower, ".", call. = FALSE)
  as.integer(value)
}

.frs_rho <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  allowed <- c(0, 0.5, 1, 1.5, 2, 2.5)
  if (length(value) != 1L || is.na(value) || !value %in% allowed)
    stop("RHO must be one of 0,0.5,1,1.5,2,2.5.", call. = FALSE)
  value
}

.frs_logical <- function(x, name) {
  if (is.logical(x) && length(x) == 1L && !is.na(x)) return(x)
  normalized <- toupper(trimws(as.character(x)))
  if (length(normalized) != 1L || is.na(normalized) ||
      !normalized %in% c("TRUE", "FALSE"))
    stop(name, " must be one strict logical.", call. = FALSE)
  normalized == "TRUE"
}

.frs_verify_bundle <- function(bundle, sha256_file) {
  checksum_path <- file.path(bundle, "sha256.txt")
  lines <- readLines(checksum_path, warn = FALSE)
  matches <- regexec("^([0-9a-f]{64})  ([^/]+)$", lines)
  pieces <- regmatches(lines, matches)
  if (length(lines) != 2L || any(lengths(pieces) != 3L))
    stop("malformed draw checksum manifest: ", bundle, call. = FALSE)
  hashes <- vapply(pieces, `[[`, character(1L), 2L)
  files <- vapply(pieces, `[[`, character(1L), 3L)
  if (anyDuplicated(files) || !setequal(files, c("result.csv", "draw.rds")))
    stop("draw checksum manifest must name result.csv and draw.rds exactly.",
         call. = FALSE)
  observed <- vapply(file.path(bundle, files), sha256_file, character(1L))
  if (!identical(unname(observed), unname(hashes)))
    stop("draw payload checksum mismatch: ", bundle, call. = FALSE)
  list(
    payload_hashes = stats::setNames(hashes, files),
    checksum_manifest_sha256 = sha256_file(checksum_path)
  )
}

.frs_read_candidate <- function(bundle) {
  path <- file.path(bundle, "result.csv")
  if (!file.exists(path)) return(NULL)
  result <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                     na.strings = c("", "NA"))
  required <- c("sim_id", "rho", "draw_id")
  if (nrow(result) != 1L || length(setdiff(required, names(result))) > 0L)
    stop("candidate draw result has invalid identifying schema: ", bundle,
         call. = FALSE)
  list(bundle = bundle, result = result)
}

.frs_compare_saved_summary <- function(csv, saved, bundle) {
  if (!is.data.frame(saved) || nrow(saved) != 1L ||
      length(setdiff(names(csv), names(saved))) > 0L) {
    stop("draw.rds summary schema disagrees with result.csv: ", bundle,
         call. = FALSE)
  }
  for (name in names(csv)) {
    left <- csv[[name]]
    right <- saved[[name]]
    if (is.numeric(left)) {
      left <- as.numeric(left)
      right <- suppressWarnings(as.numeric(right))
      equal <- length(right) == 1L &&
        identical(is.na(left), is.na(right)) &&
        (is.na(left) || isTRUE(all.equal(
          left, right, tolerance = 1e-12, check.attributes = FALSE
        )))
    } else if (is.logical(left)) {
      equal <- identical(left, as.logical(right))
    } else {
      normalize <- function(x) {
        x <- as.character(x)
        x[is.na(x)] <- ""
        x
      }
      equal <- identical(normalize(left), normalize(right))
    }
    if (!isTRUE(equal)) {
      stop("draw.rds summary disagrees with result.csv field '", name,
           "': ", bundle, call. = FALSE)
    }
  }
  invisible(TRUE)
}

.frs_validate_selected <- function(selected, sim_id, rho, draw_ids,
                                   sha256_file) {
  required <- c(
    "sim_id", "rho", "draw_id", "resample_seed", "identity", "status",
    "failure_message", "failure_stage",
    "estimate_refit_relearned_weights", "estimate_refit_original_weights",
    "estimate_fixed_nuisance_original_weights",
    "estimate_fixed_nuisance_relearned_weights", "estimate_target_refit",
    "estimate_reference", "se_analytic_reference", "nuisance_refit",
    "inference_validated", "conditional_on_saved_design_and_outer_partition",
    "nuisance_cv_groups_duplicate_origins",
    "nuisance_cv_duplicate_origin_leakage_possible", "resampling_scheme",
    "package_fingerprint", "workflow_fingerprint",
    "source_bundle_fingerprint"
  )
  rows <- vector("list", length(selected))
  audit <- vector("list", length(selected))
  for (i in seq_along(selected)) {
    entry <- selected[[i]]
    hash <- .frs_verify_bundle(entry$bundle, sha256_file)
    row <- entry$result
    saved <- readRDS(file.path(entry$bundle, "draw.rds"))
    if (!is.list(saved) || !is.list(saved$result) ||
        is.null(saved$summary) ||
        !is.character(saved$source_bundle) ||
        length(saved$source_bundle) != 1L || is.na(saved$source_bundle) ||
        !dir.exists(saved$source_bundle) ||
        !file.exists(file.path(saved$source_bundle, "sha256.txt"))) {
      stop("draw.rds is incomplete or has an invalid source bundle: ",
           entry$bundle, call. = FALSE)
    }
    .frs_compare_saved_summary(row, saved$summary, entry$bundle)
    missing <- setdiff(required, names(row))
    if (length(missing) > 0L)
      stop("draw result schema missing: ", paste(missing, collapse = ", "),
           call. = FALSE)
    observed_sim <- .frs_integer(row$sim_id, "sim_id", 1L)
    observed_draw <- .frs_integer(row$draw_id, "draw_id", 0L)
    if (observed_sim != sim_id || as.numeric(row$rho) != rho ||
        !observed_draw %in% draw_ids)
      stop("selected draw identifiers changed during validation.",
           call. = FALSE)
    identity <- .frs_logical(row$identity, "identity")
    if (!identical(identity, observed_draw == 0L))
      stop("identity flag must equal draw_id == 0.", call. = FALSE)
    status <- as.character(row$status)
    if (length(status) != 1L || is.na(status) ||
        !status %in% c("completed", "failed"))
      stop("draw status must be completed or failed.", call. = FALSE)
    succeeded <- status == "completed"
    grouped_flag <- .frs_logical(
      row$nuisance_cv_groups_duplicate_origins,
      "nuisance_cv_groups_duplicate_origins"
    )
    leakage_flag <- .frs_logical(
      row$nuisance_cv_duplicate_origin_leakage_possible,
      "nuisance_cv_duplicate_origin_leakage_possible"
    )
    protocol <- if ("nuisance_cv_partition" %in% names(row)) {
      as.character(row$nuisance_cv_partition)
    } else if (!grouped_flag && leakage_flag) {
      "row_level"
    } else {
      stop("draw lacks an explicit valid nuisance CV partition protocol.",
           call. = FALSE)
    }
    if (identical(protocol, "origin_grouped")) {
      if (!grouped_flag || leakage_flag) {
        stop("origin-grouped protocol booleans are inconsistent.",
             call. = FALSE)
      }
    } else if (!identical(protocol, "row_level") || grouped_flag ||
               !leakage_flag) {
      stop("archived row-level protocol booleans are inconsistent.",
           call. = FALSE)
    }
    estimate_columns <- c(
      "estimate_refit_relearned_weights", "estimate_refit_original_weights",
      "estimate_fixed_nuisance_original_weights",
      "estimate_fixed_nuisance_relearned_weights", "estimate_target_refit"
    )
    values <- suppressWarnings(as.numeric(unlist(row[estimate_columns],
                                                 use.names = FALSE)))
    if (succeeded && (any(!is.finite(values)) ||
        is.na(row$failure_message) == FALSE ||
        is.na(row$failure_stage) == FALSE))
      stop("completed draw has nonfinite estimates or a failure record.",
           call. = FALSE)
    if (!succeeded && (is.na(row$failure_message) ||
        !nzchar(trimws(row$failure_message)) || is.na(row$failure_stage) ||
        !nzchar(trimws(row$failure_stage))))
      stop("failed draw must retain stage and message.", call. = FALSE)
    reference_se <- suppressWarnings(as.numeric(row$se_analytic_reference))
    reference_estimate <- suppressWarnings(as.numeric(row$estimate_reference))
    if (length(reference_se) != 1L || !is.finite(reference_se) ||
        reference_se <= 0 || length(reference_estimate) != 1L ||
        !is.finite(reference_estimate)) {
      stop("reference estimate/analytic SE must be finite with positive SE.",
           call. = FALSE)
    }
    result <- saved$result
    if (succeeded) {
      completed_scalars <- c(
        estimate_refit_relearned_weights =
          result$estimate_refit_relearned_weights,
        estimate_refit_original_weights =
          result$estimate_refit_original_weights,
        estimate_fixed_nuisance_original_weights =
          result$matched_fixed_nuisance$estimate_original_weights,
        estimate_fixed_nuisance_relearned_weights =
          result$matched_fixed_nuisance$estimate_relearned_weights,
        estimate_target_refit = result$fitted$target_only$estimate
      )
      fitted_se <- suppressWarnings(as.numeric(result$fitted$se))
      if (length(completed_scalars) != length(estimate_columns) ||
          any(!is.finite(completed_scalars)) ||
          max(abs(completed_scalars - values)) > 1e-12 ||
          abs(as.numeric(result$fitted$estimate) -
              row$estimate_refit_relearned_weights) > 1e-12 ||
          length(fitted_se) != 1L || !is.finite(fitted_se) || fitted_se <= 0) {
        stop("completed draw.rds scalars disagree with result.csv.",
             call. = FALSE)
      }
      if (identity && abs(fitted_se - reference_se) > 1e-10) {
        stop("identity draw fitted SE does not reproduce the reference SE.",
             call. = FALSE)
      }
    } else {
      if (!identical(as.character(result$failure_message),
                     as.character(row$failure_message)) ||
          !identical(as.character(result$failure_stage),
                     as.character(row$failure_stage))) {
        stop("failed draw.rds record disagrees with result.csv.",
             call. = FALSE)
      }
      fixed_csv <- as.numeric(row[c(
        "estimate_fixed_nuisance_original_weights",
        "estimate_fixed_nuisance_relearned_weights"
      )])
      if (is.null(result$matched_fixed_nuisance)) {
        if (any(!is.na(fixed_csv))) {
          stop("failed draw CSV retains unmatched fixed-nuisance estimates.",
               call. = FALSE)
        }
      } else {
        fixed_rds <- c(
          result$matched_fixed_nuisance$estimate_original_weights,
          result$matched_fixed_nuisance$estimate_relearned_weights
        )
        if (length(fixed_rds) != 2L || any(!is.finite(fixed_rds)) ||
            any(!is.finite(fixed_csv)) ||
            max(abs(fixed_rds - fixed_csv)) > 1e-12) {
          stop("failed draw fixed-nuisance scalars disagree with draw.rds.",
               call. = FALSE)
        }
      }
    }
    source_bundle_hash <- sha256_file(file.path(
      saved$source_bundle, "sha256.txt"
    ))
    if (!identical(source_bundle_hash,
                   as.character(row$source_bundle_fingerprint))) {
      stop("draw source-bundle checksum disagrees with result provenance.",
           call. = FALSE)
    }
    partition_audit <- if (succeeded) {
      roce_audit_nuisance_cv_partitions(result, protocol)
    } else {
      list(protocol = protocol, n_records = 0L, n_explicit = 0L)
    }
    if (identical(protocol, "origin_grouped")) {
      if (!"source_package_fingerprint" %in% names(row) ||
          !grepl("^[0-9a-f]{64}$",
                 as.character(row$source_package_fingerprint))) {
        stop("origin-grouped draw lacks a source-package fingerprint.",
             call. = FALSE)
      }
      source_metadata_path <- file.path(saved$source_bundle, "metadata.txt")
      if (!file.exists(source_metadata_path)) {
        stop("origin-grouped source bundle lacks metadata.txt.",
             call. = FALSE)
      }
      source_manifest_lines <- readLines(file.path(
        saved$source_bundle, "sha256.txt"
      ), warn = FALSE)
      metadata_manifest_line <- grep(
        "^[0-9a-f]{64}  metadata[.]txt$", source_manifest_lines,
        value = TRUE
      )
      if (length(metadata_manifest_line) != 1L ||
          !identical(
            substr(metadata_manifest_line, 1L, 64L),
            sha256_file(source_metadata_path)
          )) {
        stop("origin-grouped source metadata checksum is invalid.",
             call. = FALSE)
      }
      metadata_lines <- readLines(source_metadata_path, warn = FALSE)
      metadata <- stats::setNames(
        sub("^[^=]*=", "", metadata_lines), sub("=.*$", "", metadata_lines)
      )
      if (!identical(
        unname(metadata["package_fingerprint"]),
        as.character(row$source_package_fingerprint)
      )) {
        stop("source-package fingerprint disagrees with source metadata.",
             call. = FALSE)
      }
    }
    expected_seed <- 190000L + sim_id * 10000L + observed_draw
    if (.frs_integer(row$resample_seed, "resample_seed", 1L) != expected_seed)
      stop("draw resample_seed does not match the frozen mapping.",
           call. = FALSE)
    logical_expectations <- c(
      nuisance_refit = TRUE, inference_validated = FALSE,
      conditional_on_saved_design_and_outer_partition = TRUE,
      nuisance_cv_groups_duplicate_origins = grouped_flag,
      nuisance_cv_duplicate_origin_leakage_possible = leakage_flag
    )
    for (name in names(logical_expectations)) {
      if (!identical(.frs_logical(row[[name]], name),
                     logical_expectations[[name]]))
        stop("draw logical contract mismatch: ", name, call. = FALSE)
    }
    if (!identical(as.character(row$resampling_scheme),
                   "site_by_original_outer_fold_multinomial"))
      stop("draw resampling scheme mismatch.", call. = FALSE)
    digests <- c("package_fingerprint", "workflow_fingerprint",
                 "source_bundle_fingerprint")
    if (any(!grepl("^[0-9a-f]{64}$", unlist(row[digests],
                                             use.names = FALSE))))
      stop("draw provenance digest is invalid.", call. = FALSE)
    rows[[i]] <- row
    audit[[i]] <- data.frame(
      draw_id = observed_draw, bundle = normalizePath(entry$bundle),
      sha256_verified = TRUE, rds_crosscheck_verified = TRUE,
      nuisance_cv_partition = protocol,
      nuisance_cv_partition_records = partition_audit$n_records,
      nuisance_cv_explicit_partition_records = partition_audit$n_explicit,
      checksum_manifest_sha256 = hash$checksum_manifest_sha256,
      stringsAsFactors = FALSE
    )
  }
  audit <- do.call(rbind, audit)
  if (length(unique(audit$nuisance_cv_partition)) != 1L) {
    stop("draw bundles mix nuisance CV partition protocols.", call. = FALSE)
  }
  rows <- do.call(rbind, rows)
  rows <- rows[order(rows$draw_id), , drop = FALSE]
  audit <- audit[match(rows$draw_id, audit$draw_id), , drop = FALSE]
  attr(rows, "nuisance_cv_partition") <- audit$nuisance_cv_partition[[1L]]
  provenance <- c("package_fingerprint", "workflow_fingerprint",
                  "source_bundle_fingerprint", "estimate_reference",
                  "se_analytic_reference", "resampling_scheme")
  if (identical(audit$nuisance_cv_partition[[1L]], "origin_grouped")) {
    provenance <- c(provenance, "source_package_fingerprint",
                    "nuisance_cv_partition")
  }
  for (name in provenance) {
    values <- rows[[name]]
    if (anyNA(values) || length(unique(as.character(values))) != 1L)
      stop("draw bundles mix provenance/reference field: ", name,
           call. = FALSE)
  }
  identity <- rows[rows$draw_id == 0L, , drop = FALSE]
  identity_estimates <- as.numeric(identity[1L, c(
    "estimate_refit_relearned_weights", "estimate_refit_original_weights",
    "estimate_fixed_nuisance_original_weights",
    "estimate_fixed_nuisance_relearned_weights"
  )])
  if (nrow(identity) != 1L || identity$status != "completed" ||
      any(!is.finite(identity_estimates)) ||
      max(abs(identity_estimates - identity$estimate_reference)) > 1e-10)
    stop("draw 0 did not pass the same-seed/rho identity gate.",
         call. = FALSE)
  list(rows = rows, audit = audit)
}

.frs_covariance_decomposition <- function(success) {
  A <- success$estimate_fixed_nuisance_original_weights
  B <- success$estimate_fixed_nuisance_relearned_weights
  C <- success$estimate_refit_original_weights
  D <- success$estimate_refit_relearned_weights
  weight <- B - A
  nuisance <- C - A
  interaction <- D - C - B + A
  total <- D - A
  variance <- function(x) if (length(x) > 1L) stats::var(x) else NA_real_
  covariance <- function(x, y) if (length(x) > 1L) stats::cov(x, y) else NA_real_
  change_terms <- c(
    variance_weight_change = variance(weight),
    variance_nuisance_change = variance(nuisance),
    variance_interaction = variance(interaction),
    twice_cov_weight_nuisance = 2 * covariance(weight, nuisance),
    twice_cov_weight_interaction = 2 * covariance(weight, interaction),
    twice_cov_nuisance_interaction = 2 * covariance(nuisance, interaction)
  )
  reconstructed_change <- if (all(is.finite(change_terms))) {
    sum(change_terms)
  } else NA_real_
  observed_change <- variance(total)
  change <- data.frame(
    decomposition_scope = "change_D_minus_A",
    term = c(
      names(change_terms), "reconstructed_variance_D_minus_A",
      "observed_variance_D_minus_A", "change_decomposition_error"
    ),
    value = c(
      change_terms, reconstructed_change, observed_change,
      abs(reconstructed_change - observed_change)
    ),
    stringsAsFactors = FALSE
  )

  full_terms <- c(
    variance_baseline_A = variance(A),
    variance_weight_change = variance(weight),
    variance_nuisance_change = variance(nuisance),
    variance_interaction = variance(interaction),
    twice_cov_A_weight = 2 * covariance(A, weight),
    twice_cov_A_nuisance = 2 * covariance(A, nuisance),
    twice_cov_A_interaction = 2 * covariance(A, interaction),
    twice_cov_weight_nuisance = 2 * covariance(weight, nuisance),
    twice_cov_weight_interaction = 2 * covariance(weight, interaction),
    twice_cov_nuisance_interaction = 2 * covariance(nuisance, interaction)
  )
  reconstructed_full <- if (all(is.finite(full_terms))) {
    sum(full_terms)
  } else NA_real_
  observed_full <- variance(D)
  full <- data.frame(
    decomposition_scope = "full_D",
    term = c(
      names(full_terms), "reconstructed_variance_D",
      "observed_variance_D", "full_decomposition_error"
    ),
    value = c(
      full_terms, reconstructed_full, observed_full,
      abs(reconstructed_full - observed_full)
    ),
    stringsAsFactors = FALSE
  )
  rbind(change, full)
}

.frs_summary <- function(rows) {
  attempted <- rows[rows$draw_id != 0L, , drop = FALSE]
  success <- attempted[attempted$status == "completed", , drop = FALSE]
  estimate_columns <- c(
    fixed_nuisance_original_weights =
      "estimate_fixed_nuisance_original_weights",
    fixed_nuisance_relearned_weights =
      "estimate_fixed_nuisance_relearned_weights",
    refit_nuisance_original_weights = "estimate_refit_original_weights",
    refit_nuisance_relearned_weights = "estimate_refit_relearned_weights"
  )
  sd_values <- vapply(estimate_columns, function(name) {
    if (nrow(success) > 1L) stats::sd(success[[name]]) else NA_real_
  }, numeric(1L))
  A <- success$estimate_fixed_nuisance_original_weights
  B <- success$estimate_fixed_nuisance_relearned_weights
  C <- success$estimate_refit_original_weights
  D <- success$estimate_refit_relearned_weights
  safe_mean <- function(x) if (length(x) > 0L) mean(x) else NA_real_
  paired_identity_error <- if (length(A) > 0L) {
    max(abs((D - A) - ((B - A) + (C - A) + (D - C - B + A))))
  } else {
    NA_real_
  }
  summary <- data.frame(
    sim_id = rows$sim_id[[1L]], rho = rows$rho[[1L]],
    n_attempted_nonidentity = nrow(attempted),
    n_successful_nonidentity = nrow(success),
    n_failed_nonidentity = sum(attempted$status == "failed"),
    identity_point_and_se_gate_passed = TRUE,
    diagnostic_draw_set_complete =
      nrow(attempted) >= 2L && all(attempted$status == "completed"),
    inference_validated = FALSE,
    small_B_warning = nrow(success) < 200L,
    nuisance_cv_partition = if (is.null(attr(
      rows, "nuisance_cv_partition"
    ))) "row_level" else attr(rows, "nuisance_cv_partition"),
    sd_fixed_nuisance_original_weights = sd_values[[1L]],
    sd_fixed_nuisance_relearned_weights = sd_values[[2L]],
    sd_refit_nuisance_original_weights = sd_values[[3L]],
    sd_refit_nuisance_relearned_weights = sd_values[[4L]],
    mean_weight_change_fixed_nuisance = safe_mean(B - A),
    mean_nuisance_change_original_weights = safe_mean(C - A),
    mean_total_change = safe_mean(D - A),
    mean_weight_nuisance_interaction = safe_mean(D - C - B + A),
    paired_identity_error = paired_identity_error,
    interpretation = paste0(
      "diagnostic only; draw 0 excluded from SD; no coverage claim; ",
      "small-B Monte Carlo variability may be material"
    ),
    stringsAsFactors = FALSE
  )
  list(summary = summary, success = success,
       covariance = .frs_covariance_decomposition(success))
}

summarize_full_refit_draws_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) stop(paste(
    "usage: summarize_full_refit_draws.R",
    "BUNDLE_ROOT SIM_ID RHO DRAW_IDS OUTPUT_DIR"), call. = FALSE)
  root <- normalizePath(args[[1L]], mustWork = TRUE)
  sim_id <- .frs_integer(args[[2L]], "SIM_ID", 1L)
  rho <- .frs_rho(args[[3L]])
  draw_ids <- .frs_parse_draw_ids(args[[4L]])
  output <- file.path(normalizePath(dirname(args[[5L]]), mustWork = FALSE),
                      basename(args[[5L]]))
  if (file.exists(output) || dir.exists(output))
    stop("summary output already exists: ", output, call. = FALSE)
  project_root <- normalizePath(Sys.getenv("ROCE_PROJECT_ROOT", getwd()),
                                mustWork = TRUE)
  provenance <- new.env(parent = baseenv())
  sys.source(file.path(project_root, "scripts/slurm/result_provenance.R"),
             provenance)
  atomic <- new.env(parent = baseenv())
  sys.source(file.path(project_root, "scripts/slurm/atomic_output.R"), atomic)
  script_path <- file.path(project_root, "diagnosis", "tate_common_weight",
                           "summarize_full_refit_draws.R")
  script_hash <- provenance$roce_sha256_file(script_path)
  audit_dependency_path <- file.path(
    project_root, "diagnosis", "tate_common_weight",
    "audit_full_refit_intermediates.R"
  )
  partition_validator_path <- file.path(
    project_root, "diagnosis", "tate_common_weight",
    "nuisance_cv_partition_audit.R"
  )
  audit_dependency_hash <- provenance$roce_sha256_file(audit_dependency_path)
  partition_validator_hash <- provenance$roce_sha256_file(
    partition_validator_path
  )

  directories <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  directories <- directories[!grepl("\\.lock$", directories)]
  candidates <- Filter(Negate(is.null), lapply(directories,
                                               .frs_read_candidate))
  selected <- Filter(function(x) {
    row <- x$result
    suppressWarnings(as.numeric(row$sim_id)) == sim_id &&
      suppressWarnings(as.numeric(row$rho)) == rho &&
      suppressWarnings(as.numeric(row$draw_id)) %in% draw_ids
  }, candidates)
  observed_ids <- vapply(selected, function(x)
    suppressWarnings(as.numeric(x$result$draw_id)), numeric(1L))
  if (length(selected) != length(draw_ids) || anyDuplicated(observed_ids) ||
      !setequal(observed_ids, draw_ids))
    stop("expected draw IDs are missing or duplicated.", call. = FALSE)
  validated <- .frs_validate_selected(
    selected, sim_id, rho, draw_ids, provenance$roce_sha256_file
  )
  summarized <- .frs_summary(validated$rows)
  draw_audit <- cbind(
    validated$rows,
    validated$audit[, setdiff(names(validated$audit), "draw_id"), drop = FALSE]
  )
  writer <- function(stage) {
    if (!identical(script_hash, provenance$roce_sha256_file(script_path)) ||
        !identical(audit_dependency_hash, provenance$roce_sha256_file(
          audit_dependency_path
        )) ||
        !identical(partition_validator_hash, provenance$roce_sha256_file(
          partition_validator_path
        ))) {
      stop("summary or partition-validation dependency changed before atomic write.",
           call. = FALSE)
    }
    utils::write.csv(draw_audit, file.path(stage, "draw_audit.csv"),
                     row.names = FALSE, na = "")
    utils::write.csv(summarized$summary,
      file.path(stage, "resampling_summary.csv"), row.names = FALSE, na = "")
    utils::write.csv(summarized$covariance,
      file.path(stage, "covariance_decomposition.csv"), row.names = FALSE,
      na = "")
    writeLines(c(
      "full_refit_draw_summary=complete",
      paste0("bundle_root=", root), paste0("sim_id=", sim_id),
      paste0("rho=", rho),
      paste0("draw_ids=", paste(draw_ids, collapse = ";")),
      paste0("summary_script_fingerprint=", script_hash),
      paste0("audit_dependency_fingerprint=", audit_dependency_hash),
      paste0("partition_validator_fingerprint=", partition_validator_hash),
      paste0("input_checksum_manifest_sha256=", paste(
        paste0(draw_audit$draw_id, ":", draw_audit$checksum_manifest_sha256),
        collapse = ";")),
      paste0("diagnostic_draw_set_complete=",
             summarized$summary$diagnostic_draw_set_complete),
      paste0("nuisance_cv_partition=",
             summarized$summary$nuisance_cv_partition),
      "inference_validated=FALSE",
      "identity_gate_scope=point estimates and analytic SE; fold-weight identity is validated by the draw runner and is not independently recomputed here",
      "interpretation=diagnostic small-B resampling summary only; not a validated interval or coverage claim"
    ), file.path(stage, "metadata.txt"), useBytes = TRUE)
    payloads <- file.path(stage, c("draw_audit.csv",
      "resampling_summary.csv", "covariance_decomposition.csv",
      "metadata.txt"))
    hashes <- vapply(payloads, provenance$roce_sha256_file, character(1L))
    writeLines(paste(hashes, basename(payloads), sep = "  "),
               file.path(stage, "sha256.txt"), useBytes = TRUE)
  }
  atomic$roce_write_atomic_directory(
    output, writer, caller = "full-refit draw summary"
  )
  message("[done] wrote full-refit draw summary: ", output)
  invisible(output)
}

if (sys.nframe() == 0L) summarize_full_refit_draws_main()
