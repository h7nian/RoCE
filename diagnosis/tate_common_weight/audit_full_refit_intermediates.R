#!/usr/bin/env Rscript

# Read-only numerical audit for a saved full-refit resampling draw.  This
# replays only the row resampling; it never refits a nuisance model.

source("diagnosis/tate_common_weight/nuisance_cv_partition_audit.R")


roce_audit_full_refit_resampling <- function(
    saved_row_ids, saved_multiplicities, replay, fold_indices) {
  if (!identical(replay$original_row_ids, saved_row_ids) ||
      !identical(replay$multiplicities, saved_multiplicities)) {
    stop("replayed row maps or multiplicities do not match the draw.",
         call. = FALSE)
  }
  sites <- names(fold_indices)
  if (is.null(sites) || !identical(names(saved_row_ids), sites) ||
      !identical(names(saved_multiplicities), sites)) {
    stop("row-map, multiplicity, and fold site names must match exactly.",
         call. = FALSE)
  }
  rows <- list()
  outside_total <- 0L
  for (site in sites) {
    map <- saved_row_ids[[site]]
    multiplicity <- saved_multiplicities[[site]]
    if (!is.numeric(map) || !is.numeric(multiplicity) ||
        length(map) != length(multiplicity) || any(!is.finite(map)) ||
        any(!is.finite(multiplicity)) || any(map != floor(map)) ||
        any(multiplicity != floor(multiplicity)) || any(multiplicity < 0) ||
        !identical(as.integer(tabulate(as.integer(map), nbins = length(map))),
                   as.integer(multiplicity))) {
      stop("invalid row map or multiplicity counts at site ", site, ".",
           call. = FALSE)
    }
    for (k in seq_along(fold_indices[[site]])) {
      positions <- as.integer(fold_indices[[site]][[k]])
      origins <- as.integer(map[positions])
      outside <- sum(!origins %in% positions)
      outside_total <- outside_total + outside
      rows[[length(rows) + 1L]] <- data.frame(
        site = site,
        outer_fold = k,
        n_positions = length(positions),
        n_unique_origins = length(unique(origins)),
        n_duplicate_copies = length(origins) - length(unique(origins)),
        origins_outside_fold = outside,
        multiplicity_sum = sum(multiplicity[positions]),
        min_multiplicity = min(multiplicity[positions]),
        max_multiplicity = max(multiplicity[positions]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (outside_total != 0L) {
    stop("at least one resampled origin crossed its original outer fold.",
         call. = FALSE)
  }
  list(diagnostics = do.call(rbind, rows), cross_fold_origins = outside_total)
}

roce_audit_full_refit_identity <- function(saved, reference, fitted,
                                            tolerance = 1e-10) {
  if (!isTRUE(saved$summary$identity)) return(invisible(TRUE))
  fields <- c(
    "estimate_refit_relearned_weights", "estimate_refit_original_weights",
    "estimate_fixed_nuisance_original_weights",
    "estimate_fixed_nuisance_relearned_weights"
  )
  values <- suppressWarnings(as.numeric(saved$summary[1L, fields]))
  errors <- c(
    values - reference$estimate,
    fitted$estimate - reference$estimate,
    fitted$se - reference$se,
    fitted$variance - reference$variance,
    as.numeric(fitted$fold_weights - reference$fold_weights)
  )
  if (length(values) != length(fields) || any(!is.finite(errors)) ||
      max(abs(errors)) > tolerance) {
    stop("identity draw does not reproduce all four estimates, fitted SE, ",
         "variance, and fold weights.", call. = FALSE)
  }
  for (site in names(saved$result$original_row_ids)) {
    ids <- saved$result$original_row_ids[[site]]
    multiplicity <- saved$result$multiplicities[[site]]
    if (!identical(as.integer(ids), seq_along(ids)) ||
        any(as.integer(multiplicity) != 1L)) {
      stop("identity draw has a non-identity row map or multiplicity.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

roce_audit_nuisance_cv_partitions <- function(result, protocol) {
  if (!is.character(protocol) || length(protocol) != 1L || is.na(protocol) ||
      !protocol %in% c("row_level", "origin_grouped")) {
    stop("nuisance CV partition protocol must be row_level or origin_grouped.",
         call. = FALSE)
  }
  if (identical(protocol, "row_level")) {
    return(list(protocol = protocol, n_records = 0L, n_explicit = 0L))
  }
  if (!identical(result$nuisance_cv_grouping, "origin")) {
    stop("origin-grouped draw must declare nuisance_cv_grouping='origin'.",
         call. = FALSE)
  }
  records <- result$nuisance_cv_partitions
  if (!is.list(records) || length(records) == 0L) {
    stop("origin-grouped draw must retain nonempty nuisance CV partitions.",
         call. = FALSE)
  }
  explicit <- 0L
  for (index in seq_along(records)) {
    record <- records[[index]]
    required <- c("caller", "n_folds", "cv_group_id", "cv_fold_id", "valid")
    if (!is.list(record) || !all(required %in% names(record)) ||
        !isTRUE(record$valid)) {
      stop("invalid nuisance CV partition record ", index, ".",
           call. = FALSE)
    }
    if (!isTRUE(roce_cv_partition_record_valid(record))) {
      stop("nuisance CV partition record ", index,
           " fails independent group/fold validation.", call. = FALSE)
    }
    explicit <- explicit + as.integer(!is.null(record$cv_fold_id))
  }
  list(protocol = protocol, n_records = length(records), n_explicit = explicit)
}

audit_full_refit_intermediates_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop(
      paste0(
        "Usage: audit_full_refit_intermediates.R ",
        "DRAW_DIR INSTALLED_R_LIBRARY AUDIT_OUTPUT_DIR"
      ),
      call. = FALSE
    )
  }
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/full_refit_resampling.R")
  partition_validator_path <- normalizePath(
    "diagnosis/tate_common_weight/nuisance_cv_partition_audit.R",
    mustWork = TRUE
  )
  partition_validator_fingerprint <- roce_sha256_file(
    partition_validator_path
  )

  draw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
  project_library <- normalizePath(args[[2L]], mustWork = TRUE)
  output_dir <- args[[3L]]
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("audit output already exists; refusing to overwrite.", call. = FALSE)
  }
  output_parent <- normalizePath(
    dirname(output_dir), mustWork = dir.exists(dirname(output_dir))
  )
  protected <- c(draw_dir)
  is_within <- function(path, root) {
    identical(path, root) || startsWith(path, paste0(root, .Platform$file.sep))
  }
  if (any(vapply(protected, function(root) is_within(output_parent, root),
                 logical(1L)))) {
    stop("audit output must be outside the draw bundle.", call. = FALSE)
  }

  .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  if (!identical(
    normalizePath(find.package("RoCE"), mustWork = TRUE),
    normalizePath(file.path(project_library, "RoCE"), mustWork = TRUE)
  )) {
    stop("the explicitly selected RoCE installation was not loaded.",
         call. = FALSE)
  }

  read_checksum_manifest <- function(directory, expected_names, label) {
    manifest_path <- file.path(directory, "sha256.txt")
    lines <- readLines(manifest_path, warn = FALSE)
    hashes <- substr(lines, 1L, 64L)
    names <- substring(lines, 67L)
    if (length(lines) != length(expected_names) ||
        any(substr(lines, 65L, 66L) != "  ") ||
        any(!grepl("^[0-9a-f]{64}$", hashes)) || anyDuplicated(names) ||
        !setequal(names, expected_names)) {
      stop("malformed ", label, " checksum manifest.", call. = FALSE)
    }
    observed <- vapply(
      file.path(directory, names), roce_sha256_file, character(1L)
    )
    if (!identical(unname(observed), hashes)) {
      stop(label, " checksum mismatch.", call. = FALSE)
    }
    list(path = manifest_path, hashes = stats::setNames(hashes, names))
  }

  draw_manifest <- read_checksum_manifest(
    draw_dir, c("result.csv", "draw.rds"), "draw bundle"
  )
  saved <- readRDS(file.path(draw_dir, "draw.rds"))
  required <- c("result", "summary", "warnings", "source_bundle")
  if (!is.list(saved) || !all(required %in% names(saved)) ||
      !is.data.frame(saved$summary) || nrow(saved$summary) != 1L) {
    stop("draw.rds has an invalid top-level schema.", call. = FALSE)
  }
  csv_summary <- read.csv(
    file.path(draw_dir, "result.csv"), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  summaries_agree <- identical(names(csv_summary), names(saved$summary)) &&
    nrow(csv_summary) == 1L && all(vapply(names(csv_summary), function(field) {
      csv_value <- csv_summary[[field]]
      rds_value <- saved$summary[[field]]
      if (is.numeric(csv_value) && is.numeric(rds_value)) {
        isTRUE(all.equal(csv_value, rds_value, check.attributes = FALSE,
                         tolerance = 1e-14))
      } else {
        csv_text <- as.character(csv_value)
        rds_text <- as.character(rds_value)
        csv_text[is.na(csv_text)] <- ""
        rds_text[is.na(rds_text)] <- ""
        identical(csv_text, rds_text)
      }
    }, logical(1L)))
  if (!summaries_agree) {
    stop("result.csv and draw.rds summary disagree.", call. = FALSE)
  }
  if (!identical(as.character(saved$summary$status), "completed") ||
      nzchar(as.character(saved$summary$failure_message)) ||
      nzchar(as.character(saved$summary$failure_stage))) {
    stop("only a completed draw without a recorded failure can be audited.",
         call. = FALSE)
  }
  grouped_flag <- isTRUE(saved$summary$nuisance_cv_groups_duplicate_origins)
  leakage_flag <- isTRUE(
    saved$summary$nuisance_cv_duplicate_origin_leakage_possible
  )
  protocol <- if ("nuisance_cv_partition" %in% names(saved$summary)) {
    as.character(saved$summary$nuisance_cv_partition)
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

  source_bundle <- normalizePath(saved$source_bundle, mustWork = TRUE)
  if (is_within(output_parent, source_bundle)) {
    stop("audit output must be outside the source bundle.", call. = FALSE)
  }
  source_manifest <- read_checksum_manifest(
    source_bundle,
    c("results.csv", "diagnostic_qc.csv", "artifacts.rds", "metadata.txt"),
    "source bundle"
  )
  source_bundle_fingerprint <- roce_sha256_file(source_manifest$path)
  if (!identical(
    tolower(source_bundle_fingerprint),
    tolower(as.character(saved$summary$source_bundle_fingerprint))
  )) {
    stop("source-bundle fingerprint does not match the draw summary.",
         call. = FALSE)
  }
  source_metadata_lines <- readLines(
    file.path(source_bundle, "metadata.txt"), warn = FALSE
  )
  source_metadata <- stats::setNames(
    sub("^[^=]*=", "", source_metadata_lines),
    sub("=.*$", "", source_metadata_lines)
  )
  source_package_fingerprint <- unname(source_metadata["package_fingerprint"])
  if (length(source_package_fingerprint) != 1L ||
      is.na(source_package_fingerprint) ||
      !grepl("^[0-9a-f]{64}$", source_package_fingerprint)) {
    stop("source bundle lacks a valid package fingerprint.", call. = FALSE)
  }
  if (identical(protocol, "origin_grouped")) {
    if (!"source_package_fingerprint" %in% names(saved$summary) ||
        !identical(
          as.character(saved$summary$source_package_fingerprint),
          source_package_fingerprint
        )) {
      stop("origin-grouped draw source-package fingerprint is invalid.",
           call. = FALSE)
    }
  }
  installed_fingerprint <- .weight_calibration_installed_package_fingerprint(
    project_library, roce_sha256_file
  )
  if (!identical(
    installed_fingerprint,
    as.character(saved$summary$package_fingerprint)
  )) {
    stop("installed-package fingerprint does not match the draw.",
         call. = FALSE)
  }

  source_saved <- readRDS(file.path(source_bundle, "artifacts.rds"))
  rho <- as.numeric(saved$summary$rho)
  rho_key <- format(rho, scientific = FALSE, trim = TRUE)
  entry <- source_saved$group_result$artifacts[[rho_key]]
  if (is.null(entry) || is.null(entry$data_split) ||
      is.null(entry$direct_tate_results$one_round_crossfit)) {
    stop("selected source artifact is incomplete.", call. = FALSE)
  }
  reference <- entry$direct_tate_results$one_round_crossfit
  fitted <- saved$result$fitted
  reference_structure <- RoCE:::.validate_weight_bootstrap_result(reference)
  fitted_structure <- RoCE:::.validate_weight_bootstrap_result(fitted)
  if (!identical(reference_structure, fitted_structure)) {
    stop("reference and fitted site/fold structures differ.", call. = FALSE)
  }
  partition_audit <- roce_audit_nuisance_cv_partitions(
    saved$result, protocol
  )


  source_names <- names(reference$weights)
  sites <- c("t", source_names)
  reference_fold_info <- reference$intermediates$fold_info
  fold_indices <- c(
    list(t = lapply(reference_fold_info, `[[`, "target_idx")),
    stats::setNames(lapply(seq_along(source_names), function(j) {
      lapply(reference_fold_info, function(info) info$source_idx[[j]])
    }), source_names)
  )
  replay <- roce_refit_resample(
    entry$data_split, fold_indices,
    seed = as.integer(saved$summary$resample_seed),
    identity = isTRUE(saved$summary$identity)
  )
  resampling_audit <- roce_audit_full_refit_resampling(
    saved$result$original_row_ids, saved$result$multiplicities,
    replay, fold_indices
  )
  resampling_diagnostics <- resampling_audit$diagnostics
  duplicates_cross_fold <- resampling_audit$cross_fold_origins
  roce_audit_full_refit_identity(saved, reference, fitted)

  # Both treatment arms must expose exactly the same outer evaluation IDs.
  shared_arm_ids <- TRUE
  for (k in seq_len(fitted$n_folds)) {
    tate_info <- fitted$intermediates$fold_info[[k]]
    for (arm in fitted$arm_results[c("mu1", "mu0")]) {
      arm_info <- arm$intermediates$fold_info[[k]]
      shared_arm_ids <- shared_arm_ids &&
        identical(arm_info$target_idx, tate_info$target_idx) &&
        all(vapply(seq_along(source_names), function(j) {
          identical(arm_info$source_idx[[j]], tate_info$source_idx[[j]])
        }, logical(1L)))
    }
  }
  if (!shared_arm_ids) {
    stop("treated/control/TATE outer-fold observation IDs are not shared.",
         call. = FALSE)
  }

  tolerance <- 1e-10
  phi <- as.numeric(fitted$all_phi_tau)
  sizes <- c(
    fitted$intermediates$sample_sizes$n_t,
    fitted$intermediates$sample_sizes$n_source
  )
  if (length(phi) != sum(sizes) || any(!is.finite(phi))) {
    stop("saved all_phi_tau is incomplete or nonfinite.", call. = FALSE)
  }
  point_reconstructed <- mean(phi)
  starts <- cumsum(c(1L, head(as.integer(sizes), -1L)))
  within_site_ss <- vapply(seq_along(sizes), function(g) {
    values <- phi[starts[[g]]:(starts[[g]] + sizes[[g]] - 1L)]
    sum((values - mean(values))^2)
  }, numeric(1L))
  variance_reconstructed <- sum(within_site_ss) / length(phi)^2
  point_error <- abs(point_reconstructed - fitted$estimate)
  variance_error <- abs(variance_reconstructed - fitted$variance)
  if (point_error > tolerance || variance_error > tolerance) {
    stop("point or site-centered variance reconstruction failed.",
         call. = FALSE)
  }

  centered_moment <- function(x) mean((x - mean(x))^2)
  centered_covariance <- function(x, y) {
    mean((x - mean(x)) * (y - mean(y)))
  }
  moment_rows <- list()
  append_component_rows <- function(component, layer, outer_fold, inner_fold) {
    K <- length(source_names)
    V_ot_error <- abs(centered_moment(component$varphi_ot) - component$V_ot)
    for (j in seq_len(K)) {
      cross_error <- 0
      if (K > 1L) {
        for (l in seq_len(K)) {
          expected <- if (j == l) 0 else centered_covariance(
            component$zeta_components[[j]], component$zeta_components[[l]]
          )
          cross_error <- max(
            cross_error, abs(expected - component$C_cross[j, l])
          )
        }
      }
      moment_rows[[length(moment_rows) + 1L]] <<- data.frame(
        layer = layer,
        outer_fold = outer_fold,
        inner_fold = inner_fold,
        source = source_names[[j]],
        n_target = length(component$varphi_ot),
        n_source = length(component$xi_components[[j]]),
        V_ot_error = V_ot_error,
        V_t_error = abs(centered_moment(component$zeta_components[[j]]) -
                          component$V_t[[j]]),
        V_s_error = abs(centered_moment(component$xi_components[[j]]) -
                          component$V_s[[j]]),
        C_ot_error = abs(centered_covariance(
          component$varphi_ot, component$zeta_components[[j]]
        ) - component$C_ot[[j]]),
        C_cross_max_error = cross_error,
        stringsAsFactors = FALSE
      )
    }
  }
  for (k in seq_len(fitted$n_folds)) {
    append_component_rows(
      fitted$intermediates$fold_info[[k]], "outer", k, NA_integer_
    )
    components <- fitted$intermediates$inner_fold_info[[k]]$components
    for (component in components) {
      append_component_rows(
        component, "inner", k, as.integer(component$inner_fold)
      )
    }
  }
  moment_diagnostics <- do.call(rbind, moment_rows)
  max_moment_error <- max(unlist(moment_diagnostics[c(
    "V_ot_error", "V_t_error", "V_s_error", "C_ot_error",
    "C_cross_max_error"
  )]))
  if (!is.finite(max_moment_error) || max_moment_error > tolerance) {
    stop("raw influence moments do not reconstruct saved components.",
         call. = FALSE)
  }

  weight_rows <- do.call(rbind, lapply(seq_len(fitted$n_folds), function(k) {
    data.frame(
      outer_fold = k,
      source = source_names,
      weight = as.numeric(fitted$fold_weights[k, ]),
      included = as.logical(fitted$fold_source_included[k, ]),
      wald = as.numeric(fitted$fold_wald_statistics[k, ]),
      source_weight_sum = sum(fitted$fold_weights[k, ]),
      target_anchor_weight = 1 - sum(fitted$fold_weights[k, ]),
      stringsAsFactors = FALSE
    )
  }))
  if (any(!is.finite(weight_rows$weight)) || any(!is.finite(weight_rows$wald))) {
    stop("saved weights or Wald statistics are nonfinite.", call. = FALSE)
  }

  nuisance <- fitted$nuisance_fit_diagnostics
  if (!is.data.frame(nuisance) || nrow(nuisance) != 2L * fitted$n_folds ||
      !all(c("fold", "A_val") %in% names(nuisance))) {
    stop("fitted nuisance diagnostics have an invalid row schema.",
         call. = FALSE)
  }
  nuisance_numeric <- setdiff(names(nuisance), c("fold", "A_val"))
  nuisance_values <- suppressWarnings(as.numeric(unlist(
    nuisance[nuisance_numeric], use.names = FALSE
  )))
  if (length(nuisance_values) == 0L || any(!is.finite(nuisance_values))) {
    stop("fitted nuisance diagnostic values are missing or nonfinite.",
         call. = FALSE)
  }
  sum_columns <- function(pattern) {
    columns <- grep(pattern, names(nuisance), value = TRUE)
    if (length(columns) == 0L) return(NA_real_)
    sum(as.numeric(unlist(nuisance[columns], use.names = FALSE)))
  }
  max_columns <- function(pattern) {
    columns <- grep(pattern, names(nuisance), value = TRUE)
    if (length(columns) == 0L) return(NA_real_)
    max(as.numeric(unlist(nuisance[columns], use.names = FALSE)))
  }
  nuisance_summary <- data.frame(
    rows = nrow(nuisance),
    recorded_nonconverged = sum_columns("_nonconverged$"),
    recorded_line_search_failures = sum_columns("_line_search_failures$"),
    outcome_degenerate = sum_columns("outcome_degenerate$"),
    cv_invalid_fold_fits = sum_columns("_cv_invalid_fold_fits$"),
    cv_invalid_lambdas = sum_columns("_cv_invalid_lambdas$"),
    cv_path_tail_skipped_fold_fits =
      sum_columns("_cv_path_tail_skipped_fold_fits$"),
    support_floor_applied = sum_columns("_support_floor_applied$"),
    max_initial_dr_abs_coefficient =
      max_columns("^max_initial_dr_abs_coefficient$"),
    max_calibrated_dr_abs_coefficient =
      max_columns("^max_calibrated_dr_abs_coefficient$"),
    max_calibrated_outcome_abs_coefficient =
      max_columns("^max_calibrated_outcome_abs_coefficient$"),
    max_initial_dr_update_ratio =
      max_columns("^max_initial_dr_update_ratio$"),
    max_calibrated_dr_update_ratio =
      max_columns("^max_calibrated_dr_update_ratio$"),
    max_calibrated_outcome_update_ratio =
      max_columns("^max_calibrated_outcome_update_ratio$"),
    selected_loss_fields_available = any(grepl(
      "selected.*loss|loss.*selected", names(nuisance), ignore.case = TRUE
    )),
    interpretation = paste0(
      "nonconvergence and line-search counts are recorded solver status; ",
      "degeneracy, CV candidate exclusions, and support floors are retained ",
      "diagnostics and are not converted into invented validity cutoffs"
    ),
    stringsAsFactors = FALSE
  )

  gate <- data.frame(
    check = c(
      "draw_payload_hashes", "source_payload_hashes",
      "source_bundle_fingerprint", "installed_package_fingerprint",
      "row_map_replay", "duplicates_stay_within_outer_fold",
      "shared_arm_outer_ids", "point_reconstruction",
      "site_centered_variance_reconstruction", "raw_moment_reconstruction",
      "nuisance_diagnostics_schema", "nuisance_cv_partition_records"
    ),
    passed = TRUE,
    detail = c(
      paste(names(draw_manifest$hashes), collapse = ";"),
      paste(names(source_manifest$hashes), collapse = ";"),
      source_bundle_fingerprint,
      installed_fingerprint,
      paste0("identity=", saved$summary$identity),
      paste0("cross_fold_origins=", duplicates_cross_fold),
      paste0("shared=", shared_arm_ids),
      sprintf("abs_error=%.3g", point_error),
      sprintf("abs_error=%.3g", variance_error),
      sprintf("max_abs_error=%.3g", max_moment_error),
      sprintf(
        "rows=%d;nonconverged=%g;line_search_failures=%g;warnings=%d",
        nrow(nuisance), nuisance_summary$recorded_nonconverged,
        nuisance_summary$recorded_line_search_failures,
        length(saved$warnings)
      ),
      sprintf(
        "protocol=%s;records=%d;explicit=%d",
        protocol, partition_audit$n_records, partition_audit$n_explicit
      )
    ),
    stringsAsFactors = FALSE
  )

  summary <- data.frame(
    sim_id = saved$summary$sim_id,
    rho = rho,
    draw_id = saved$summary$draw_id,
    identity = saved$summary$identity,
    point_reconstructed = point_reconstructed,
    point_saved = fitted$estimate,
    point_abs_error = point_error,
    variance_reconstructed = variance_reconstructed,
    variance_saved = fitted$variance,
    variance_abs_error = variance_error,
    max_raw_moment_abs_error = max_moment_error,
    min_fold_source_weight = min(weight_rows$weight),
    max_fold_source_weight = max(weight_rows$weight),
    min_fold_source_weight_sum = min(weight_rows$source_weight_sum),
    max_fold_source_weight_sum = max(weight_rows$source_weight_sum),
    min_target_anchor_weight = min(weight_rows$target_anchor_weight),
    max_target_anchor_weight = max(weight_rows$target_anchor_weight),
    negative_weight_count = sum(weight_rows$weight < 0),
    weight_above_one_count = sum(weight_rows$weight > 1),
    warnings_recorded = length(saved$warnings),
    nuisance_cv_partition = protocol,
    nuisance_cv_partition_records = partition_audit$n_records,
    nuisance_cv_explicit_partition_records = partition_audit$n_explicit,
    source_package_fingerprint = source_package_fingerprint,
    refit_package_fingerprint = installed_fingerprint,
    gate = "passed",
    stringsAsFactors = FALSE
  )

  roce_write_atomic_directory(output_dir, function(staging) {
    write.csv(summary, file.path(staging, "summary.csv"), row.names = FALSE)
    write.csv(gate, file.path(staging, "gate.csv"), row.names = FALSE)
    write.csv(
      resampling_diagnostics,
      file.path(staging, "site_fold_resampling.csv"), row.names = FALSE
    )
    write.csv(
      moment_diagnostics,
      file.path(staging, "fold_source_moments.csv"), row.names = FALSE
    )
    write.csv(
      weight_rows,
      file.path(staging, "fold_source_weights.csv"), row.names = FALSE
    )
    write.csv(
      nuisance,
      file.path(staging, "nuisance_fit_diagnostics.csv"), row.names = FALSE
    )
    write.csv(
      nuisance_summary,
      file.path(staging, "nuisance_summary.csv"), row.names = FALSE
    )
    writeLines(
      c(
        "full_refit_intermediate_audit=passed",
        paste0("draw_bundle=", draw_dir),
        paste0("source_bundle=", source_bundle),
        paste0("installed_package_fingerprint=", installed_fingerprint),
        paste0("source_bundle_fingerprint=", source_bundle_fingerprint),
        paste0("source_package_fingerprint=", source_package_fingerprint),
        paste0("refit_package_fingerprint=", installed_fingerprint),
        paste0("nuisance_cv_partition=", protocol),
        paste0("nuisance_cv_partition_records=", partition_audit$n_records),
        paste0("partition_validator_fingerprint=",
               partition_validator_fingerprint),
        paste0("point_abs_error=", format(point_error, scientific = TRUE)),
        paste0("variance_abs_error=", format(variance_error, scientific = TRUE)),
        paste0("max_raw_moment_abs_error=",
               format(max_moment_error, scientific = TRUE))
      ),
      file.path(staging, "audit_passed.txt")
    )
    payloads <- c(
      "summary.csv", "gate.csv", "site_fold_resampling.csv",
      "fold_source_moments.csv", "fold_source_weights.csv",
      "nuisance_fit_diagnostics.csv", "nuisance_summary.csv",
      "audit_passed.txt"
    )
    hashes <- vapply(
      file.path(staging, payloads), roce_sha256_file, character(1L)
    )
    writeLines(
      paste(hashes, payloads, sep = "  "),
      file.path(staging, "sha256.txt")
    )
  }, caller = "full-refit intermediate audit")
  message("[passed] full-refit intermediate audit: ", output_dir)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  audit_full_refit_intermediates_main()
}
