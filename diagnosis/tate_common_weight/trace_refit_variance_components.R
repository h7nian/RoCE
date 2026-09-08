#!/usr/bin/env Rscript
# Four plug-in variance reconstructions on the SAME resampled observations.
# These are within-draw diagnostics, not variances across bootstrap draws.

trace_refit_variance_components_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("Usage: trace_refit_variance_components.R DRAW_DIR ROCE_LIBRARY OUTPUT_DIR")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  draw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
  project_library <- normalizePath(args[[2L]], mustWork = TRUE)
  output_dir <- args[[3L]]
  if (file.exists(output_dir)) stop("Refusing to overwrite existing output.")
  verify <- function(directory, expected_files) {
    lines <- readLines(file.path(directory, "sha256.txt"))
    files <- substring(lines, 67L)
    hashes <- substr(lines, 1L, 64L)
    if (length(lines) != length(expected_files) || anyDuplicated(files) ||
        !setequal(files, expected_files) ||
        any(substr(lines, 65L, 66L) != "  ") ||
        any(!grepl("^[0-9a-f]{64}$", hashes)) ||
        !identical(unname(vapply(file.path(directory, files), roce_sha256_file,
                                 character(1L))), hashes)) stop("Input checksum mismatch.")
    roce_sha256_file(file.path(directory, "sha256.txt"))
  }
  draw_hash <- verify(draw_dir, c("result.csv", "draw.rds"))
  saved <- readRDS(file.path(draw_dir, "draw.rds"))
  if (!identical(saved$summary$status, "completed")) stop("A completed draw is required.")
  source_dir <- normalizePath(saved$source_bundle, mustWork = TRUE)
  for (protected in c(draw_dir, source_dir)) {
    parent <- dirname(output_dir)
    destination <- file.path(normalizePath(parent, mustWork = TRUE), basename(output_dir))
    if (identical(destination, protected) || startsWith(destination, paste0(protected, "/"))) {
      stop("Output must be outside scientific input bundles.")
    }
  }
  source_hash <- verify(source_dir, c("results.csv", "diagnostic_qc.csv", "artifacts.rds", "metadata.txt"))
  if (!identical(source_hash, saved$summary$source_bundle_fingerprint)) stop("Source provenance mismatch.")
  .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  package_hash <- .weight_calibration_installed_package_fingerprint(project_library, roce_sha256_file)
  if (!identical(normalizePath(find.package("RoCE")), normalizePath(file.path(project_library, "RoCE"))) ||
      !identical(package_hash, saved$summary$package_fingerprint)) stop("Installed package mismatch.")
  rho_key <- format(saved$summary$rho, scientific = FALSE, trim = TRUE)
  entry <- readRDS(file.path(source_dir, "artifacts.rds"))$group_result$artifacts[[rho_key]]
  reference <- entry$direct_tate_results$one_round_crossfit
  fitted <- saved$result$fitted
  RoCE:::.validate_weight_bootstrap_result(reference)
  RoCE:::.validate_weight_bootstrap_result(fitted)
  sizes <- reference$intermediates$sample_sizes
  sites <- c("t", names(reference$weights))
  group_sizes <- c(sizes$n_t, sizes$n_source)
  offsets <- c(0L, head(cumsum(group_sizes), -1L))
  maps <- saved$result$original_row_ids
  if (!identical(names(maps), sites)) stop("Row-map site order mismatch.")
  for (j in seq_along(sites)) {
    idx <- maps[[j]]
    if (!is.numeric(idx) || length(idx) != group_sizes[j] ||
        any(!is.finite(idx)) || any(idx != floor(idx)) ||
        any(idx < 1) || any(idx > group_sizes[j]) ||
        !identical(tabulate(idx, nbins = group_sizes[j]), saved$result$multiplicities[[j]])) {
      stop("Invalid row map or multiplicities.")
    }
  }
  original_positions <- unlist(Map(function(idx, offset) idx + offset, maps, offsets), use.names = FALSE)
  target_rows <- list()
  for (mode in c("fixed_nuisance", "refitted_nuisance")) {
    object <- if (mode == "fixed_nuisance") reference else fitted
    for (arm in c("mu1", "mu0")) {
      A_val <- if (arm == "mu1") 1L else 0L
      propensity <- prediction <- saved_phi <- rep(NA_real_, sizes$n_t)
      fold_id <- integer(sizes$n_t)
      for (k in seq_len(object$n_folds)) {
        idx <- object$intermediates$fold_info[[k]]$target_idx
        values <- object$arm_results[[arm]]$fold_results[[k]]$target_only
        if (any(lengths(values[c("prop_scores", "m_pred", "varphi_ot")]) != length(idx))) {
          stop("Target nuisance prediction length mismatch.")
        }
        propensity[idx] <- values$prop_scores
        prediction[idx] <- values$m_pred
        saved_phi[idx] <- values$varphi_ot + values$estimate
        fold_id[idx] <- k
      }
      if (mode == "fixed_nuisance") {
        propensity <- propensity[maps$t]
        prediction <- prediction[maps$t]
        saved_phi <- saved_phi[maps$t]
        fold_id <- fold_id[maps$t]
      }
      A <- entry$data_split$t$A[maps$t]
      Y <- entry$data_split$t$Y[maps$t]
      if (any(!is.finite(c(propensity, prediction, saved_phi))) ||
          any(propensity <= 0 | propensity >= 1) ||
          any(prediction < 0 | prediction > 1) || any(fold_id == 0L)) {
        stop("Invalid binary-outcome target nuisance predictions.")
      }
      arm_probability <- pmax(if (A_val == 1L) propensity else 1 - propensity,
                              RoCE:::PROP_SCORE_LOWER)
      raw <- prediction + as.numeric(A == A_val) * (Y - prediction) / arm_probability
      if (max(abs(raw - saved_phi)) > 1e-10) stop("Target AIPW formula does not reproduce saved values.")
      target_rows[[length(target_rows) + 1L]] <- data.frame(
        mode = mode, arm = arm, resampled_position = seq_len(sizes$n_t),
        original_row_id = maps$t, fold = fold_id, A = A, Y = Y,
        propensity_A1 = propensity, arm_probability = arm_probability,
        outcome_prediction = prediction, aipw_raw = raw,
        at_propensity_lower_bound = propensity == RoCE:::PROP_SCORE_LOWER,
        at_propensity_upper_bound = propensity == RoCE:::PROP_SCORE_UPPER
      )
    }
  }
  target_predictions <- do.call(rbind, target_rows)
  evaluate <- function(object, weights, arm = NULL) {
    info <- if (is.null(arm)) object$intermediates$fold_info else
      object$arm_results[[arm]]$intermediates$fold_info
    RoCE:::.compute_phase3_all_phi(
      object$n_folds, weights, info, sizes$n_t, sizes$n_source, sizes$N_all,
      names(reference$weights), length(reference$weights), verbose = FALSE
    )
  }
  specifications <- list(
    A_fixed_nuisance_original_weights = list(reference, reference$fold_weights, TRUE),
    B_fixed_nuisance_relearned_weights = list(reference, saved$result$matched_fixed_nuisance$fold_weights, TRUE),
    C_refit_nuisance_original_weights = list(fitted, reference$fold_weights, FALSE),
    D_refit_nuisance_relearned_weights = list(fitted, fitted$fold_weights, FALSE)
  )
  expected <- c(saved$result$matched_fixed_nuisance$estimate_original_weights,
                saved$result$matched_fixed_nuisance$estimate_relearned_weights,
                saved$result$estimate_refit_original_weights,
                saved$result$estimate_refit_relearned_weights)
  rows <- list()
  totals <- list()
  tails <- list()
  for (stage in names(specifications)) {
    spec <- specifications[[stage]]
    phi <- evaluate(spec[[1L]], spec[[2L]])
    phi1 <- evaluate(spec[[1L]], spec[[2L]], "mu1")
    phi0 <- evaluate(spec[[1L]], spec[[2L]], "mu0")
    if (spec[[3L]]) {
      phi <- phi[original_positions]
      phi1 <- phi1[original_positions]
      phi0 <- phi0[original_positions]
    }
    contrast_error <- max(abs(phi - (phi1 - phi0)))
    if (!is.finite(contrast_error) || contrast_error > 1e-10) stop("Arm contrast reconstruction failed.")
    stage_variance <- 0
    for (j in seq_along(sites)) {
      idx <- seq_len(group_sizes[j]) + offsets[j]
      centered <- phi[idx] - mean(phi[idx])
      centered1 <- phi1[idx] - mean(phi1[idx])
      centered0 <- phi0[idx] - mean(phi0[idx])
      v1 <- sum(centered1^2) / sizes$N_all^2
      v0 <- sum(centered0^2) / sizes$N_all^2
      covariance <- sum(centered1 * centered0) / sizes$N_all^2
      variance <- sum(centered^2) / sizes$N_all^2
      error <- abs(variance - (v1 + v0 - 2 * covariance))
      if (!is.finite(error) || error > 1e-12) stop("Cross-arm covariance identity failed.")
      rank_index <- order(abs(centered), decreasing = TRUE)
      squared_total <- sum(centered^2)
      rows[[length(rows) + 1L]] <- data.frame(
        stage = stage, site = sites[j], n = group_sizes[j],
        variance_contribution = variance, variance_mu1 = v1, variance_mu0 = v0,
        covariance_mu1_mu0 = covariance, contrast_variance_error = error,
        max_abs_centered_scaled_phi = max(abs(centered)),
        largest_observation_variance_share = if (squared_total > 0) centered[rank_index[1L]]^2 / squared_total else NA_real_,
        largest_five_variance_share = if (squared_total > 0) sum(centered[head(rank_index, 5L)]^2) / squared_total else NA_real_
      )
      selected <- head(rank_index, 5L)
      tails[[length(tails) + 1L]] <- data.frame(
        stage = stage, site = sites[j], resampled_position = selected,
        original_row_id = maps[[j]][selected],
        A = entry$data_split[[sites[j]]]$A[maps[[j]][selected]],
        Y = entry$data_split[[sites[j]]]$Y[maps[[j]][selected]],
        centered_scaled_phi = centered[selected]
      )
      stage_variance <- stage_variance + variance
    }
    expected_estimate <- expected[match(stage, names(specifications))]
    if (abs(mean(phi) - expected_estimate) > 1e-12) stop("Saved point estimate disagrees with reconstruction.")
    totals[[stage]] <- data.frame(stage = stage, estimate = mean(phi),
      plugin_variance = stage_variance, plugin_se = sqrt(stage_variance),
      arm_contrast_max_error = contrast_error)
  }
  totals <- do.call(rbind, totals)
  rows <- do.call(rbind, rows)
  tails <- do.call(rbind, tails)
  if (abs(tail(totals$plugin_variance, 1L) - fitted$variance) > 1e-12) stop("Final variance mismatch.")
  script_hash <- roce_sha256_file("diagnosis/tate_common_weight/trace_refit_variance_components.R")
  if (!identical(draw_hash, verify(draw_dir, c("result.csv", "draw.rds"))) ||
      !identical(source_hash, verify(source_dir, c("results.csv", "diagnostic_qc.csv", "artifacts.rds", "metadata.txt")))) {
    stop("Inputs changed during trace.")
  }
  roce_write_atomic_directory(output_dir, function(staging) {
    write.csv(totals, file.path(staging, "plugin_variance_stages.csv"), row.names = FALSE)
    write.csv(rows, file.path(staging, "site_arm_variance.csv"), row.names = FALSE)
    write.csv(tails, file.path(staging, "largest_contributions.csv"), row.names = FALSE)
    write.csv(target_predictions, file.path(staging, "target_nuisance_predictions.csv"), row.names = FALSE)
    writeLines(c(paste0("trace_script_fingerprint=", script_hash),
      paste0("draw_checksum_manifest_fingerprint=", draw_hash),
      paste0("source_checksum_manifest_fingerprint=", source_hash),
      paste0("package_fingerprint=", package_hash), "reconstruction_gate=passed",
      "interpretation=within-draw plugin variances, not variances across resampling draws",
      "inference_validated=FALSE"), file.path(staging, "metadata.txt"))
    files <- c("plugin_variance_stages.csv", "site_arm_variance.csv", "largest_contributions.csv", "target_nuisance_predictions.csv", "metadata.txt")
    hashes <- vapply(file.path(staging, files), roce_sha256_file, character(1L))
    writeLines(paste(hashes, files, sep = "  "), file.path(staging, "sha256.txt"))
  })
  print(totals)
  print(rows[c("stage", "site", "variance_contribution", "variance_mu1", "variance_mu0", "covariance_mu1_mu0")])
  invisible(totals)
}

if (sys.nframe() == 0L) trace_refit_variance_components_main()
