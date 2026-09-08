#!/usr/bin/env Rscript
# Inner-score development check. No common weights or confidence intervals.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) stop("usage: run_nested_target_projection.R LIBRARY SEED_BUNDLE ORIGINAL_STATES NESTED_STATES OUTPUT")
  lib <- normalizePath(args[1], mustWork = TRUE)
  if (file.exists(args[5])) stop("output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/sparse_moment_projection.R")
  source("diagnosis/tate_common_weight/target_projection_validation.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_target_nuisance_system.R", equations)
  code_files <- paste0("diagnosis/tate_common_weight/", c("run_nested_target_projection.R",
    "sparse_moment_projection.R", "target_projection_validation.R", "probe_target_nuisance_system.R",
    "NESTED_TARGET_PROJECTION_PROTOCOL.md"))
  code_hashes <- vapply(code_files, roce_sha256_file, "")
  stopifnot(.weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(args[2], "sha256.txt")) ==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  checked_rds <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  artifact <- checked_rds(args[2], "artifacts.rds")$group_result$artifacts[["0"]]
  original <- checked_rds(args[3], "target_cv_fit_states.rds")
  nested <- checked_rds(args[4], "nested_target_states.rds")
  stopifnot(length(nested) == 30L, all(vapply(nested, function(x) isTRUE(x$eligible), logical(1L))))
  fit <- artifact$direct_tate_results$one_round_crossfit
  ids <- lapply(fit$intermediates$fold_info, `[[`, "target_idx")
  target <- artifact$data_split$t
  stopifnot(identical(sort(unlist(ids, use.names = FALSE)), seq_len(target$n)))
  block <- function(index) list(W = cbind(1, target$W_outcome[index, , drop = FALSE]),
    Z = cbind(1, target$Z_site[index, , drop = FALSE]), A = target$A[index], Y = target$Y[index])
  limits <- list(ps = c(RoCE:::PROP_SCORE_LOWER, RoCE:::PROP_SCORE_UPPER),
    outcome = c(RoCE:::OUTCOME_PRED_LOWER, RoCE:::OUTCOME_PRED_UPPER))
  candidates <- c(null = NA_real_, c0.5 = .5, c1 = 1, c2 = 2)
  candidate_fit <- function(system, rows, candidate) {
    if (candidate == "null") return(list(coefficients = numeric(length(system$gradient)),
      penalties = rep(Inf, length(system$gradient)), block_fits = list(), failure = NULL))
    .fit_target_projection_candidate(system, rows, candidates[[candidate]])
  }
  records <- reports <- list()
  for (outer_fold in 1:5) for (evaluation_fold in setdiff(1:5, outer_fold)) {
    validation_folds <- setdiff(1:5, c(outer_fold, evaluation_fold))
    losses <- attempts <- splits <- list()
    for (validation_fold in validation_folds) {
      excluded <- sort(c(outer_fold, evaluation_fold, validation_fold))
      train_ids <- unlist(ids[setdiff(1:5, excluded)], use.names = FALSE)
      keys <- paste0("exclude_", paste(excluded, collapse = "_"), ":",
        c("propensity", "outcome_treated", "outcome_control"))
      stopifnot(all(keys %in% names(nested)), all(vapply(nested[keys], function(x)
        identical(x$training_ids, train_ids) && identical(x$excluded_folds, excluded), logical(1L))))
      theta <- unlist(lapply(nested[keys], `[[`, "coefficients"), use.names = FALSE)
      training <- block(train_ids)
      training_system <- equations$.target_tate_system(training, theta, limits, TRUE)
      validation_system <- equations$.target_tate_system(block(ids[[validation_fold]]), theta, limits, TRUE)
      derivative_rows <- .target_score_derivative_rows(training, theta, limits, training_system)
      splits[[as.character(validation_fold)]] <- list(training_ids = train_ids,
        validation_ids = ids[[validation_fold]], nuisance_keys = keys,
        training_clipping = training_system$clipping, validation_clipping = validation_system$clipping)
      for (candidate in names(candidates)) {
        answer <- candidate_fit(training_system, derivative_rows, candidate)
        failed <- !is.null(answer$failure)
        risk <- if (failed) NA_real_ else -.5*sum(answer$coefficients*
          drop(validation_system$jacobian %*% answer$coefficients))+
          sum(validation_system$gradient*answer$coefficients)
        losses[[length(losses)+1L]] <- data.frame(candidate, fold = validation_fold,
          n_validation = length(ids[[validation_fold]]), loss = risk, failed)
        attempts[[paste(validation_fold, candidate, sep = ":")]] <- answer
      }
    }
    losses <- do.call(rbind, losses)
    selection <- .select_target_projection_candidate(losses, validation_folds)
    train_ids <- unlist(ids[validation_folds], use.names = FALSE)
    keys <- paste0(c("PS_1_", "OR_1_", "OR_0_"), outer_fold, "_", evaluation_fold)
    stopifnot(all(keys %in% names(original)))
    theta <- unlist(lapply(original[keys], function(x) as.numeric(stats::coef(x, s = "lambda.min"))), use.names = FALSE)
    training <- block(train_ids)
    system <- equations$.target_tate_system(training, theta, limits, TRUE)
    final <- candidate_fit(system, .target_score_derivative_rows(training, theta, limits, system), selection$selected)
    heldout <- equations$.target_tate_system(block(ids[[evaluation_fold]]), theta, limits)
    saved1 <- fit$arm_results$mu1$fold_results[[outer_fold]]$target_only_inner[[paste0("k2_", evaluation_fold)]]
    saved0 <- fit$arm_results$mu0$fold_results[[outer_fold]]$target_only_inner[[paste0("k2_", evaluation_fold)]]
    prediction_error <- max(abs(heldout$predictions$p-saved1$prop_scores),
      abs(heldout$predictions$p-saved0$prop_scores), abs(heldout$predictions$m1-saved1$m_pred),
      abs(heldout$predictions$m0-saved0$m_pred))
    stopifnot(prediction_error < 1e-12, abs(heldout$estimate-saved1$estimate+saved0$estimate) < 1e-12,
      !any(train_ids %in% c(ids[[outer_fold]], ids[[evaluation_fold]])))
    success <- is.null(final$failure)
    adjusted <- if (success) heldout$score-drop(heldout$moment_rows %*% final$coefficients) else
      rep(NA_real_, length(heldout$score))
    key <- paste(outer_fold, evaluation_fold, sep = ":")
    records[[key]] <- list(splits = splits, losses = losses, attempts = attempts, selection = selection,
      final_fit = final, final_training_ids = train_ids, evaluation_ids = ids[[evaluation_fold]],
      final_nuisance_keys = keys, original_score = heldout$score, adjusted_score = adjusted,
      training_clipping = system$clipping, evaluation_clipping = heldout$clipping)
    reports[[key]] <- data.frame(outer_fold, evaluation_fold, selected_candidate = selection$selected,
      selected_penalty_scale = unname(candidates[[selection$selected]]), success,
      all_nonnull_ineligible = selection$all_nonnull_ineligible,
      coefficient_l1 = sum(abs(final$coefficients)), prediction_error,
      original_estimate = mean(heldout$score), adjusted_estimate = mean(adjusted),
      original_score_sd = sd(heldout$score), adjusted_score_sd = sd(adjusted))
  }
  reports <- do.call(rbind, reports)
  stopifnot(nrow(reports) == 20L, identical(code_hashes, vapply(code_files, roce_sha256_file, "")))
  roce_write_atomic_directory(args[5], function(stage) {
    saveRDS(records, file.path(stage, "inner_projection_states.rds"))
    write.csv(reports, file.path(stage, "inner_projection_reports.csv"), row.names = FALSE)
    write.csv(data.frame(path = code_files, sha256 = code_hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c("inference_validated=FALSE", "common_weights_recomputed=FALSE", "nuisance_refits=0",
      "resampling_draws=0", "new_mc_replications=0",
      paste0("input_", 2:4, "_manifest_sha256=", vapply(args[2:4], function(root)
        roce_sha256_file(file.path(root, "sha256.txt")), ""))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "nested target projection validation")
  print(reports, row.names = FALSE)
  if (!all(reports$success)) stop("selected projection failures retained; do not aggregate successful pairs only")
}
if (sys.nframe() == 0L) main()
