#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(4L, 5L)) stop("usage: run_target_projection_validation.R V19_LIBRARY FIT_STATES SEED_BUNDLE OUTPUT [OUTER_FOLD]")
  outer_fold <- if (length(args) == 5L) suppressWarnings(as.numeric(args[5])) else 1L
  if (length(outer_fold) != 1L || is.na(outer_fold) || !outer_fold %in% 1:5) stop("OUTER_FOLD must be an integer from 1 to 5")
  outer_fold <- as.integer(outer_fold)
  validation_folds <- setdiff(1:5, outer_fold)
  lib <- normalizePath(args[1], mustWork = TRUE)
  states_root <- normalizePath(args[2], mustWork = TRUE)
  seed_root <- normalizePath(args[3], mustWork = TRUE); output <- args[4]
  if (file.exists(output)) stop("projection-validation output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/sparse_moment_projection.R")
  source("diagnosis/tate_common_weight/target_projection_validation.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_target_nuisance_system.R", equations)
  code_files <- c("diagnosis/tate_common_weight/run_target_projection_validation.R",
    "diagnosis/tate_common_weight/target_projection_validation.R",
    "diagnosis/tate_common_weight/sparse_moment_projection.R",
    "diagnosis/tate_common_weight/probe_target_nuisance_system.R",
    "diagnosis/tate_common_weight/TARGET_PROJECTION_VALIDATION_PROTOCOL.md")
  code_hashes <- vapply(code_files, roce_sha256_file, "")
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"),
    .weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(seed_root, "sha256.txt")) ==
      "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  checked_rds <- function(directory, name) {
    lines <- readLines(file.path(directory, "sha256.txt"))
    expected <- substr(lines[substring(lines, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, name)) == expected)
    readRDS(file.path(directory, name))
  }
  fits <- checked_rds(states_root, "target_cv_fit_states.rds")
  bundle <- checked_rds(seed_root, "artifacts.rds")
  artifact <- bundle$group_result$artifacts[["0"]]
  f <- artifact$direct_tate_results$one_round_crossfit
  ids <- lapply(f$intermediates$fold_info, `[[`, "target_idx")
  stopifnot(identical(sort(unlist(ids, use.names = FALSE)), seq_len(artifact$data_split$t$n)))
  limits <- list(ps = c(RoCE:::PROP_SCORE_LOWER, RoCE:::PROP_SCORE_UPPER),
                 outcome = c(RoCE:::OUTCOME_PRED_LOWER, RoCE:::OUTCOME_PRED_UPPER))
  used_keys <- character()
  model_theta <- function(inner = NULL) {
    suffix <- if (is.null(inner)) "NA" else as.character(inner)
    keys <- paste0(c("PS_1_", "OR_1_", "OR_0_"), outer_fold, "_", suffix)
    stopifnot(all(keys %in% names(fits)))
    used_keys <<- c(used_keys, keys)
    unlist(lapply(fits[keys], function(x) as.numeric(stats::coef(x, s = "lambda.min"))), use.names = FALSE)
  }
  block <- function(index) list(
    W = cbind(1, artifact$data_split$t$W_outcome[index, , drop = FALSE]),
    Z = cbind(1, artifact$data_split$t$Z_site[index, , drop = FALSE]),
    A = artifact$data_split$t$A[index], Y = artifact$data_split$t$Y[index], ids = index)
  candidates <- c(null = NA_real_, c0.5 = .5, c1 = 1, c2 = 2)
  losses <- attempts <- split_records <- list()
  for (inner in validation_folds) {
    train_ids <- unlist(ids[setdiff(validation_folds, inner)], use.names = FALSE)
    validation_ids <- ids[[inner]]
    stopifnot(!any(train_ids %in% c(ids[[outer_fold]], validation_ids)),
              !any(validation_ids %in% ids[[outer_fold]]))
    theta <- model_theta(inner)
    training <- block(train_ids); validation <- block(validation_ids)
    training_system <- equations$.target_tate_system(training, theta, limits, TRUE)
    validation_system <- equations$.target_tate_system(validation, theta, limits, TRUE)
    derivative_rows <- .target_score_derivative_rows(training, theta, limits, training_system)
    saved1 <- f$arm_results$mu1$fold_results[[outer_fold]]$target_only_inner[[paste0("k2_", inner)]]
    saved0 <- f$arm_results$mu0$fold_results[[outer_fold]]$target_only_inner[[paste0("k2_", inner)]]
    prediction_error <- max(abs(validation_system$predictions$p-saved1$prop_scores),
      abs(validation_system$predictions$p-saved0$prop_scores),
      abs(validation_system$predictions$m1-saved1$m_pred),
      abs(validation_system$predictions$m0-saved0$m_pred))
    stopifnot(prediction_error < 1e-12)
    split_records[[as.character(inner)]] <- list(training_ids = train_ids,
      validation_ids = validation_ids, nuisance_model_keys = tail(used_keys, 3L),
      prediction_error = prediction_error)
    for (candidate in names(candidates)) {
      answer <- if (candidate == "null") list(coefficients = numeric(length(theta)),
        penalties = rep(Inf, length(theta)), block_fits = list(), failure = NULL) else
        .fit_target_projection_candidate(training_system, derivative_rows, candidates[[candidate]])
      failed <- !is.null(answer$failure)
      risk <- if (failed) NA_real_ else -.5*sum(answer$coefficients*
        drop(validation_system$jacobian %*% answer$coefficients))+
        sum(validation_system$gradient*answer$coefficients)
      stopifnot(failed || is.finite(risk))
      losses[[length(losses)+1L]] <- data.frame(candidate, fold = inner,
        n_validation = length(validation_ids), loss = risk, failed,
        coefficient_l1 = sum(abs(answer$coefficients)),
        failure = if (failed) conditionMessage(answer$failure) else "")
      attempts[[paste(inner, candidate, sep = ":")]] <- answer
    }
  }
  losses <- do.call(rbind, losses)
  selection <- .select_target_projection_candidate(losses, validation_folds)
  # The outer evaluation data are not accessed until the selection is fixed.
  train_ids <- unlist(ids[validation_folds], use.names = FALSE)
  training <- block(train_ids); theta <- model_theta()
  training_system <- equations$.target_tate_system(training, theta, limits, TRUE)
  derivative_rows <- .target_score_derivative_rows(training, theta, limits, training_system)
  final <- if (selection$selected == "null") list(coefficients = numeric(length(theta)),
    penalties = rep(Inf, length(theta)), block_fits = list(), failure = NULL) else
    .fit_target_projection_candidate(training_system, derivative_rows, candidates[[selection$selected]])
  success <- is.null(final$failure)
  evaluation <- block(ids[[outer_fold]])
  heldout <- equations$.target_tate_system(evaluation, theta, limits)
  outer1 <- f$arm_results$mu1$fold_results[[outer_fold]]$target_only
  outer0 <- f$arm_results$mu0$fold_results[[outer_fold]]$target_only
  stopifnot(max(abs(heldout$predictions$p-outer1$prop_scores),
                abs(heldout$predictions$p-outer0$prop_scores),
                abs(heldout$predictions$m1-outer1$m_pred),
                abs(heldout$predictions$m0-outer0$m_pred)) < 1e-12,
            abs(heldout$estimate-outer1$estimate+outer0$estimate) < 1e-12)
  adjusted <- if (success) heldout$score-drop(heldout$moment_rows %*% final$coefficients) else
    rep(NA_real_, length(heldout$score))
  stopifnot(length(unique(used_keys)) == 15L, !any(train_ids %in% evaluation$ids))
  report <- data.frame(outer_fold = outer_fold, selected_candidate = selection$selected,
    selected_penalty_scale = unname(candidates[[selection$selected]]), selected_fit_succeeded = success,
    all_nonnull_ineligible = selection$all_nonnull_ineligible,
    coefficient_l1 = sum(abs(final$coefficients)), original_target_fold_tate = heldout$estimate,
    adjusted_target_fold_tate = mean(adjusted), original_score_sd = sd(heldout$score),
    adjusted_score_sd = sd(adjusted),
    failure = if (success) "" else conditionMessage(final$failure))
  stopifnot(identical(code_hashes, vapply(code_files, roce_sha256_file, "")))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(losses, file.path(stage, "validation_losses.csv"), row.names = FALSE)
    write.csv(selection$scores, file.path(stage, "candidate_risks.csv"), row.names = FALSE)
    write.csv(report, file.path(stage, "selected_holdout_report.csv"), row.names = FALSE)
    write.csv(data.frame(path = code_files, sha256 = code_hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    saveRDS(list(splits = split_records, attempts = attempts, selection = selection,
      final_fit = final, final_training_ids = train_ids, outer_evaluation_ids = evaluation$ids,
      original_score = heldout$score, adjusted_score = adjusted), file.path(stage, "selection_state.rds"))
    writeLines(c("target_projection_validation=complete", "sim_id=10013", paste0("outer_fold=", outer_fold),
      "procedure_status=exploratory_method_development", "candidate_grid=null;0.5;1;2",
      "selection_uses_only_inner_validation_risk=TRUE", "outer_records_excluded_from_used_models=TRUE",
      "outer_fold_previously_inspected_during_method_development=TRUE",
      "nuisance_refits=0", "resampling_draws=0", "new_mc_replications=0",
      "inference_validated=FALSE", paste0("selected_candidate=", selection$selected),
      paste0("selected_fit_succeeded=", success),
      paste0("source_states_manifest_sha256=", roce_sha256_file(file.path(states_root, "sha256.txt")))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "training-only target projection validation")
  print(selection$scores, row.names = FALSE); print(report, row.names = FALSE)
  if (!success) stop("selected projection fit failed; retained without switching candidates")
}

if (sys.nframe() == 0L) main()
