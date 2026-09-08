#!/usr/bin/env Rscript
# Isolated three-fold-excluded nuisance fits for inner projection validation.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: fit_nested_target_states.R V19_LIBRARY SEED_BUNDLE OUTPUT")
  lib <- normalizePath(args[1], mustWork = TRUE)
  seed_root <- normalizePath(args[2], mustWork = TRUE)
  if (file.exists(args[3])) stop("output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/plan_nested_target_projection.R")
  code_files <- c("diagnosis/tate_common_weight/fit_nested_target_states.R",
    "diagnosis/tate_common_weight/plan_nested_target_projection.R",
    "diagnosis/tate_common_weight/NESTED_TARGET_FIT_PROTOCOL.md")
  code_hashes <- vapply(code_files, roce_sha256_file, "")
  stopifnot(.weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(seed_root, "sha256.txt")) ==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  manifest <- readLines(file.path(seed_root, "sha256.txt"))
  expected <- substr(manifest[substring(manifest, 67L) == "artifacts.rds"], 1L, 64L)
  stopifnot(length(expected) == 1L, roce_sha256_file(file.path(seed_root, "artifacts.rds")) == expected)
  bundle <- readRDS(file.path(seed_root, "artifacts.rds"))
  artifact <- bundle$group_result$artifacts[["0"]]
  target <- artifact$data_split$t
  ids <- lapply(artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info, `[[`, "target_idx")
  stopifnot(identical(sort(unlist(ids, use.names = FALSE)), seq_len(target$n)))
  plan <- .nested_target_projection_plan()
  subsets <- plan[!duplicated(plan$nuisance_training_key), ]
  states <- diagnostics <- list()
  for (row in seq_len(nrow(subsets))) {
    excluded <- as.integer(strsplit(subsets$excluded_folds[row], ":", fixed = TRUE)[[1]])
    training_ids <- unlist(ids[setdiff(1:5, excluded)], use.names = FALSE)
    excluded_ids <- unlist(ids[excluded], use.names = FALSE)
    stopifnot(!any(training_ids %in% excluded_ids))
    for (model_index in 1:3) {
      model <- c("propensity", "outcome_treated", "outcome_control")[model_index]
      model_ids <- if (model_index == 1L) training_ids else
        training_ids[target$A[training_ids] == as.integer(model_index == 2L)]
      x <- if (model_index == 1L) target$Z_site[model_ids, , drop = FALSE] else
        target$W_outcome[model_ids, , drop = FALSE]
      y <- if (model_index == 1L) target$A[model_ids] else target$Y[model_ids]
      stopifnot(all(y %in% 0:1), min(tabulate(y+1L, nbins = 2L)) >= 2L)
      seed <- as.integer(730000L + 100L*sum(excluded*c(100L, 10L, 1L)) + model_index)
      n_cv <- RoCE:::get_cv_fold_count(nrow(x), min_per_fold = 10L)
      cv_fold_id <- RoCE:::with_seed(seed, sample(rep(seq_len(n_cv), length.out = nrow(x))))
      min_training_class <- min(vapply(seq_len(n_cv), function(k)
        min(tabulate(y[cv_fold_id != k]+1L, nbins = 2L)), integer(1L)))
      warnings <- character()
      answer <- tryCatch(withCallingHandlers({
        if (min_training_class < 2L) stop("insufficient class support in a CV training split")
        glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1, nlambda = 100L,
          maxit = RoCE:::GLMNET_MAX_ITER, foldid = cv_fold_id, keep = TRUE, parallel = FALSE)
      }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = identity)
      failed <- inherits(answer, "error")
      key <- paste(subsets$nuisance_training_key[row], model, sep = ":")
      coefficient <- if (failed) NULL else as.numeric(stats::coef(answer, s = "lambda.min"))
      prediction_error <- if (failed) NA_real_ else max(abs(
        as.numeric(predict(answer, newx = x, s = "lambda.min", type = "response"))-
        plogis(drop(cbind(1, x) %*% coefficient))))
      eligible <- !failed && is.finite(prediction_error) && prediction_error < 1e-12 &&
        answer$glmnet.fit$jerr == 0L && all(is.finite(coefficient))
      states[[key]] <- list(fit = if (failed) NULL else answer, coefficients = coefficient,
        excluded_folds = excluded, training_ids = training_ids, model_ids = model_ids,
        seed = seed, cv_fold_id = cv_fold_id, warnings = warnings,
        failure = if (failed) conditionMessage(answer) else NULL, eligible = eligible)
      diagnostics[[key]] <- data.frame(key, model, seed, n_training = length(model_ids),
        n_cv, min_training_class, failed, eligible, warning_count = length(warnings),
        glmnet_jerr = if (failed) NA_integer_ else answer$glmnet.fit$jerr,
        prediction_error, failure = if (failed) conditionMessage(answer) else "")
    }
  }
  diagnostics <- do.call(rbind, diagnostics)
  stopifnot(length(states) == 30L, !anyDuplicated(diagnostics$seed),
    identical(code_hashes, vapply(code_files, roce_sha256_file, "")))
  roce_write_atomic_directory(args[3], function(stage) {
    saveRDS(states, file.path(stage, "nested_target_states.rds"))
    write.csv(diagnostics, file.path(stage, "fit_diagnostics.csv"), row.names = FALSE)
    write.csv(plan, file.path(stage, "nested_split_plan.csv"), row.names = FALSE)
    write.csv(data.frame(path = code_files, sha256 = code_hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c("inference_validated=FALSE", "resampling_draws=0", "new_mc_replications=0",
      "original_saved_fit_identity_claimed=FALSE", paste0("all_fits_eligible=", all(diagnostics$eligible)),
      paste0("source_manifest_sha256=", roce_sha256_file(file.path(seed_root, "sha256.txt")))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "nested target nuisance state capture")
  print(diagnostics, row.names = FALSE)
  if (!all(diagnostics$eligible)) stop("ineligible nuisance fits retained; do not proceed to selection")
}
if (sys.nframe() == 0L) main()
