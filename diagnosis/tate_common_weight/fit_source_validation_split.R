#!/usr/bin/env Rscript
# Source validation splits: global outer 1, validation 2:5, source s1.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% 3:5) stop("usage: fit_source_validation_split.R LIBRARY SEED_BUNDLE OUTPUT [VALIDATION_FOLD] [CONFIG]")
  working_config <- if (length(args) == 5L) args[5] else "C1"
  validation_fold <- if (length(args) >= 4L) suppressWarnings(as.numeric(args[4])) else 2L
  if (length(validation_fold) != 1L || is.na(validation_fold) || !validation_fold %in% 2:5) stop("VALIDATION_FOLD must be in 2:5")
  validation_fold <- as.integer(validation_fold)
  fold_map <- c(validation_fold, setdiff(2:5, validation_fold))
  lib <- normalizePath(args[1], mustWork = TRUE)
  if (file.exists(args[3])) stop("output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  code <- c("diagnosis/tate_common_weight/fit_source_validation_split.R",
    "diagnosis/tate_common_weight/SOURCE_VALIDATION_FIT_PROTOCOL.md",
    "diagnosis/tate_common_weight/source_working_basis.R")
  code_hashes <- vapply(code, roce_sha256_file, "")
  stopifnot(.weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(args[2], "sha256.txt")) ==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  manifest <- readLines(file.path(args[2], "sha256.txt"))
  expected <- substr(manifest[substring(manifest, 67L) == "artifacts.rds"], 1L, 64L)
  stopifnot(length(expected) == 1L, roce_sha256_file(file.path(args[2], "artifacts.rds")) == expected)
  artifact <- readRDS(file.path(args[2], "artifacts.rds"))$group_result$artifacts[["0"]]
  data <- .source_working_basis(artifact$data_split, working_config)
  info <- artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  global_ids <- list(target = lapply(info, `[[`, "target_idx"),
    source = lapply(info, function(x) x$source_idx[[1L]]))
  make_views <- function(site_data, ids) {
    views <- lapply(ids[fold_map], function(index) list(original_idx = index, n = length(index)))
    attr(views, ".data_ref") <- site_data
    views
  }
  target_folds <- make_views(data$t, global_ids$target)
  source_folds <- list(s1 = make_views(data$s1, global_ids$source))
  splits <- lapply(2:4, function(k) {
    training <- setdiff(2:4, k)
    list(local_calibration_fold = k, global_calibration_fold = fold_map[k],
      target_training_ids = RoCE:::combine_folds(target_folds, training)$original_idx,
      source_training_ids = RoCE:::combine_folds(source_folds$s1, training)$original_idx,
      target_calibration_ids = target_folds[[k]]$original_idx,
      source_calibration_ids = source_folds$s1[[k]]$original_idx)
  })
  for (split in splits) for (site in names(global_ids)) {
    forbidden <- unlist(global_ids[[site]][c(1L, validation_fold)], use.names = FALSE)
    training <- split[[paste0(site, "_training_ids")]]
    calibration <- split[[paste0(site, "_calibration_ids")]]
    stopifnot(!any(c(training, calibration) %in% forbidden), !any(training %in% calibration))
  }
  results <- list()
  for (arm in 0:1) {
    warnings <- character(); initial <- inputs <- list()
    answer <- tryCatch(RoCE:::with_seed(880102L+10L*(validation_fold-2L)+arm, withCallingHandlers({
      cached_lambda <- NULL
      for (k in 2:4) {
        training <- RoCE:::combine_folds(target_folds, setdiff(2:4, k))
        calibration <- RoCE:::materialize_fold(target_folds, k)
        alpha <- RoCE:::fit_initial_outcome(training$W_outcome, training$Y, training$A, arm,
          lambda = cached_lambda, nlambda = 100L, family = "binomial", lambda_rule = "min",
          cv_group_id = training$cv_group_id)
        selected_lambda <- attr(alpha, "lambda_used")
        if (is.null(cached_lambda) && length(selected_lambda) == 1L && is.finite(selected_lambda))
          cached_lambda <- selected_lambda
        key <- paste0("k2_", k)
        initial[[key]] <- alpha
        inputs[[key]] <- list(alpha_init = alpha, mean_phi = c(1, colMeans(training$Z_site)),
          mean_grad_psi_init = RoCE:::.mean_glm_gradient_site_basis(
            calibration$W_outcome, calibration$Z_site, alpha, 1L, 1L))
      }
      RoCE:::process_source_site(s = "s1", source_folds = source_folds, target_folds = target_folds,
        k1 = 1L, n_folds = 4L, A_val = arm, M_tau = 5, M_tau_inference = 5,
        data_split = data, get_fold_inputs = function(site, k) inputs[[paste0("k2_", k)]],
        family_int = 1L, link_int = 1L, use_lambda_cache = TRUE,
        nuisance_nlambda = 100L, nuisance_lambda_rule = "min")
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })), error = identity)
    failed <- inherits(answer, "error")
    results[[paste0("mu", arm)]] <- list(failed = failed, warnings = warnings,
      failure = if (failed) conditionMessage(answer) else NULL,
      result = if (failed) NULL else answer, initial_outcomes = initial,
      target_summaries = inputs, seed = 880102L+10L*(validation_fold-2L)+arm)
  }
  stopifnot(identical(code_hashes, vapply(code, roce_sha256_file, "")))
  roce_write_atomic_directory(args[3], function(stage) {
    saveRDS(list(results = results, splits = splits, global_fold_map = fold_map, working_config = working_config,
      validation_ids = lapply(global_ids, `[[`, validation_fold), excluded_outer_ids = lapply(global_ids, `[[`, 1L)),
      file.path(stage, "source_validation_states.rds"))
    write.csv(data.frame(arm = 0:1, failed = vapply(results, `[[`, logical(1L), "failed"),
      warning_count = vapply(results, function(x) length(x$warnings), integer(1L))),
      file.path(stage, "fit_status.csv"), row.names = FALSE)
    write.csv(data.frame(path = code, sha256 = code_hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c("inference_validated=FALSE", "projection_coefficients_fitted=FALSE",
      "original_saved_fit_identity_claimed=FALSE", "resampling_draws=0",
      paste0("seed_manifest_sha256=", roce_sha256_file(file.path(args[2], "sha256.txt")))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "source validation nuisance fit")
  if (any(vapply(results, `[[`, logical(1L), "failed"))) stop("source fit failures retained for review")
  cat("Both source arms fitted; numerical/KKT audit still required.\n")
}
if (sys.nframe() == 0L) main()
