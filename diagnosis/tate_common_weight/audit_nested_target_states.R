#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: audit_nested_target_states.R STATES_ROOT SEED_BUNDLE OUTPUT")
  suppressPackageStartupMessages(library(glmnet))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  checked_rds <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  states <- checked_rds(args[1], "nested_target_states.rds")
  artifact <- checked_rds(args[2], "artifacts.rds")$group_result$artifacts[["0"]]
  ids <- lapply(artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info, `[[`, "target_idx")
  target <- artifact$data_split$t
  stopifnot(length(states) == 30L, !anyDuplicated(names(states)),
    identical(sort(unlist(ids, use.names = FALSE)), seq_len(target$n)))
  rows <- list()
  triples <- combn(1:5, 3, simplify = FALSE)
  expected_keys <- unlist(lapply(triples, function(triple) paste0(
    "exclude_", paste(triple, collapse = "_"), ":",
    c("propensity", "outcome_treated", "outcome_control"))), use.names = FALSE)
  stopifnot(setequal(names(states), expected_keys))
  for (triple in triples) {
    expected_training <- unlist(ids[setdiff(1:5, triple)], use.names = FALSE)
    for (model_index in 1:3) {
      model <- c("propensity", "outcome_treated", "outcome_control")[model_index]
      key <- paste0("exclude_", paste(triple, collapse = "_"), ":", model)
      state <- states[[key]]
      expected_model_ids <- if (model_index == 1L) expected_training else
        expected_training[target$A[expected_training] == as.integer(model_index == 2L)]
      expected_seed <- as.integer(730000L+100L*sum(triple*c(100L, 10L, 1L))+model_index)
      stopifnot(identical(state$excluded_folds, triple),
        identical(state$training_ids, expected_training), identical(state$model_ids, expected_model_ids),
        identical(state$seed, expected_seed), is.null(state$failure), isTRUE(state$eligible),
        !is.null(state$fit), state$fit$glmnet.fit$jerr == 0L)
      y <- if (model_index == 1L) target$A[expected_model_ids] else target$Y[expected_model_ids]
      fold_id <- state$cv_fold_id
      n_cv <- max(fold_id)
      stopifnot(length(fold_id) == length(y), !anyNA(fold_id),
        setequal(fold_id, seq_len(n_cv)), diff(range(tabulate(fold_id))) <= 1L,
        all(vapply(seq_len(n_cv), function(k)
          min(tabulate(y[fold_id != k]+1L, nbins = 2L)) >= 2L, logical(1L))))
      x <- if (model_index == 1L) target$Z_site else target$W_outcome
      coefficients <- as.numeric(stats::coef(state$fit, s = "lambda.min"))
      stopifnot(identical(coefficients, state$coefficients), all(is.finite(coefficients)))
      predicted <- as.numeric(predict(state$fit, newx = x, s = "lambda.min", type = "response"))
      error <- max(abs(predicted-plogis(drop(cbind(1, x) %*% coefficients))))
      stopifnot(is.finite(error), error < 1e-12)
      rows[[key]] <- data.frame(key, n_training = length(expected_model_ids), n_cv,
        warning_count = length(state$warnings), prediction_error = error,
        minimum_prediction = min(predicted), maximum_prediction = max(predicted))
    }
  }
  rows <- do.call(rbind, rows)
  roce_write_atomic_directory(args[3], function(stage) {
    write.csv(rows, file.path(stage, "state_checks.csv"), row.names = FALSE)
    writeLines(c("all_30_states_checked=TRUE", "training_exclusions_verified=TRUE",
      "coefficient_prediction_identity_verified=TRUE", "inference_validated=FALSE",
      paste0("audit_code_sha256=", roce_sha256_file("diagnosis/tate_common_weight/audit_nested_target_states.R")),
      paste0("states_manifest_sha256=", roce_sha256_file(file.path(args[1], "sha256.txt")))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "nested target state audit")
  print(rows, row.names = FALSE)
}
if (sys.nframe() == 0L) main()
