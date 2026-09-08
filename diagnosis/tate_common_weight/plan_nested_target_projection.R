#!/usr/bin/env Rscript
# Dependency audit, not a nuisance fit or a change to the estimator.
.nested_target_projection_plan <- function(n_folds = 5L) {
  if (!identical(n_folds, 5L)) stop("this development protocol requires five folds")
  rows <- list()
  for (outer_fold in seq_len(n_folds)) {
    for (evaluation_fold in setdiff(seq_len(n_folds), outer_fold)) {
      for (validation_fold in setdiff(seq_len(n_folds), c(outer_fold, evaluation_fold))) {
        excluded <- sort(c(outer_fold, evaluation_fold, validation_fold))
        training <- setdiff(seq_len(n_folds), excluded)
        rows[[length(rows)+1L]] <- data.frame(outer_fold, evaluation_fold, validation_fold,
          excluded_folds = paste(excluded, collapse = ":"),
          training_folds = paste(training, collapse = ":"),
          nuisance_training_key = paste0("exclude_", paste(excluded, collapse = "_")))
      }
    }
  }
  do.call(rbind, rows)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: plan_nested_target_projection.R FIVE_FOLD_ROOT OUTPUT")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  plan <- .nested_target_projection_plan()
  audit <- list()
  for (outer_fold in 1:5) {
    root <- file.path(args[1], paste0("target_projection_fold_", outer_fold, "_v2"))
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == "selection_state.rds"], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, "selection_state.rds")) == expected)
    state <- readRDS(file.path(root, "selection_state.rds"))
    for (evaluation_fold in setdiff(1:5, outer_fold)) {
      evaluation_ids <- state$splits[[as.character(evaluation_fold)]]$validation_ids
      validation_overlap <- sum(vapply(state$splits, function(split)
        sum(evaluation_ids %in% split$validation_ids), integer(1L)))
      training_overlap <- sum(vapply(state$splits, function(split)
        sum(evaluation_ids %in% split$training_ids), integer(1L)))
      stopifnot(validation_overlap == length(evaluation_ids),
                training_overlap == 3L*length(evaluation_ids))
      audit[[length(audit)+1L]] <- data.frame(outer_fold, evaluation_fold,
        n_evaluation = length(evaluation_ids),
        outer_selection_validation_overlap = validation_overlap,
        outer_selection_training_overlap = training_overlap,
        outer_selected_scale_is_inner_heldout = FALSE)
    }
  }
  audit <- do.call(rbind, audit)
  stopifnot(nrow(plan) == 60L, length(unique(plan$nuisance_training_key)) == 10L)
  roce_write_atomic_directory(args[2], function(stage) {
    write.csv(plan, file.path(stage, "nested_split_plan.csv"), row.names = FALSE)
    write.csv(audit, file.path(stage, "outer_selection_reuse_audit.csv"), row.names = FALSE)
    writeLines(c("status=dependency_audit_only", "new_nuisance_fits=0",
      "ordered_validation_splits=60", "unique_training_subsets=10",
      "potential_shared_nuisance_models=30", "shared_fit_seed_policy_not_yet_frozen=TRUE",
      "existing_estimator_bug_established=FALSE", "inference_validated=FALSE",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/plan_nested_target_projection.R"))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "nested target projection dependency audit")
  cat("Verified 20 outer/evaluation pairs; planned 60 validation splits on 10 training subsets.\n")
}
if (sys.nframe() == 0L) main()
