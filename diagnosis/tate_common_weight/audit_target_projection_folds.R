#!/usr/bin/env Rscript
# Read-only input audit; publish only after all five fold bundles pass.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: audit_target_projection_folds.R INPUT_ROOT OUTPUT")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/target_projection_validation.R")
  root <- normalizePath(args[1], mustWork = TRUE)
  states <- reports <- hashes <- vector("list", 5L)
  payload <- c("candidate_risks.csv", "code_hashes.csv", "metadata.txt",
               "selected_holdout_report.csv", "selection_state.rds", "validation_losses.csv")
  for (outer_fold in 1:5) {
    directory <- file.path(root, paste0("target_projection_fold_", outer_fold, "_v2"))
    manifest <- readLines(file.path(directory, "sha256.txt"))
    payload_names <- substring(manifest, 67L)
    stopifnot(length(payload_names) == length(payload), !anyDuplicated(payload_names), setequal(payload_names, payload))
    actual <- vapply(file.path(directory, payload_names), roce_sha256_file, "")
    stopifnot(identical(unname(actual), substr(manifest, 1L, 64L)))
    hashes[[outer_fold]] <- read.csv(file.path(directory, "code_hashes.csv"))
    if (outer_fold > 1L) stopifnot(identical(hashes[[outer_fold]], hashes[[1L]]))
    state <- readRDS(file.path(directory, "selection_state.rds"))
    report <- read.csv(file.path(directory, "selected_holdout_report.csv"))
    losses <- read.csv(file.path(directory, "validation_losses.csv"))
    selection <- .select_target_projection_candidate(losses, setdiff(1:5, outer_fold))
    stopifnot(isTRUE(all.equal(selection, state$selection, tolerance = 1e-14)), nrow(report) == 1L,
      report$outer_fold == outer_fold, report$selected_fit_succeeded,
      report$selected_candidate == selection$selected, is.null(state$final_fit$failure),
      isTRUE(all.equal(selection$scores, read.csv(file.path(directory, "candidate_risks.csv")), tolerance = 1e-14)))
    stopifnot(length(state$outer_evaluation_ids) == length(state$original_score),
      length(state$adjusted_score) == length(state$original_score),
      all(is.finite(state$original_score)), all(is.finite(state$adjusted_score)),
      abs(mean(state$original_score)-report$original_target_fold_tate) < 1e-12,
      abs(mean(state$adjusted_score)-report$adjusted_target_fold_tate) < 1e-12)
    states[[outer_fold]] <- state
    reports[[outer_fold]] <- report
  }
  ids <- lapply(states, `[[`, "outer_evaluation_ids")
  all_ids <- unlist(ids, use.names = FALSE)
  stopifnot(!anyDuplicated(all_ids), identical(sort(all_ids), seq_along(all_ids)))
  for (outer_fold in 1:5) {
    state <- states[[outer_fold]]
    stopifnot(setequal(state$final_training_ids, setdiff(all_ids, ids[[outer_fold]])),
      !anyDuplicated(state$final_training_ids),
      setequal(names(state$splits), as.character(setdiff(1:5, outer_fold))))
    for (inner in setdiff(1:5, outer_fold)) {
      split <- state$splits[[as.character(inner)]]
      expected_keys <- paste0(c("PS_1_", "OR_1_", "OR_0_"), outer_fold, "_", inner)
      stopifnot(identical(split$nuisance_model_keys, expected_keys),
        identical(split$validation_ids, ids[[inner]]),
        !anyDuplicated(split$training_ids),
        setequal(split$training_ids, setdiff(all_ids, c(ids[[outer_fold]], ids[[inner]]))),
        is.finite(split$prediction_error), split$prediction_error < 1e-12)
    }
  }
  scores <- do.call(rbind, lapply(1:5, function(k) data.frame(
    observation_id = ids[[k]], outer_fold = k,
    original_score = states[[k]]$original_score, adjusted_score = states[[k]]$adjusted_score)))
  scores <- scores[order(scores$observation_id), ]
  reports <- do.call(rbind, reports)
  summary <- data.frame(n_target = nrow(scores),
    original_target_tate = mean(scores$original_score), adjusted_target_tate = mean(scores$adjusted_score),
    inference_validated = FALSE)
  stopifnot(abs(summary$original_target_tate-weighted.mean(reports$original_target_fold_tate, lengths(ids))) < 1e-12,
    abs(summary$adjusted_target_tate-weighted.mean(reports$adjusted_target_fold_tate, lengths(ids))) < 1e-12)
  roce_write_atomic_directory(args[2], function(stage) {
    write.csv(scores, file.path(stage, "observation_scores.csv"), row.names = FALSE)
    write.csv(reports, file.path(stage, "fold_reports.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "target_summary.csv"), row.names = FALSE)
    writeLines(c("all_five_fold_payloads_verified=TRUE", "split_exclusions_verified=TRUE",
      "selection_reproduced=TRUE", "inference_validated=FALSE",
      paste0("audit_code_sha256=", roce_sha256_file("diagnosis/tate_common_weight/audit_target_projection_folds.R")),
      paste0("fold_", 1:5, "_manifest_sha256=", vapply(1:5, function(k)
        roce_sha256_file(file.path(root, paste0("target_projection_fold_", k, "_v2"), "sha256.txt")), ""))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "five-fold target projection audit")
  print(summary, row.names = FALSE)
}
if (sys.nframe() == 0L) main()
