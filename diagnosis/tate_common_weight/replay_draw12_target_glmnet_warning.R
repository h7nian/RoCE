#!/usr/bin/env Rscript

# Bounded diagnostic replay of target PS/OR CV calls only. This does not run a
# source nuisance fit, aggregation, full refit, or modify the frozen package.

library(RoCE)
source("scripts/slurm/result_provenance.R")

bundle <- "results/direct_tate_mc500_b5000/full_refit_grouped_cv_v19/seed_000001_rho_0_draw_0012/draw.rds"
draw <- readRDS(bundle)
source_saved <- readRDS(file.path(draw$source_bundle, "artifacts.rds"))
entry <- source_saved$group_result$artifacts[["0"]]
reference <- entry$direct_tate_results$one_round_crossfit
rows <- draw$result$original_row_ids$t
target <- entry$data_split$t
for (field in names(target)) {
  if (field == "n") next
  value <- target[[field]]
  target[[field]] <- if (is.matrix(value)) value[rows, , drop = FALSE] else value[rows]
}
target$cv_group_id <- rows
folds <- lapply(reference$intermediates$fold_info, function(info) {
  list(original_idx = info$target_idx, n = length(info$target_idx))
})
attr(folds, ".data_ref") <- target

run_cv <- function(model, arm, k1, k2, data, eval_data, expected_prediction) {
  x <- if (model == "PS") data$Z_site else data$W_outcome[data$A == arm, , drop = FALSE]
  y <- if (model == "PS") data$A else data$Y[data$A == arm]
  groups <- if (model == "PS") data$cv_group_id else data$cv_group_id[data$A == arm]
  outcome_spec <- RoCE:::resolve_glm_family(reference$family)
  family <- if (model == "PS") "binomial" else outcome_spec$glmnet_family
  n_folds <- RoCE:::get_cv_fold_count(nrow(x), min_per_fold = 10L)
  caller <- paste("estimate_complement_fold_aipw", model)
  seed <- RoCE:::.target_only_cv_seed(
    k1, k2, if (model == "PS") "propensity" else "outcome",
    if (model == "OR") arm else 0L
  )
  warnings <- character()
  fit <- RoCE:::with_seed(seed, withCallingHandlers({
    fold_id <- RoCE:::.make_nuisance_cv_fold_id(groups, n_folds, caller)
    glmnet::cv.glmnet(
      x = x, y = y, family = family, alpha = 1,
      foldid = fold_id, nlambda = 100L, maxit = RoCE:::GLMNET_MAX_ITER
    )
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
  }))
  prediction <- as.numeric(predict(
    fit,
    newx = if (model == "PS") eval_data$Z_site else eval_data$W_outcome,
    s = "lambda.min", type = "response"
  ))
  prediction <- if (model == "PS") {
    RoCE:::clip_propensity(prediction)
  } else {
    RoCE:::clip_outcome_pred(prediction, reference$family)
  }
  matching_records <- Filter(function(record) {
    identical(record$caller, caller) && identical(record$cv_group_id, groups) &&
      identical(record$cv_fold_id, fold_id)
  }, draw$result$nuisance_cv_partitions)
  if (length(matching_records) != 1L) {
    stop(sprintf(
      "saved partition match count=%d for model=%s arm=%d k1=%d k2=%s",
      length(matching_records), model, arm, k1,
      if (is.null(k2)) "NA" else as.character(k2)
    ))
  }
  data.frame(
    model = model, arm = arm, k1 = k1,
    k2 = if (is.null(k2)) NA_integer_ else k2,
    n_train = nrow(x), n_cv_folds = n_folds,
    lambda_path_length = length(fit$lambda),
    lambda_min = fit$lambda.min, lambda_1se = fit$lambda.1se,
    selected_lambda = fit$lambda.min,
    prediction_length = length(prediction),
    saved_prediction_length = length(expected_prediction),
    saved_prediction_max_error = max(abs(prediction - expected_prediction)),
    saved_partition_match_count = length(matching_records),
    saved_fold_id_max_error = max(abs(
      fold_id - matching_records[[1L]]$cv_fold_id
    )),
    warning_count = length(warnings),
    warning = paste(unique(warnings), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

rows_out <- list()
index <- 0L
for (arm in c(1L, 0L)) {
  for (k1 in seq_len(5L)) {
    for (k2 in c(NA_integer_, setdiff(seq_len(5L), k1))) {
      excluded <- if (is.na(k2)) k1 else c(k1, k2)
      training <- RoCE:::combine_folds(folds, setdiff(seq_len(5L), excluded))
      eval_fold <- if (is.na(k2)) k1 else k2
      eval_data <- RoCE:::materialize_fold(folds, eval_fold)
      saved_target <- if (is.na(k2)) {
        draw$result$fitted$arm_results[[paste0("mu", arm)]]$
          fold_results[[k1]]$target_only
      } else {
        draw$result$fitted$arm_results[[paste0("mu", arm)]]$
          fold_results[[k1]]$target_only_inner[[paste0("k2_", k2)]]
      }
      models <- if (arm == 0L) "OR" else c("PS", "OR")
      for (model in models) {
        index <- index + 1L
        rows_out[[index]] <- run_cv(
          model, arm, k1, if (is.na(k2)) NULL else k2, training, eval_data,
          if (model == "PS") saved_target$prop_scores else saved_target$m_pred
        )
      }
    }
  }
}
result <- do.call(rbind, rows_out)
output <- "results/direct_tate_mc500_b5000/full_refit_grouped_cv_v19_audits/draw12_target_glmnet_warning_replay_v2"
if (file.exists(output) || dir.exists(output)) stop("diagnostic output exists")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(output, "target_cv_replay.csv"), row.names = FALSE)
writeLines(c(
  "target_cv_warning_replay=completed",
  paste0("reference_family=", reference$family),
  paste0("script_fingerprint=", roce_sha256_file(
    "diagnosis/tate_common_weight/replay_draw12_target_glmnet_warning.R"
  ))
), file.path(output, "metadata.txt"))
files <- c("target_cv_replay.csv", "metadata.txt")
writeLines(paste(
  vapply(file.path(output, files), roce_sha256_file, character(1L)),
  files, sep = "  "
), file.path(output, "sha256.txt"))
print(result[result$warning_count > 0L, ], row.names = FALSE)
cat("calls=", nrow(result), " warnings=", sum(result$warning_count), "\n", sep = "")
