#!/usr/bin/env Rscript

# Isolated target-only fold-4 nuisance diagnostic. No source models are fit.

audit_target_fold4_tuning_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop("Usage: audit_target_fold4_tuning.R DRAW_DIR R_LIBRARY OUTPUT_DIR",
         call. = FALSE)
  }
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", ROCE_NUISANCE_CV_THREADS = "1"
  )
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/full_refit_resampling.R")
  draw_dir <- normalizePath(args[[1L]], mustWork = TRUE)
  project_library <- normalizePath(args[[2L]], mustWork = TRUE)
  output_dir <- args[[3L]]
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("target-fold audit output already exists.", call. = FALSE)
  }
  .libPaths(c(project_library, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  if (!identical(
    normalizePath(find.package("RoCE"), mustWork = TRUE),
    normalizePath(file.path(project_library, "RoCE"), mustWork = TRUE)
  )) stop("wrong RoCE installation loaded.", call. = FALSE)

  draw <- readRDS(file.path(draw_dir, "draw.rds"))
  source_saved <- readRDS(file.path(draw$source_bundle, "artifacts.rds"))
  rho_key <- format(draw$summary$rho, scientific = FALSE, trim = TRUE)
  entry <- source_saved$group_result$artifacts[[rho_key]]
  reference <- entry$direct_tate_results$one_round_crossfit
  refitted <- draw$result$fitted
  if (is.null(reference) || is.null(refitted) ||
      !identical(reference$n_folds, 5L) ||
      !identical(reference$nuisance_lambda_rule, "min")) {
    stop("diagnostic requires the frozen five-fold lambda.min fit.",
         call. = FALSE)
  }

  sources <- names(reference$weights)
  info <- reference$intermediates$fold_info
  fold_indices <- c(
    list(t = lapply(info, `[[`, "target_idx")),
    setNames(lapply(seq_along(sources), function(j) {
      lapply(info, function(x) x$source_idx[[j]])
    }), sources)
  )
  original <- roce_refit_resample(
    entry$data_split, fold_indices, seed = draw$summary$resample_seed,
    identity = TRUE
  )
  replay <- roce_refit_resample(
    entry$data_split, fold_indices, seed = draw$summary$resample_seed,
    identity = FALSE
  )
  if (!identical(replay$original_row_ids, draw$result$original_row_ids)) {
    stop("replayed resampling map does not match the saved draw.",
         call. = FALSE)
  }

  k1 <- 4L
  A_val <- 1L
  fit_one <- function(label, folds, saved_fit, requested_position) {
    complement <- RoCE::estimate_target_only_from_complement(
      target_folds = folds,
      k1 = k1,
      n_folds = 5L,
      family = "binomial",
      A_val = A_val,
      propensity_cache = NULL,
      nuisance_lambda_rule = "min"
    )
    train <- RoCE:::combine_folds(folds, setdiff(seq_len(5L), k1))
    evaluation <- RoCE:::materialize_fold(folds, k1)
    X_ps <- as.matrix(train$Z_site)
    X_or <- as.matrix(train$W_outcome)
    X_ps_eval <- as.matrix(evaluation$Z_site)
    X_or_eval <- as.matrix(evaluation$W_outcome)

    ps_cv <- RoCE:::with_seed(
      RoCE:::.target_only_cv_seed(k1, NULL, "propensity"),
      glmnet::cv.glmnet(
        x = X_ps, y = as.numeric(train$A), family = "binomial",
        alpha = 1,
        nfolds = RoCE:::get_cv_fold_count(nrow(X_ps), min_per_fold = 10L),
        nlambda = RoCE:::LAMBDA_GRID_SIZE_STANDARD,
        maxit = RoCE:::GLMNET_MAX_ITER
      )
    )
    treated <- which(train$A == A_val)
    or_cv <- RoCE:::with_seed(
      RoCE:::.target_only_cv_seed(k1, NULL, "outcome", A_val),
      glmnet::cv.glmnet(
        x = X_or[treated, , drop = FALSE], y = train$Y[treated],
        family = "binomial", alpha = 1,
        nfolds = RoCE:::get_cv_fold_count(length(treated), min_per_fold = 10L),
        nlambda = RoCE:::LAMBDA_GRID_SIZE_STANDARD,
        maxit = RoCE:::GLMNET_MAX_ITER
      )
    )
    ps_direct <- RoCE:::clip_propensity(as.numeric(predict(
      ps_cv, newx = X_ps_eval, s = "lambda.min", type = "response"
    )))
    or_direct <- RoCE:::clip_outcome_pred(as.numeric(predict(
      or_cv, newx = X_or_eval, s = "lambda.min", type = "response"
    )), "binomial")
    phi_direct <- RoCE:::calculate_aipw_pseudo_outcome(
      as.numeric(evaluation$Y), evaluation$A, or_direct, ps_direct,
      A_val = A_val
    )
    saved_phi <- saved_fit$varphi_ot + saved_fit$estimate
    prediction_errors <- c(
      complement_ps = max(abs(complement$prop_scores - ps_direct)),
      complement_or = max(abs(complement$m_pred - or_direct)),
      complement_phi = max(abs(
        complement$varphi_ot + complement$estimate - phi_direct
      )),
      saved_ps = max(abs(saved_fit$prop_scores - ps_direct)),
      saved_or = max(abs(saved_fit$m_pred - or_direct)),
      saved_phi = max(abs(saved_phi - phi_direct))
    )
    if (any(!is.finite(prediction_errors)) || max(prediction_errors) > 1e-10) {
      stop(label, " predictions do not reproduce saved/complement vectors.",
           call. = FALSE)
    }
    eval_index <- match(requested_position, evaluation$original_idx)
    if (is.na(eval_index)) stop("requested row is absent from fold 4.")

    cv_row <- function(model, object) {
      coefficient <- as.matrix(stats::coef(object, s = "lambda.min"))[, 1L]
      selected <- which(abs(coefficient) > 0)
      min_index <- which.min(abs(object$lambda - object$lambda.min))
      one_se_index <- which.min(abs(object$lambda - object$lambda.1se))
      data.frame(
        sample = label,
        model = model,
        n_train = if (model == "PS") nrow(X_ps) else length(treated),
        n_cv_folds = if (model == "PS") {
          RoCE:::get_cv_fold_count(nrow(X_ps), min_per_fold = 10L)
        } else {
          RoCE:::get_cv_fold_count(length(treated), min_per_fold = 10L)
        },
        lambda_min = object$lambda.min,
        lambda_1se = object$lambda.1se,
        cv_loss_at_lambda_min = object$cvm[[min_index]],
        cv_se_at_lambda_min = object$cvsd[[min_index]],
        cv_loss_at_lambda_1se = object$cvm[[one_se_index]],
        selected_nonzero_including_intercept = length(selected),
        selected_nonzero_slopes = sum(selected != 1L),
        selected_coefficient_names = paste(names(coefficient)[selected],
                                           collapse = ";"),
        max_abs_selected_coefficient = max(abs(coefficient[selected])),
        stringsAsFactors = FALSE
      )
    }
    list(
      tuning = rbind(cv_row("PS", ps_cv), cv_row("OR", or_cv)),
      row = data.frame(
        sample = label,
        requested_position = requested_position,
        evaluation_vector_index = eval_index,
        origin_id = if (label == "original") requested_position else
          replay$original_row_ids$t[[requested_position]],
        A = evaluation$A[[eval_index]],
        Y = evaluation$Y[[eval_index]],
        propensity = ps_direct[[eval_index]],
        m1 = or_direct[[eval_index]],
        aipw = phi_direct[[eval_index]],
        stringsAsFactors = FALSE
      ),
      errors = data.frame(
        sample = label,
        quantity = names(prediction_errors),
        max_abs_error = as.numeric(prediction_errors),
        stringsAsFactors = FALSE
      )
    )
  }

  original_saved <- reference$arm_results$mu1$fold_results[[k1]]$target_only
  refitted_saved <- refitted$arm_results$mu1$fold_results[[k1]]$target_only
  original_result <- fit_one(
    "original", original$precomputed_folds$target_folds,
    original_saved, 295L
  )
  refitted_result <- fit_one(
    "resampled", replay$precomputed_folds$target_folds,
    refitted_saved, 773L
  )
  tuning <- rbind(original_result$tuning, refitted_result$tuning)
  row_scores <- rbind(original_result$row, refitted_result$row)
  prediction_errors <- rbind(original_result$errors, refitted_result$errors)
  if (!identical(row_scores$origin_id, c(295L, 295L))) {
    stop("the requested original/resampled rows do not share origin 295.",
         call. = FALSE)
  }

  summary <- data.frame(
    quantity = c(
      "PS lambda.min", "PS lambda.1se", "PS CV loss at min",
      "PS nonzero slopes", "OR lambda.min", "OR lambda.1se",
      "OR CV loss at min", "OR nonzero slopes", "row propensity",
      "row m1", "row AIPW"
    ),
    original = c(
      tuning$lambda_min[tuning$sample == "original" & tuning$model == "PS"],
      tuning$lambda_1se[tuning$sample == "original" & tuning$model == "PS"],
      tuning$cv_loss_at_lambda_min[
        tuning$sample == "original" & tuning$model == "PS"
      ],
      tuning$selected_nonzero_slopes[
        tuning$sample == "original" & tuning$model == "PS"
      ],
      tuning$lambda_min[tuning$sample == "original" & tuning$model == "OR"],
      tuning$lambda_1se[tuning$sample == "original" & tuning$model == "OR"],
      tuning$cv_loss_at_lambda_min[
        tuning$sample == "original" & tuning$model == "OR"
      ],
      tuning$selected_nonzero_slopes[
        tuning$sample == "original" & tuning$model == "OR"
      ],
      row_scores$propensity[row_scores$sample == "original"],
      row_scores$m1[row_scores$sample == "original"],
      row_scores$aipw[row_scores$sample == "original"]
    ),
    resampled = c(
      tuning$lambda_min[tuning$sample == "resampled" & tuning$model == "PS"],
      tuning$lambda_1se[tuning$sample == "resampled" & tuning$model == "PS"],
      tuning$cv_loss_at_lambda_min[
        tuning$sample == "resampled" & tuning$model == "PS"
      ],
      tuning$selected_nonzero_slopes[
        tuning$sample == "resampled" & tuning$model == "PS"
      ],
      tuning$lambda_min[tuning$sample == "resampled" & tuning$model == "OR"],
      tuning$lambda_1se[tuning$sample == "resampled" & tuning$model == "OR"],
      tuning$cv_loss_at_lambda_min[
        tuning$sample == "resampled" & tuning$model == "OR"
      ],
      tuning$selected_nonzero_slopes[
        tuning$sample == "resampled" & tuning$model == "OR"
      ],
      row_scores$propensity[row_scores$sample == "resampled"],
      row_scores$m1[row_scores$sample == "resampled"],
      row_scores$aipw[row_scores$sample == "resampled"]
    ),
    stringsAsFactors = FALSE
  )
  summary$difference <- summary$resampled - summary$original

  roce_write_atomic_directory(output_dir, function(staging) {
    write.csv(summary, file.path(staging, "summary.csv"), row.names = FALSE)
    write.csv(tuning, file.path(staging, "tuning.csv"), row.names = FALSE)
    write.csv(row_scores, file.path(staging, "row295_scores.csv"),
              row.names = FALSE)
    write.csv(prediction_errors, file.path(staging, "prediction_identity.csv"),
              row.names = FALSE)
    writeLines(c(
      "target_fold4_tuning_audit=passed",
      "scope=target-only fold4 A_val1; no source or full RoCE refit",
      paste0("draw_bundle=", draw_dir),
      paste0("max_prediction_abs_error=",
             format(max(prediction_errors$max_abs_error), scientific = TRUE)),
      "nuisance_cv_duplicate_origins_are_possible=true",
      "causal_attribution_to_duplicate_leakage=false"
    ), file.path(staging, "audit_passed.txt"))
    payloads <- c(
      "summary.csv", "tuning.csv", "row295_scores.csv",
      "prediction_identity.csv", "audit_passed.txt"
    )
    hashes <- vapply(file.path(staging, payloads), roce_sha256_file,
                     character(1L))
    writeLines(paste(hashes, payloads, sep = "  "),
               file.path(staging, "sha256.txt"))
  }, caller = "target fold-4 tuning audit")
  message("[passed] target fold-4 tuning audit: ", output_dir)
  invisible(summary)
}

if (sys.nframe() == 0L) audit_target_fold4_tuning_main()
