#!/usr/bin/env Rscript

# Controlled grouped-origin CV sensitivity for the saved target fold-4 draw.
# It fits only the two target glmnet CV paths (PS and treated outcome).

audit_target_fold4_grouped_cv_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) {
    stop("Usage: audit_target_fold4_grouped_cv.R DRAW_DIR R_LIBRARY OUTPUT_DIR",
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
    stop("grouped-CV audit output already exists.", call. = FALSE)
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
  fitted <- draw$result$fitted
  if (is.null(reference) || is.null(fitted) ||
      !identical(reference$nuisance_lambda_rule, "min")) {
    stop("requires a completed frozen lambda.min draw.", call. = FALSE)
  }
  sources <- names(reference$weights)
  info <- reference$intermediates$fold_info
  fold_indices <- c(
    list(t = lapply(info, `[[`, "target_idx")),
    setNames(lapply(seq_along(sources), function(j) {
      lapply(info, function(x) x$source_idx[[j]])
    }), sources)
  )
  replay <- roce_refit_resample(
    entry$data_split, fold_indices, seed = draw$summary$resample_seed,
    identity = FALSE
  )
  if (!identical(replay$original_row_ids, draw$result$original_row_ids)) {
    stop("replayed row map differs from saved draw.", call. = FALSE)
  }

  k1 <- 4L
  folds <- replay$precomputed_folds$target_folds
  train <- RoCE:::combine_folds(folds, setdiff(seq_len(5L), k1))
  evaluation <- RoCE:::materialize_fold(folds, k1)
  training_origins <- replay$original_row_ids$t[train$original_idx]
  treated <- which(train$A == 1L)
  model_inputs <- list(
    PS = list(
      x = as.matrix(train$Z_site), y = as.numeric(train$A),
      newx = as.matrix(evaluation$Z_site), origins = training_origins,
      seed = RoCE:::.target_only_cv_seed(k1, NULL, "propensity")
    ),
    OR = list(
      x = as.matrix(train$W_outcome)[treated, , drop = FALSE],
      y = train$Y[treated], newx = as.matrix(evaluation$W_outcome),
      origins = training_origins[treated],
      seed = RoCE:::.target_only_cv_seed(k1, NULL, "outcome", 1L)
    )
  )

  make_grouped_foldid <- function(origins, nfolds, seed) {
    unique_origins <- unique(origins)
    assignment <- RoCE:::with_seed(seed, sample(rep(
      seq_len(nfolds), length.out = length(unique_origins)
    )))
    as.integer(assignment[match(origins, unique_origins)])
  }
  leakage <- function(origins, foldid, scheme, model) {
    split_folds <- split(foldid, origins)
    origin_fold_count <- lengths(lapply(split_folds, unique))
    leaking_origins <- names(origin_fold_count)[origin_fold_count > 1L]
    data.frame(
      model = model,
      scheme = scheme,
      n_rows = length(origins),
      n_unique_origins = length(unique(origins)),
      n_duplicate_copies = length(origins) - length(unique(origins)),
      n_origins_with_multiple_copies = sum(table(origins) > 1L),
      n_origins_crossing_cv_folds = length(leaking_origins),
      n_rows_from_crossing_origins = sum(origins %in% leaking_origins),
      max_origin_multiplicity = max(tabulate(origins)),
      stringsAsFactors = FALSE
    )
  }
  fit_cv <- function(input, model) {
    nfolds <- RoCE:::get_cv_fold_count(nrow(input$x), min_per_fold = 10L)
    default <- RoCE:::with_seed(input$seed, glmnet::cv.glmnet(
      x = input$x, y = input$y, family = "binomial", alpha = 1,
      nfolds = nfolds, nlambda = RoCE:::LAMBDA_GRID_SIZE_STANDARD,
      maxit = RoCE:::GLMNET_MAX_ITER, keep = TRUE
    ))
    grouped_foldid <- make_grouped_foldid(input$origins, nfolds, input$seed)
    grouped <- glmnet::cv.glmnet(
      x = input$x, y = input$y, family = "binomial", alpha = 1,
      foldid = grouped_foldid, nlambda = RoCE:::LAMBDA_GRID_SIZE_STANDARD,
      maxit = RoCE:::GLMNET_MAX_ITER, keep = TRUE
    )
    summarize <- function(object, scheme) {
      coefficients <- as.matrix(stats::coef(
        object, s = "lambda.min"
      ))[, 1L]
      selected <- which(abs(coefficients) > 0)
      min_index <- which.min(abs(object$lambda - object$lambda.min))
      prediction <- as.numeric(predict(
        object, newx = input$newx, s = "lambda.min", type = "response"
      ))
      data.frame(
        model = model,
        scheme = scheme,
        n_train = nrow(input$x),
        n_cv_folds = nfolds,
        lambda_min = object$lambda.min,
        lambda_1se = object$lambda.1se,
        cv_loss_at_lambda_min = object$cvm[[min_index]],
        cv_se_at_lambda_min = object$cvsd[[min_index]],
        nonzero_slopes = sum(selected != 1L),
        max_abs_selected_coefficient = max(abs(coefficients[selected])),
        row295_prediction = prediction[
          match(773L, evaluation$original_idx)
        ],
        stringsAsFactors = FALSE
      )
    }
    list(
      summary = rbind(
        summarize(default, "row_level_default"),
        summarize(grouped, "grouped_by_origin")
      ),
      leakage = rbind(
        leakage(input$origins, default$foldid, "row_level_default", model),
        leakage(input$origins, grouped_foldid, "grouped_by_origin", model)
      ),
      default = default,
      grouped = grouped
    )
  }
  fitted_cv <- lapply(names(model_inputs), function(model) {
    fit_cv(model_inputs[[model]], model)
  })
  names(fitted_cv) <- names(model_inputs)
  tuning <- do.call(rbind, lapply(fitted_cv, `[[`, "summary"))
  leakage_table <- do.call(rbind, lapply(fitted_cv, `[[`, "leakage"))
  if (any(leakage_table$n_origins_crossing_cv_folds[
    leakage_table$scheme == "grouped_by_origin"
  ] != 0L)) stop("grouped-origin fold IDs still leak origins.", call. = FALSE)

  # keep=TRUE must reproduce the already audited default row-level fits.
  prior_tuning <- read.csv(file.path(
    dirname(draw_dir), "..", "full_refit_calibration_v18_audits",
    paste0(basename(draw_dir), "_target_fold4"), "tuning.csv"
  ), stringsAsFactors = FALSE)
  prior_tuning <- prior_tuning[prior_tuning$sample == "resampled", ]
  default_rows <- tuning[tuning$scheme == "row_level_default", ]
  default_rows <- default_rows[match(prior_tuning$model, default_rows$model), ]
  default_error <- max(abs(c(
    default_rows$lambda_min - prior_tuning$lambda_min,
    default_rows$lambda_1se - prior_tuning$lambda_1se,
    default_rows$cv_loss_at_lambda_min - prior_tuning$cv_loss_at_lambda_min,
    default_rows$nonzero_slopes - prior_tuning$selected_nonzero_slopes
  )))
  saved_ps <- fitted$arm_results$mu1$fold_results[[k1]]$target_only$prop_scores
  saved_or <- fitted$arm_results$mu1$fold_results[[k1]]$target_only$m_pred
  default_prediction_error <- max(abs(c(
    default_rows$row295_prediction[default_rows$model == "PS"] -
      saved_ps[match(773L, evaluation$original_idx)],
    default_rows$row295_prediction[default_rows$model == "OR"] -
      saved_or[match(773L, evaluation$original_idx)]
  )))
  if (!is.finite(default_error) || default_error > 1e-12 ||
      default_prediction_error > 1e-12) {
    stop("keep=TRUE changed the default row-level tuning or predictions.",
         call. = FALSE)
  }
  row_index <- match(773L, evaluation$original_idx)
  row_sensitivity <- do.call(rbind, lapply(
    c("row_level_default", "grouped_by_origin"), function(scheme) {
      propensity <- tuning$row295_prediction[
        tuning$model == "PS" & tuning$scheme == scheme
      ]
      m1 <- tuning$row295_prediction[
        tuning$model == "OR" & tuning$scheme == scheme
      ]
      data.frame(
        scheme = scheme,
        resampled_position = 773L,
        origin_id = replay$original_row_ids$t[[773L]],
        A = evaluation$A[[row_index]],
        Y = evaluation$Y[[row_index]],
        propensity = propensity,
        m1 = m1,
        aipw = m1 + (evaluation$Y[[row_index]] - m1) / propensity,
        stringsAsFactors = FALSE
      )
    }
  ))

  roce_write_atomic_directory(output_dir, function(staging) {
    write.csv(tuning, file.path(staging, "tuning.csv"), row.names = FALSE)
    write.csv(leakage_table, file.path(staging, "cv_origin_leakage.csv"),
              row.names = FALSE)
    write.csv(row_sensitivity, file.path(staging, "row295_sensitivity.csv"),
              row.names = FALSE)
    writeLines(c(
      "target_fold4_grouped_cv_sensitivity=passed",
      "scope=same resampled target PS/OR matrices; foldid changed only",
      paste0("default_keep_max_error=", default_error),
      paste0("default_prediction_max_error=", default_prediction_error),
      "grouped_origin_leakage=0",
      "causal_attribution_to_duplicate_leakage=false"
    ), file.path(staging, "audit_passed.txt"))
    payloads <- c(
      "tuning.csv", "cv_origin_leakage.csv", "row295_sensitivity.csv",
      "audit_passed.txt"
    )
    hashes <- vapply(file.path(staging, payloads), roce_sha256_file,
                     character(1L))
    writeLines(paste(hashes, payloads, sep = "  "),
               file.path(staging, "sha256.txt"))
  }, caller = "target fold-4 grouped-origin CV sensitivity")
  message("[passed] target fold-4 grouped-origin CV sensitivity: ", output_dir)
  invisible(list(tuning = tuning, leakage = leakage_table))
}

if (sys.nframe() == 0L) audit_target_fold4_grouped_cv_main()
