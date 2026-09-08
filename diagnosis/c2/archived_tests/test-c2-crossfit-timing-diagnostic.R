.c2_crossfit_timing_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_CROSSFIT_TIMING", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-crossfit-timing", filter, fixed = TRUE)
}

.c2_timing_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_timing_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) default else parsed
}

.c2_timing_sources_env <- function(available) {
  value <- Sys.getenv("ROCE_C2_CROSSFIT_TIMING_SOURCES", unset = "")
  if (!nzchar(value)) return(available[1L])
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  bad <- setdiff(parsed, available)
  if (length(bad) > 0L) {
    stop(sprintf(
      "ROCE_C2_CROSSFIT_TIMING_SOURCES contains unavailable source(s): %s; available: %s",
      paste(bad, collapse = ", "), paste(available, collapse = ", ")
    ), call. = FALSE)
  }
  parsed
}

.c2_timing_log <- function(path, state, detail) {
  line <- sprintf("[%s] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  state, detail)
  cat(line, "\n", sep = "")
  flush.console()
  cat(line, "\n", file = path, append = TRUE, sep = "")
}

.c2_timing_stage <- function(path, rows_env, stage, detail, expr,
                             extra_done = function(value) "") {
  .c2_timing_log(path, "stage_start", sprintf("stage=%s %s", stage, detail))
  start <- Sys.time()
  value <- force(expr)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  done_extra <- extra_done(value)
  .c2_timing_log(
    path, "stage_done",
    sprintf("stage=%s elapsed_sec=%.3f %s %s", stage, elapsed, detail, done_extra)
  )
  rows_env$rows[[length(rows_env$rows) + 1L]] <- data.frame(
    stage = stage,
    detail = detail,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
  value
}

.c2_timing_support <- function(coef) {
  coef <- as.numeric(coef)
  if (length(coef) <= 1L) return(0L)
  sum(abs(coef[-1L]) > 1e-8)
}

.c2_timing_lambda_detail <- function(coef) {
  sprintf(
    "lambda_used=%.6g lambda_min=%.6g lambda_1se=%.6g support=%d",
    as.numeric(attr(coef, "lambda_used") %||% NA_real_),
    as.numeric(attr(coef, "lambda_min") %||% NA_real_),
    as.numeric(attr(coef, "lambda_1se") %||% NA_real_),
    .c2_timing_support(coef)
  )
}

test_that("c2-crossfit-timing localizes one-round source-site bottlenecks", {
  skip_if_not(.c2_crossfit_timing_enabled(),
              message = "set ROCE_RUN_C2_CROSSFIT_TIMING=1 or run with --filter c2-crossfit-timing")

  n_total <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_N", 5000L)
  K_sites <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_K", 5L)
  p <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_P", 50L)
  n_folds <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_FOLDS", 3L)
  nlambda_init <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_NLAMBDA_INIT", 20L)
  seed <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_SEED", 270001L)
  k1 <- .c2_timing_int_env("ROCE_C2_CROSSFIT_TIMING_K1", 1L)
  shift_strength <- .c2_timing_num_env("ROCE_C2_CROSSFIT_TIMING_SHIFT", 0.5)
  nuisance_rule <- Sys.getenv("ROCE_C2_CROSSFIT_TIMING_NUISANCE_RULE", "min")
  nuisance_rule <- match.arg(nuisance_rule, c("min", "1se"))

  out_dir <- file.path("c2_crossfit_timing_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  progress_path <- file.path(out_dir, paste0(ts, "_c2_crossfit_timing_progress.log"))
  timing_path <- file.path(out_dir, paste0(ts, "_c2_crossfit_timing_stages.csv"))
  rows_env <- new.env(parent = emptyenv())
  rows_env$rows <- list()

  .c2_timing_log(
    progress_path, "config",
    sprintf(
      "n_total=%d K=%d p=%d folds=%d nlambda_init=%d seed=%d k1=%d nuisance_rule=%s",
      n_total, K_sites, p, n_folds, nlambda_init, seed, k1, nuisance_rule
    )
  )

  set.seed(seed)
  data <- .c2_timing_stage(progress_path, rows_env, "generate_data", "", {
    generate_simulation_data(
      n_total = n_total,
      K = K_sites,
      p = p,
      config = "C2",
      estimand_type = "superpopulation",
      site_allocation = "model",
      transform_type = "mild",
      outcome_type = "binary",
      heterogeneity_type = "none",
      shift_strength = shift_strength,
      dgp_type = "roce",
      warn_ignored = FALSE
    )
  })

  data_split <- split_data_by_site(data)
  source_sites <- setdiff(names(data_split), "t")
  selected_sources <- .c2_timing_sources_env(source_sites)
  if (k1 > n_folds) {
    stop(sprintf("ROCE_C2_CROSSFIT_TIMING_K1=%d exceeds n_folds=%d", k1, n_folds),
         call. = FALSE)
  }

  glm_spec <- resolve_glm_family("binomial")
  family_int <- glm_spec$family_int
  link_int <- glm_spec$link_int
  A_val <- 1L

  folds <- .c2_timing_stage(progress_path, rows_env, "build_folds", "", {
    build_crossfit_folds(data_split, n_folds)
  })
  target_folds <- folds$target_folds
  source_folds <- folds$source_folds
  secondary_folds <- setdiff(seq_len(n_folds), k1)

  target_initial_models <- list()
  target_summaries <- list()
  cached_lambda <- NULL

  for (k2 in secondary_folds) {
    k2_key <- paste0("k2_", k2)
    training_folds <- setdiff(seq_len(n_folds), c(k1, k2))
    target_train <- .c2_timing_stage(
      progress_path, rows_env, "target_combine_folds",
      sprintf("k1=%d k2=%d train_folds=%s", k1, k2, paste(training_folds, collapse = ",")),
      combine_folds(target_folds, training_folds)
    )
    target_calib_k2 <- .c2_timing_stage(
      progress_path, rows_env, "target_materialize_calib",
      sprintf("k1=%d k2=%d", k1, k2),
      materialize_fold(target_folds, k2)
    )

    alpha_init <- .c2_timing_stage(
      progress_path, rows_env, "target_fit_initial_outcome",
      sprintf("k1=%d k2=%d n_train=%d n_arm=%d cached=%s",
              k1, k2, target_train$n, sum(target_train$A == A_val),
              if (is.null(cached_lambda)) "no" else "yes"),
      fit_initial_outcome(
        target_train$W_outcome, target_train$Y, target_train$A, A_val,
        lambda = cached_lambda,
        nlambda = nlambda_init,
        family = "binomial",
        lambda_rule = nuisance_rule
      ),
      extra_done = .c2_timing_lambda_detail
    )
    if (is.null(cached_lambda)) {
      cached_lambda <- attr(alpha_init, "lambda_used")
    }
    target_initial_models[[k2_key]] <- alpha_init

    mean_grad_psi_init <- .c2_timing_stage(
      progress_path, rows_env, "target_mean_grad",
      sprintf("k1=%d k2=%d n_calib=%d", k1, k2, target_calib_k2$n),
      .mean_glm_gradient_site_basis(
        target_calib_k2$W_outcome, target_calib_k2$Z_site,
        alpha_init, family_int, link_int
      )
    )
    target_summaries[[k2_key]] <- list(
      mean_grad_psi_init = mean_grad_psi_init,
      mean_phi = c(1, colMeans(target_train$Z_site))
    )
  }

  for (s in selected_sources) {
    .c2_timing_log(progress_path, "source_start", sprintf("source=%s", s))
    gamma_init_list <- list()
    alpha_init_list <- list()
    source_calib_list <- list()
    mean_grad_psi_list <- list()
    cached_lambda_init_dr <- NULL
    prev_gamma_init <- NULL

    for (k2 in secondary_folds) {
      k2_key <- paste0("k2_", k2)
      training_folds <- setdiff(seq_len(n_folds), c(k1, k2))
      source_calib <- .c2_timing_stage(
        progress_path, rows_env, "source_materialize_calib",
        sprintf("source=%s k1=%d k2=%d", s, k1, k2),
        materialize_fold(source_folds[[s]], k2)
      )
      source_train <- .c2_timing_stage(
        progress_path, rows_env, "source_combine_folds",
        sprintf("source=%s k1=%d k2=%d train_folds=%s",
                s, k1, k2, paste(training_folds, collapse = ",")),
        combine_folds(source_folds[[s]], training_folds)
      )

      gamma_init <- .c2_timing_stage(
        progress_path, rows_env, "source_fit_initial_density_ratio",
        sprintf("source=%s k1=%d k2=%d n_train=%d n_arm=%d cached=%s",
                s, k1, k2, source_train$n, sum(source_train$A == A_val),
                if (is.null(cached_lambda_init_dr)) "no" else "yes"),
        fit_initial_density_ratio(
          source_train$Z_site, source_train$A,
          target_summaries[[k2_key]]$mean_phi,
          lambda = cached_lambda_init_dr,
          A_val = A_val,
          warm_start = prev_gamma_init,
          lambda_rule = nuisance_rule
        ),
        extra_done = .c2_timing_lambda_detail
      )
      if (is.null(cached_lambda_init_dr)) {
        cached_lambda_init_dr <- attr(gamma_init, "lambda_used")
      }
      prev_gamma_init <- as.numeric(gamma_init)

      gamma_init_list[[k2_key]] <- gamma_init
      alpha_init_list[[k2_key]] <- target_initial_models[[k2_key]]
      source_calib_list[[k2_key]] <- source_calib
      mean_grad_psi_list[[k2_key]] <- target_summaries[[k2_key]]$mean_grad_psi_init
    }

    Z_cal_stack <- .c2_timing_stage(
      progress_path, rows_env, "source_stack_Z",
      sprintf("source=%s k1=%d", s, k1),
      .stack_fold_field(source_calib_list, "Z_site", "c2-crossfit-timing")
    )
    W_cal_stack <- .c2_timing_stage(
      progress_path, rows_env, "source_stack_W",
      sprintf("source=%s k1=%d", s, k1),
      .stack_fold_field(source_calib_list, "W_outcome", "c2-crossfit-timing")
    )
    A_cal_stack <- .stack_fold_field(source_calib_list, "A", "c2-crossfit-timing")
    Y_cal_stack <- .stack_fold_field(source_calib_list, "Y", "c2-crossfit-timing")
    mean_grad_psi_avg <- .average_numeric_list(mean_grad_psi_list, "c2-crossfit-timing")

    W_plugin_block <- .c2_timing_stage(
      progress_path, rows_env, "source_make_W_plugin_block",
      sprintf("source=%s k1=%d", s, k1),
      .make_plugin_block_design(
        lapply(source_calib_list, `[[`, "W_outcome"),
        caller = "c2-crossfit-timing(gamma calibration)"
      )
    )
    alpha_plugin_block <- c(0, unlist(alpha_init_list, use.names = FALSE))

    gamma_final <- .c2_timing_stage(
      progress_path, rows_env, "source_fit_unified_density_ratio",
      sprintf("source=%s k1=%d n_cal=%d p_z=%d p_w_plugin=%d",
              s, k1, length(A_cal_stack), ncol(Z_cal_stack), ncol(W_plugin_block)),
      fit_unified_density_ratio(
        Z_site = Z_cal_stack,
        A = A_cal_stack,
        mean_grad_psi = mean_grad_psi_avg,
        alpha_init = alpha_plugin_block,
        lambda = NULL,
        calibrated = TRUE,
        M_tau = M_TAU_DEFAULT,
        W_outcome = W_plugin_block,
        A_val = A_val,
        family_int = family_int,
        link_int = link_int,
        lambda_rule = nuisance_rule
      ),
      extra_done = .c2_timing_lambda_detail
    )

    Z_plugin_block <- .c2_timing_stage(
      progress_path, rows_env, "source_make_Z_plugin_block",
      sprintf("source=%s k1=%d", s, k1),
      .make_plugin_block_design(
        lapply(source_calib_list, `[[`, "Z_site"),
        caller = "c2-crossfit-timing(alpha calibration)"
      )
    )
    gamma_plugin_block <- c(0, unlist(gamma_init_list, use.names = FALSE))

    alpha_final <- .c2_timing_stage(
      progress_path, rows_env, "source_fit_unified_outcome",
      sprintf("source=%s k1=%d n_cal=%d p_w=%d p_z_plugin=%d",
              s, k1, length(Y_cal_stack), ncol(W_cal_stack), ncol(Z_plugin_block)),
      fit_unified_outcome(
        W_outcome = W_cal_stack,
        Y = Y_cal_stack,
        A = A_cal_stack,
        A_val = A_val,
        gamma_s = gamma_plugin_block,
        lambda = NULL,
        family_int = family_int,
        link_int = link_int,
        calibrated = TRUE,
        M_tau = M_TAU_DEFAULT,
        Z_site = Z_plugin_block,
        lambda_rule = nuisance_rule
      ),
      extra_done = .c2_timing_lambda_detail
    )

    source_main_fold <- .c2_timing_stage(
      progress_path, rows_env, "source_materialize_main",
      sprintf("source=%s k1=%d", s, k1),
      materialize_fold(source_folds[[s]], k1)
    )
    correction <- .c2_timing_stage(
      progress_path, rows_env, "source_correction_term",
      sprintf("source=%s k1=%d n_main=%d", s, k1, source_main_fold$n),
      calculate_correction_term_cpp(
        source_main_fold$Z_site, source_main_fold$A, source_main_fold$Y,
        gamma_final, alpha_final, source_main_fold$W_outcome,
        M_TAU_INFERENCE_DEFAULT, family_int, link_int, A_val
      )
    )

    target_main_fold <- materialize_fold(target_folds, k1)
    mu_pred <- .c2_timing_stage(
      progress_path, rows_env, "target_predict_main",
      sprintf("source=%s k1=%d n_target_main=%d", s, k1, target_main_fold$n),
      mean(predict_glm_cpp(target_main_fold$W_outcome, alpha_final,
                           family_int, link_int))
    )
    mu_ts <- mu_pred + correction$delta_ts
    .c2_timing_log(
      progress_path, "source_done",
      sprintf("source=%s k1=%d mu_pred=%.6f delta=%.6f mu_ts=%.6f",
              s, k1, mu_pred, correction$delta_ts, mu_ts)
    )
  }

  timing_df <- do.call(rbind, rows_env$rows)
  write.csv(timing_df, timing_path, row.names = FALSE)
  cat(sprintf("\n[c2-crossfit-timing] Progress log: %s\n", progress_path))
  cat(sprintf("[c2-crossfit-timing] Stage timing CSV: %s\n", timing_path))
  cat("\n[c2-crossfit-timing] Slowest completed stages:\n")
  print(timing_df[order(-timing_df$elapsed_sec), , drop = FALSE][seq_len(min(12L, nrow(timing_df))), ],
        row.names = FALSE, digits = 4)

  expect_true(file.exists(progress_path))
  expect_true(file.exists(timing_path))
  expect_true(nrow(timing_df) > 0L)
})
