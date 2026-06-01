.c2_dr_cv_timing_enabled <- function() {
  env_enabled <- Sys.getenv("FACEHD_RUN_C2_DR_CV_TIMING", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-dr-cv-timing", filter, fixed = TRUE)
}

.c2_dr_cv_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_dr_cv_num_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || !is.finite(parsed)) default else parsed
}

.c2_dr_cv_int_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  parsed <- parsed[is.finite(parsed) & parsed > 0L]
  if (length(parsed) == 0L) default else parsed
}

.c2_dr_cv_source_env <- function(available) {
  value <- Sys.getenv("FACEHD_C2_DR_CV_TIMING_SOURCE", unset = "")
  if (!nzchar(value)) return(available[1L])
  if (!(value %in% available)) {
    stop(sprintf("FACEHD_C2_DR_CV_TIMING_SOURCE=%s unavailable; available: %s",
                 value, paste(available, collapse = ", ")), call. = FALSE)
  }
  value
}

.c2_dr_cv_log <- function(path, state, detail) {
  line <- sprintf("[%s] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                  state, detail)
  cat(line, "\n", sep = "")
  flush.console()
  cat(line, "\n", file = path, append = TRUE, sep = "")
}

test_that("c2-dr-cv-timing times initial density-ratio lambda grids", {
  skip_if_not(.c2_dr_cv_timing_enabled(),
              message = "set FACEHD_RUN_C2_DR_CV_TIMING=1 or run with --filter c2-dr-cv-timing")

  n_total <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_N", 5000L)
  K_sites <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_K", 5L)
  p <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_P", 50L)
  n_folds <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_FOLDS", 3L)
  seed <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_SEED", 270001L)
  k1 <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_K1", 1L)
  k2 <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_K2", 2L)
  max_iter <- .c2_dr_cv_int_env("FACEHD_C2_DR_CV_TIMING_MAX_ITER", MAX_ITER_DEFAULT)
  tol <- .c2_dr_cv_num_env("FACEHD_C2_DR_CV_TIMING_TOL", TOL_DEFAULT)
  shift_strength <- .c2_dr_cv_num_env("FACEHD_C2_DR_CV_TIMING_SHIFT", 0.5)
  grid_sizes <- .c2_dr_cv_int_vector_env("FACEHD_C2_DR_CV_TIMING_GRIDS", c(5L, 10L, 20L))

  out_dir <- file.path("c2_dr_cv_timing_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  progress_path <- file.path(out_dir, paste0(ts, "_c2_dr_cv_timing_progress.log"))
  csv_path <- file.path(out_dir, paste0(ts, "_c2_dr_cv_timing.csv"))

  .c2_dr_cv_log(
    progress_path, "config",
    sprintf("n_total=%d K=%d p=%d folds=%d seed=%d k1=%d k2=%d max_iter=%d grids=%s",
            n_total, K_sites, p, n_folds, seed, k1, k2, max_iter,
            paste(grid_sizes, collapse = ","))
  )

  set.seed(seed)
  data <- generate_simulation_data(
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
    dgp_type = "facehd",
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  source_sites <- setdiff(names(data_split), "t")
  source <- .c2_dr_cv_source_env(source_sites)

  folds <- build_crossfit_folds(data_split, n_folds)
  training_folds <- setdiff(seq_len(n_folds), c(k1, k2))
  target_train <- combine_folds(folds$target_folds, training_folds)
  source_train <- combine_folds(folds$source_folds[[source]], training_folds)
  mean_phi <- c(1, colMeans(target_train$Z_site))
  n_cv_folds <- .nuisance_cv_fold_count(source_train$A, 1L, "c2-dr-cv-timing")
  lambda_max <- compute_lambda_max_initial_dr(source_train$Z_site, source_train$A,
                                              mean_phi, A_val = 1L)
  lambda_min_ratio <- if (nrow(source_train$Z_site) > ncol(source_train$Z_site)) {
    LAMBDA_MIN_RATIO_LOW_DIM
  } else {
    LAMBDA_MIN_RATIO_HIGH_DIM
  }

  .c2_dr_cv_log(
    progress_path, "data",
    sprintf("source=%s train_folds=%s n_source_train=%d n_arm=%d p_z=%d n_cv_folds=%d lambda_max=%.6g",
            source, paste(training_folds, collapse = ","), source_train$n,
            sum(source_train$A == 1L), ncol(source_train$Z_site), n_cv_folds,
            lambda_max)
  )

  rows <- list()
  for (nlambda in grid_sizes) {
    lambda_grid <- build_lambda_grid(lambda_max = lambda_max,
                                     lambda_min_ratio = lambda_min_ratio,
                                     nlambda = nlambda)
    .c2_dr_cv_log(progress_path, "cv_start",
                  sprintf("source=%s nlambda=%d", source, nlambda))
    start <- Sys.time()
    cv_result <- select_lambda_cv_initial_density_ratio_cpp(
      source_train$Z_site, source_train$A, mean_phi, lambda_grid,
      n_cv_folds, max_iter, tol, 1L
    )
    cv_elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    .c2_dr_cv_log(
      progress_path, "cv_done",
      sprintf("source=%s nlambda=%d elapsed_sec=%.3f lambda_min=%.6g lambda_1se=%.6g idx_min=%d idx_1se=%d",
              source, nlambda, cv_elapsed, cv_result$lambda_min,
              cv_result$lambda_1se, cv_result$idx_min, cv_result$idx_1se)
    )

    fit_start <- Sys.time()
    gamma <- fit_initial_density_ratio(
      source_train$Z_site, source_train$A, mean_phi,
      lambda = cv_result$lambda_min,
      max_iter = max_iter,
      tol = tol,
      A_val = 1L
    )
    fit_elapsed <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))
    gamma_support <- if (length(gamma) <= 1L) 0L else sum(abs(as.numeric(gamma[-1L])) > 1e-8)
    rows[[length(rows) + 1L]] <- data.frame(
      source = source,
      nlambda = nlambda,
      cv_elapsed_sec = cv_elapsed,
      fit_elapsed_sec = fit_elapsed,
      lambda_min = as.numeric(cv_result$lambda_min),
      lambda_1se = as.numeric(cv_result$lambda_1se),
      idx_min = as.integer(cv_result$idx_min),
      idx_1se = as.integer(cv_result$idx_1se),
      gamma_support = gamma_support,
      gamma_l2 = sqrt(mean(as.numeric(gamma)^2)),
      n_source_train = source_train$n,
      n_arm = sum(source_train$A == 1L),
      p_z = ncol(source_train$Z_site),
      n_cv_folds = n_cv_folds,
      max_iter = max_iter,
      tol = tol,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, rows)
  write.csv(out, csv_path, row.names = FALSE)
  cat(sprintf("\n[c2-dr-cv-timing] Progress log: %s\n", progress_path))
  cat(sprintf("[c2-dr-cv-timing] CSV: %s\n", csv_path))
  print(out, row.names = FALSE, digits = 4)

  expect_true(file.exists(progress_path))
  expect_true(file.exists(csv_path))
  expect_true(all(is.finite(out$cv_elapsed_sec)))
  expect_true(all(is.finite(out$lambda_min)))
})
