#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/probe_scan.R
# ----------------------------------------------------------------------------
# Quantifies the "calib fold too small / p too big" hypothesis identified in
# probe.R Phase B: mu_pred_ts swings wildly across k2 folds when the
# (k1,k2) calibration fold has too few obs to identify a p-dim nuisance.
#
# Scans across (n_total, p) at K=3, fixed config C1. For each (n, p, seed,
# source), it runs Phase-B-style probing (k1=1, k2 in secondary) and
# records per-fold magnitudes plus aggregate spread statistics.
#
# THE ALGORITHM IS NOT MODIFIED. Only observation.
#
# Output: diagnosis/bias/probe_out/probe_scan_<job>.csv
#   one row per (n, p, K, seed, source, k1, k2) with intermediate magnitudes,
#   plus a wide "spread" CSV summarizing across k2.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
job_tag <- if (length(args) >= 1L) args[[1]] else format(Sys.time(), "%Y%m%d_%H%M%S")
n_seeds <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L

cat(sprintf("[probe_scan] job_tag=%s  n_seeds=%d\n", job_tag, n_seeds))

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else if (requireNamespace("RoCE", quietly = TRUE)) {
    library(RoCE)
  } else {
    stop("Neither devtools nor RoCE available.")
  }
})

out_dir <- "diagnosis/bias/probe_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Scan grid -------------------------------------------------------------
# (n_total, p, K). K=3 throughout. n_seeds replicates per setting.
SCAN_GRID <- data.frame(
  tag = c("n1k_p10",  "n1k_p20",  "n1k_p50",
          "n2k_p20",  "n5k_p10",  "n5k_p20"),
  n   = c(1000,       1000,       1000,
          2000,       5000,       5000),
  p   = c(  10,         20,         50,
            20,         10,         20),
  K   = 3,
  stringsAsFactors = FALSE
)
cat("[probe_scan] scan grid:\n"); print(SCAN_GRID, row.names = FALSE)

# ---- Per-(setting, seed, source, k1) probe routine ------------------------
# Returns list of per-k2 records.
probe_one <- function(setting, seed) {
  set.seed(seed)
  data <- generate_simulation_data(
    n_total = setting$n, K = setting$K, p = setting$p, config = "C1",
    estimand_type = "superpopulation", site_allocation = "model",
    transform_type = "mild", outcome_type = "binary",
    heterogeneity_type = "none", shift_strength = 0.5,
    dgp_type = "roce"
  )
  data_split <- split_data_by_site(data)
  source_sites <- setdiff(names(data_split), "t")
  target_data  <- data_split[["t"]]

  n_folds <- 10L
  fold_seed_base <- target_data$n + sum(sapply(source_sites,
                                               function(s) data_split[[s]]$n))
  target_folds <- partition_into_folds(target_data, n_folds, seed = fold_seed_base)
  source_folds <- setNames(
    lapply(seq_along(source_sites), function(i) {
      partition_into_folds(data_split[[source_sites[i]]], n_folds,
                           seed = fold_seed_base + i)
    }),
    source_sites
  )

  glm_spec <- resolve_glm_family("binomial")
  family_int <- glm_spec$family_int
  link_int   <- glm_spec$link_int
  A_val <- 1L
  M_tau <- 10

  records <- list()
  k1 <- 1L
  secondary <- setdiff(1:n_folds, k1)

  for (s in source_sites) {
    for (k2 in secondary) {
      training_folds <- setdiff(1:n_folds, c(k1, k2))
      source_train <- combine_folds(source_folds[[s]], training_folds)
      source_calib <- materialize_fold(source_folds[[s]], k2)
      target_calib <- materialize_fold(target_folds, k2)
      target_main  <- materialize_fold(target_folds, k1)

      n_calib_arm <- sum(source_calib$A == A_val)
      n_train_arm <- sum(source_train$A == A_val)

      alpha_init <- tryCatch(
        fit_initial_outcome(source_train$W_outcome, source_train$Y,
                            source_train$A, A_val = A_val, family = "binomial"),
        error = function(e) NULL)
      if (is.null(alpha_init)) next
      alpha_init_abs_max <- max(abs(alpha_init))

      mean_phi <- c(1, colMeans(target_calib$Z_site))
      mean_grad_psi <- tryCatch(
        .mean_glm_gradient_site_basis(target_calib$W_outcome, target_calib$Z_site,
                                      alpha_init, family_int, link_int),
        error = function(e) NULL)
      if (is.null(mean_grad_psi)) next

      gamma_init <- tryCatch(
        fit_initial_density_ratio(source_train$Z_site, source_train$A,
                                  mean_phi, A_val = A_val),
        error = function(e) NULL)
      if (is.null(gamma_init)) next
      gamma_init_abs_max <- max(abs(gamma_init))

      Zc <- cbind(1, source_calib$Z_site)
      w_lin <- as.numeric(Zc %*% gamma_init)
      w_exp_max <- max(exp(w_lin))

      gamma_cal <- tryCatch(
        fit_unified_density_ratio(
          Z_site = source_calib$Z_site, A = source_calib$A,
          mean_grad_psi = mean_grad_psi, alpha_init = alpha_init,
          calibrated = TRUE, M_tau = M_tau,
          W_outcome = source_calib$W_outcome, A_val = A_val,
          family_int = family_int, link_int = link_int),
        error = function(e) NULL)
      gamma_cal_abs_max <- if (is.null(gamma_cal)) NA_real_ else max(abs(gamma_cal))

      alpha_cal <- tryCatch(
        fit_unified_outcome(
          W_outcome = source_calib$W_outcome, Y = source_calib$Y,
          A = source_calib$A, A_val = A_val, gamma_s = gamma_init,
          family_int = family_int, link_int = link_int,
          calibrated = TRUE, M_tau = M_tau, Z_site = source_calib$Z_site),
        error = function(e) NULL)
      alpha_cal_abs_max <- if (is.null(alpha_cal)) NA_real_ else max(abs(alpha_cal))

      mu_pred_ts <- if (!is.null(alpha_cal)) {
        tryCatch(mean(predict_glm_cpp(target_main$W_outcome, alpha_cal,
                                      family_int, link_int)),
                 error = function(e) NA_real_)
      } else NA_real_

      records[[length(records) + 1L]] <- list(
        source = s, k1 = k1, k2 = k2,
        n_train = source_train$n, n_train_arm = n_train_arm,
        n_calib = source_calib$n, n_calib_arm = n_calib_arm,
        n_target_calib = target_calib$n,
        alpha_init_abs_max = alpha_init_abs_max,
        gamma_init_abs_max = gamma_init_abs_max,
        w_exp_max = w_exp_max,
        gamma_cal_abs_max = gamma_cal_abs_max,
        alpha_cal_abs_max = alpha_cal_abs_max,
        mu_pred_ts = mu_pred_ts
      )
    }
  }
  records
}

# ---- Main loop -------------------------------------------------------------
all_rows <- list()
spread_rows <- list()
t0_all <- Sys.time()

for (ig in seq_len(nrow(SCAN_GRID))) {
  setting <- SCAN_GRID[ig, ]
  cat(sprintf("\n[probe_scan] setting %d/%d  tag=%s  n=%d K=%d p=%d\n",
              ig, nrow(SCAN_GRID), setting$tag, setting$n, setting$K, setting$p))
  for (seed in seq_len(n_seeds) + 2000L) {
    t0 <- Sys.time()
    recs <- probe_one(setting, seed)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    cat(sprintf("  seed=%d  records=%d  elapsed=%.1fs\n",
                seed, length(recs), elapsed))
    if (length(recs) == 0L) next

    # Append long rows
    for (r in recs) {
      all_rows[[length(all_rows) + 1L]] <- data.frame(
        tag = setting$tag, n = setting$n, K = setting$K, p = setting$p,
        seed = seed,
        source = r$source, k1 = r$k1, k2 = r$k2,
        n_train = r$n_train, n_train_arm = r$n_train_arm,
        n_calib = r$n_calib, n_calib_arm = r$n_calib_arm,
        n_target_calib = r$n_target_calib,
        alpha_init_abs_max = r$alpha_init_abs_max,
        gamma_init_abs_max = r$gamma_init_abs_max,
        w_exp_max = r$w_exp_max,
        gamma_cal_abs_max = r$gamma_cal_abs_max,
        alpha_cal_abs_max = r$alpha_cal_abs_max,
        mu_pred_ts = r$mu_pred_ts,
        stringsAsFactors = FALSE
      )
    }

    # Per-source spread summary (across k2)
    by_source <- split(recs, vapply(recs, function(r) r$source, character(1)))
    for (s in names(by_source)) {
      mu <- vapply(by_source[[s]], function(r) r$mu_pred_ts, numeric(1))
      gc <- vapply(by_source[[s]], function(r) r$gamma_cal_abs_max, numeric(1))
      we <- vapply(by_source[[s]], function(r) r$w_exp_max, numeric(1))
      spread_rows[[length(spread_rows) + 1L]] <- data.frame(
        tag = setting$tag, n = setting$n, K = setting$K, p = setting$p,
        seed = seed, source = s,
        n_calib_arm = by_source[[s]][[1]]$n_calib_arm,
        ratio_p_over_calib_arm = setting$p / by_source[[s]][[1]]$n_calib_arm,
        mu_pred_ts_min = min(mu, na.rm = TRUE),
        mu_pred_ts_max = max(mu, na.rm = TRUE),
        mu_pred_ts_spread = max(mu, na.rm = TRUE) - min(mu, na.rm = TRUE),
        mu_pred_ts_sd = sd(mu, na.rm = TRUE),
        gamma_cal_abs_max = max(gc, na.rm = TRUE),
        w_exp_max = max(we, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
}

elapsed_all <- as.numeric(difftime(Sys.time(), t0_all, units = "mins"))
cat(sprintf("\n[probe_scan] total elapsed: %.1f min\n", elapsed_all))

long_df   <- do.call(rbind, all_rows)
spread_df <- do.call(rbind, spread_rows)

long_path   <- file.path(out_dir, sprintf("probe_scan_long_%s.csv",   job_tag))
spread_path <- file.path(out_dir, sprintf("probe_scan_spread_%s.csv", job_tag))
write.csv(long_df,   file = long_path,   row.names = FALSE)
write.csv(spread_df, file = spread_path, row.names = FALSE)
cat(sprintf("[probe_scan] wrote %s (%d rows)\n", long_path, nrow(long_df)))
cat(sprintf("[probe_scan] wrote %s (%d rows)\n", spread_path, nrow(spread_df)))

# ---- Console summary -------------------------------------------------------
cat("\n========================  PER-SETTING SPREAD MEDIANS  ========================\n")
agg <- aggregate(
  cbind(mu_pred_ts_spread, mu_pred_ts_sd, gamma_cal_abs_max, w_exp_max,
        n_calib_arm, ratio_p_over_calib_arm) ~ tag + n + p + K,
  data = spread_df,
  FUN = function(x) round(median(x, na.rm = TRUE), 4)
)
agg <- agg[order(agg$n, agg$p), ]
print(agg, row.names = FALSE)

cat("\n========================  EXPLOSION RATE (mu_pred_ts_spread > 0.2)  ========================\n")
spread_df$is_blown <- spread_df$mu_pred_ts_spread > 0.2
blown_rate <- aggregate(is_blown ~ tag + n + p, data = spread_df,
                        FUN = function(x) mean(x, na.rm = TRUE))
blown_rate <- blown_rate[order(blown_rate$n, blown_rate$p), ]
names(blown_rate)[ncol(blown_rate)] <- "frac_seeds_sources_blown"
print(blown_rate, row.names = FALSE)
