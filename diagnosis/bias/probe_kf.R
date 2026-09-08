#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/probe_kf.R
# ----------------------------------------------------------------------------
# Sweeps K_f (cross-fitting fold count) at fixed (n_total, K, p) to test the
# hypothesis that LOWERING K_f (= LARGER calibration folds) reduces the
# per-fold mu_pred_ts spread that drives the one_round / two_round
# instability.
#
# K_f is a utility / tuning parameter of run_crossfit; reducing it does NOT
# modify the algorithm — it just changes the train/calib fold split sizes.
#
# Runs two settings:
#   (a) "catastrophic" : n=1000, K=3, p=20   — where bias was ±100
#   (b) "production"   : n=5000, K=3, p=10   — where coverage was 0.88
# For each setting, scans K_f ∈ {3, 5, 7, 10}, n_seeds per cell.
#
# Output:
#   diagnosis/bias/probe_out/probe_kf_long_<job>.csv
#   diagnosis/bias/probe_out/probe_kf_spread_<job>.csv
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
job_tag <- if (length(args) >= 1L) args[[1]] else format(Sys.time(), "%Y%m%d_%H%M%S")
n_seeds <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L

cat(sprintf("[probe_kf] job_tag=%s  n_seeds=%d\n", job_tag, n_seeds))

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

# ---- Scan grid: (setting, K_f) --------------------------------------------
SETTINGS <- data.frame(
  setting_tag = c("catastrophic_n1k_p20", "production_n5k_p10"),
  n_total     = c(1000,                   5000),
  K           = c(3,                      3),
  p           = c(20,                     10),
  stringsAsFactors = FALSE
)
KF_VALUES <- c(3L, 5L, 7L, 10L)
cat("[probe_kf] settings:\n"); print(SETTINGS, row.names = FALSE)
cat(sprintf("[probe_kf] K_f values: %s\n", paste(KF_VALUES, collapse = ",")))

# ---- Per-(setting, K_f, seed, source) probe routine -----------------------
probe_one <- function(setting, K_f, seed) {
  set.seed(seed)
  data <- generate_simulation_data(
    n_total = setting$n_total, K = setting$K, p = setting$p, config = "C1",
    estimand_type = "superpopulation", site_allocation = "model",
    transform_type = "mild", outcome_type = "binary",
    heterogeneity_type = "none", shift_strength = 0.5,
    dgp_type = "roce"
  )
  data_split   <- split_data_by_site(data)
  source_sites <- setdiff(names(data_split), "t")
  target_data  <- data_split[["t"]]

  fold_seed_base <- target_data$n + sum(sapply(source_sites,
                                               function(s) data_split[[s]]$n))
  target_folds <- partition_into_folds(target_data, K_f, seed = fold_seed_base)
  source_folds <- setNames(
    lapply(seq_along(source_sites), function(i) {
      partition_into_folds(data_split[[source_sites[i]]], K_f,
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
  secondary <- setdiff(seq_len(K_f), k1)

  for (s in source_sites) {
    for (k2 in secondary) {
      training_folds <- setdiff(seq_len(K_f), c(k1, k2))
      source_train <- combine_folds(source_folds[[s]], training_folds)
      source_calib <- materialize_fold(source_folds[[s]], k2)
      target_calib <- materialize_fold(target_folds, k2)
      target_main  <- materialize_fold(target_folds, k1)

      alpha_init <- tryCatch(
        fit_initial_outcome(source_train$W_outcome, source_train$Y,
                            source_train$A, A_val = A_val, family = "binomial"),
        error = function(e) NULL)
      if (is.null(alpha_init)) next

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

      Zc <- cbind(1, source_calib$Z_site)
      w_exp_max <- max(exp(as.numeric(Zc %*% gamma_init)))

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

      mu_pred_ts <- if (!is.null(alpha_cal)) {
        tryCatch(mean(predict_glm_cpp(target_main$W_outcome, alpha_cal,
                                      family_int, link_int)),
                 error = function(e) NA_real_)
      } else NA_real_

      records[[length(records) + 1L]] <- list(
        source = s, k1 = k1, k2 = k2,
        n_calib = source_calib$n,
        n_calib_arm = sum(source_calib$A == A_val),
        gamma_cal_abs_max = gamma_cal_abs_max,
        w_exp_max = w_exp_max,
        mu_pred_ts = mu_pred_ts
      )
    }
  }
  records
}

# ---- Main loop -------------------------------------------------------------
all_rows    <- list()
spread_rows <- list()
t0_all <- Sys.time()

for (is in seq_len(nrow(SETTINGS))) {
  setting <- SETTINGS[is, ]
  for (K_f in KF_VALUES) {
    cat(sprintf("\n[probe_kf] setting=%s K_f=%d\n", setting$setting_tag, K_f))
    for (seed in seq_len(n_seeds) + 3000L) {
      t0 <- Sys.time()
      recs <- probe_one(setting, K_f, seed)
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  K_f=%d seed=%d  records=%d  elapsed=%.1fs\n",
                  K_f, seed, length(recs), el))
      if (length(recs) == 0L) next

      for (r in recs) {
        all_rows[[length(all_rows) + 1L]] <- data.frame(
          setting_tag = setting$setting_tag,
          n_total = setting$n_total, K = setting$K, p = setting$p,
          K_f = K_f, seed = seed,
          source = r$source, k1 = r$k1, k2 = r$k2,
          n_calib = r$n_calib, n_calib_arm = r$n_calib_arm,
          gamma_cal_abs_max = r$gamma_cal_abs_max,
          w_exp_max = r$w_exp_max,
          mu_pred_ts = r$mu_pred_ts,
          stringsAsFactors = FALSE
        )
      }

      by_source <- split(recs, vapply(recs, function(r) r$source, character(1)))
      for (s in names(by_source)) {
        mu <- vapply(by_source[[s]], function(r) r$mu_pred_ts, numeric(1))
        gc <- vapply(by_source[[s]], function(r) r$gamma_cal_abs_max, numeric(1))
        we <- vapply(by_source[[s]], function(r) r$w_exp_max, numeric(1))
        spread_rows[[length(spread_rows) + 1L]] <- data.frame(
          setting_tag = setting$setting_tag,
          n_total = setting$n_total, K = setting$K, p = setting$p,
          K_f = K_f, seed = seed, source = s,
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
}

el_all <- as.numeric(difftime(Sys.time(), t0_all, units = "mins"))
cat(sprintf("\n[probe_kf] total elapsed: %.1f min\n", el_all))

long_df   <- do.call(rbind, all_rows)
spread_df <- do.call(rbind, spread_rows)

long_path   <- file.path(out_dir, sprintf("probe_kf_long_%s.csv",   job_tag))
spread_path <- file.path(out_dir, sprintf("probe_kf_spread_%s.csv", job_tag))
write.csv(long_df,   file = long_path,   row.names = FALSE)
write.csv(spread_df, file = spread_path, row.names = FALSE)
cat(sprintf("[probe_kf] wrote %s (%d rows)\n", long_path, nrow(long_df)))
cat(sprintf("[probe_kf] wrote %s (%d rows)\n", spread_path, nrow(spread_df)))

# ---- Console summary -------------------------------------------------------
cat("\n========================  MEDIAN SPREAD BY (setting, K_f)  ========================\n")
agg <- aggregate(
  cbind(mu_pred_ts_spread, mu_pred_ts_sd, gamma_cal_abs_max, w_exp_max,
        n_calib_arm, ratio_p_over_calib_arm) ~ setting_tag + K_f,
  data = spread_df,
  FUN = function(x) round(median(x, na.rm = TRUE), 4)
)
agg <- agg[order(agg$setting_tag, agg$K_f), ]
print(agg, row.names = FALSE)

cat("\n========================  EXPLOSION RATE (spread > 0.2) BY (setting, K_f)  ========================\n")
spread_df$is_blown <- spread_df$mu_pred_ts_spread > 0.2
rate <- aggregate(is_blown ~ setting_tag + K_f, data = spread_df,
                  FUN = function(x) round(mean(x, na.rm = TRUE), 3))
rate <- rate[order(rate$setting_tag, rate$K_f), ]
names(rate)[ncol(rate)] <- "frac_blown"
print(rate, row.names = FALSE)
