#!/usr/bin/env Rscript
# Summarize the rho=0 n_folds coverage pilot: bias, empirical SD, RMSE, mean SE,
# coverage and CI width per (config, K, n_folds, method).
#
# usage: summarize_pilot_nfolds.R [PILOT_RAW_DIR]
args <- commandArgs(trailingOnly = TRUE)
raw_dir <- if (length(args) >= 1L) args[[1L]] else
  "results/direct_tate_mc500_b5000/pilot_nfolds_20260917/raw"
files <- list.files(raw_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0L) stop("no pilot CSVs under ", raw_dir, call. = FALSE)

rows <- do.call(rbind, lapply(files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE)
  d[d$method %in% c("one_round_crossfit_ate", "target_only_ate"),
    c("sim_id", "config", "K", "pilot_n_folds", "method", "estimate", "se",
      "bias", "coverage", "ci_width", "pilot_elapsed_seconds")]
}))

summary_rows <- do.call(rbind, lapply(
  split(rows, list(rows$config, rows$K, rows$pilot_n_folds, rows$method), drop = TRUE),
  function(g) {
    n <- nrow(g)
    bias <- mean(g$bias)
    emp_sd <- stats::sd(g$estimate)
    data.frame(
      config = g$config[[1L]], K = g$K[[1L]], n_folds = g$pilot_n_folds[[1L]],
      method = g$method[[1L]], n_reps = n,
      bias = bias, bias_mcse = emp_sd / sqrt(n),
      emp_sd = emp_sd, rmse = sqrt(mean(g$bias^2)),
      mean_se = mean(g$se), se_ratio = mean(g$se) / emp_sd,
      abs_bias_over_sd = abs(bias) / emp_sd,
      coverage = mean(g$coverage),
      coverage_mcse = sqrt(mean(g$coverage) * (1 - mean(g$coverage)) / n),
      ci_width = mean(g$ci_width),
      mean_elapsed_min = mean(g$pilot_elapsed_seconds) / 60,
      stringsAsFactors = FALSE
    )
  }
))
summary_rows <- summary_rows[order(summary_rows$config, summary_rows$K,
                                   summary_rows$n_folds, summary_rows$method), ]
rownames(summary_rows) <- NULL
numeric_cols <- vapply(summary_rows, is.numeric, logical(1L))
printed <- summary_rows
printed[numeric_cols] <- lapply(printed[numeric_cols], signif, digits = 4)
print(printed, row.names = FALSE)
out <- file.path(dirname(raw_dir), "pilot_nfolds_summary.csv")
write.csv(summary_rows, out, row.names = FALSE)
cat("\nwritten:", out, "\n")
