#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/verify_kf3.R
# ----------------------------------------------------------------------------
# Production-scale verification of the K_f=3 hypothesis.
#
# Reruns the EXACT setting that produced the 0.88 coverage in
#   results/n5000_K3_p10_C1_superpopulation_binary_model_tfmild_htnone_ss0.5_kf10_ateFALSE_summary.csv
# but with n_folds = 3 instead of 10, and writes the new summary side-by-side
# under diagnosis/bias/verify_out/ for direct comparison.
#
# This is the experiment that closes the loop: if one_round / two_round
# coverage moves from 0.886 / 0.878 (K_f=10) to ~0.94+ (K_f=3), the
# diagnosis is confirmed end-to-end. The algorithm itself is not modified.
#
# Args (positional, optional):
#   1: n_folds  (integer, default 3)
#   2: n_sims   (integer, default 200)
#
# Output:
#   diagnosis/bias/verify_out/verify_C1_n5000_K3_p10_kf<K_f>_summary.csv
#   diagnosis/bias/verify_out/verify_C1_n5000_K3_p10_kf<K_f>.csv  (per-sim)
#   diagnosis/bias/verify_out/verify_C1_n5000_K3_p10_kf<K_f>.RData
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
n_folds <- if (length(args) >= 1L) as.integer(args[[1]]) else 3L
n_sims  <- if (length(args) >= 2L) as.integer(args[[2]]) else 200L
config  <- if (length(args) >= 3L) args[[3]]            else "C1"

stopifnot(is.finite(n_folds), n_folds >= 3L,
          is.finite(n_sims),  n_sims  >= 1L,
          config %in% c("C1", "C2", "C3", "C4"))

# Pinned setting — must match results/n5000_K3_p10_<config>_*_summary.csv
# exactly except for n_folds.
n_total <- 5000L
K       <- 3L
p       <- 10L
shift   <- 0.5

cat(sprintf("[verify_kf3] n=%d K=%d p=%d config=%s shift=%g K_f=%d n_sims=%d\n",
            n_total, K, p, config, shift, n_folds, n_sims))

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) {
    library(RoCE)
  } else if (requireNamespace("devtools", quietly = FALSE)) {
    devtools::load_all(".")
  } else {
    stop("Neither RoCE nor devtools available.")
  }
})

out_dir <- "diagnosis/bias/verify_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setting_id <- sprintf("verify_%s_n%d_K%d_p%d_kf%d", config, n_total, K, p, n_folds)
csv_path     <- file.path(out_dir, sprintf("%s.csv", setting_id))
summary_path <- file.path(out_dir, sprintf("%s_summary.csv", setting_id))
rdata_path   <- file.path(out_dir, sprintf("%s.RData", setting_id))

# ---- checkpoint config (mirrors run_setting.R pattern) --------------------
checkpoint_dir <- Sys.getenv("CHECKPOINT_DIR",
                             sprintf("diagnosis/bias/checkpoints/verify_kf%d", n_folds))
job_id <- Sys.getenv("CHECKPOINT_JOB_ID", "")
if (!nzchar(job_id)) job_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
ckpt_config <- init_checkpoint_config(checkpoint_dir, setting_id, job_id)
ckpt_config$sim_interval <- 5L   # tight save cadence
cat(sprintf("[verify_kf3] checkpoint file: %s (sim_interval=%d)\n",
            ckpt_config$file, ckpt_config$sim_interval))

slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
n_cores <- if (nzchar(slurm_cpus)) {
  max(1L, as.integer(slurm_cpus))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
n_cores <- min(n_cores, n_sims)
cat(sprintf("[verify_kf3] using %d cores\n", n_cores))

t0 <- Sys.time()
results <- run_simulation_study(
  n_sims          = n_sims,
  n_total_vec     = n_total,
  K_vec           = K,
  p_vec           = p,
  configs         = config,
  n_cores         = n_cores,
  checkpoint_config = ckpt_config,
  nested_parallel = TRUE,
  estimand_type   = "superpopulation",
  site_allocation = "model",
  transform_type  = "mild",
  outcome_type    = "binary",
  heterogeneity_type = "none",
  shift_strength  = shift,
  n_folds         = n_folds,        # <-- the only knob being changed
  use_lambda_cache = TRUE,
  verbose_every   = 50L,
  parallel_strategy = "outer_priority",
  estimate_ate    = FALSE,
  dgp_type        = "roce"
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[verify_kf3] simulation elapsed = %.2f min\n", elapsed))

summary_df <- summarize_results(results)
summary_df$n_folds     <- n_folds
summary_df$elapsed_min <- round(elapsed, 3)

write.csv(results,    file = csv_path,     row.names = FALSE)
write.csv(summary_df, file = summary_path, row.names = FALSE)
save(results, summary_df, file = rdata_path)

cleanup_checkpoint(ckpt_config$file,
                   ckpt_config$preempt_signal,
                   ckpt_config$saved_signal)

cat(sprintf("\n[verify_kf3] wrote %s\n", summary_path))
cat("\n========== Coverage / bias by method (this run) ==========\n")
print(summary_df[, c("method", "bias_mean", "bias_sd", "se_mean",
                     "coverage", "n_success", "n_folds")],
      row.names = FALSE)

# Side-by-side print with original K_f=10 production summary (if found).
old_path <- sprintf(
  "results/n5000_K3_p10_%s_superpopulation_binary_model_tfmild_htnone_ss0.5_kf10_ateFALSE_summary.csv",
  config)
if (file.exists(old_path)) {
  cat("\n========== Side-by-side: K_f=10 (production) vs K_f=", n_folds, " (this run) ==========\n", sep = "")
  old <- read.csv(old_path, stringsAsFactors = FALSE)
  merged <- merge(
    old[,    c("method", "bias_mean", "coverage")],
    summary_df[, c("method", "bias_mean", "coverage")],
    by = "method", suffixes = c("_kf10", sprintf("_kf%d", n_folds))
  )
  merged$delta_coverage <- merged[[sprintf("coverage_kf%d", n_folds)]] - merged$coverage_kf10
  merged$bias_kf10 <- round(merged$bias_mean_kf10, 5)
  merged[[sprintf("bias_kf%d", n_folds)]] <-
      round(merged[[sprintf("bias_mean_kf%d", n_folds)]], 5)
  show <- merged[, c("method",
                     "bias_kf10", "coverage_kf10",
                     sprintf("bias_kf%d", n_folds),
                     sprintf("coverage_kf%d", n_folds),
                     "delta_coverage")]
  print(show, row.names = FALSE)
}

cat("\n[verify_kf3] DONE\n")
