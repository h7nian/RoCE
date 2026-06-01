#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/run_setting.R
# ----------------------------------------------------------------------------
# Runs ONE small simulation setting (used as the per-task script of the
# bias_diagnose SLURM array job). Wraps run_simulation_study() with a fixed
# C1 config and writes the summary to diagnosis/bias/results/.
#
# Purpose: diagnose where the positive bias of one_round_crossfit /
# two_round_crossfit comes from under correctly-specified C1.
# THE ALGORITHM IS NOT MODIFIED — this is a measurement-only suite.
#
# Args (all required, positional):
#   1: tag        short label embedded in the output filename
#   2: n_total    sample size
#   3: K          number of source sites
#   4: p          number of covariates
#   5: shift      shift_strength (numeric)
#   6: n_sims     number of Monte Carlo replicates
#
# Output:
#   diagnosis/bias/results/<setting_id>.csv
#   diagnosis/bias/results/<setting_id>_summary.csv
#   diagnosis/bias/results/<setting_id>.RData
#
# Checkpoint:
#   Mirrors main.R — uses init_checkpoint_config() with
#   CHECKPOINT_DIR (env, default diagnosis/bias/checkpoints) and the
#   job id provided via CHECKPOINT_JOB_ID (env). The simulation-level
#   resume is handled by run_simulation_study(); on a SIGUSR1 the
#   shell handler in bias_diagnose.cmd touches the preempt signal, the
#   R loop flushes state, and SLURM requeues the same array task.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(sprintf(
    "run_setting.R: expected 6 args (tag,n,K,p,shift,n_sims), got %d.",
    length(args)
  ))
}

tag         <- args[[1]]
n_total     <- as.integer(args[[2]])
K           <- as.integer(args[[3]])
p           <- as.integer(args[[4]])
shift       <- as.numeric(args[[5]])
n_sims      <- as.integer(args[[6]])

stopifnot(is.finite(n_total), n_total > 0,
          is.finite(K), K >= 1,
          is.finite(p), p >= 1,
          is.finite(shift), shift >= 0,
          is.finite(n_sims), n_sims >= 1)

cat(sprintf("[diag/bias] tag=%s  n=%d  K=%d  p=%d  shift=%g  sims=%d\n",
            tag, n_total, K, p, shift, n_sims))

# ---- Load package -----------------------------------------------------------
if (requireNamespace("FACEHD", quietly = FALSE)) {
  library(FACEHD)
} else if (requireNamespace("devtools", quietly = FALSE)) {
  devtools::load_all(".")
} else {
  stop("Neither installed FACEHD nor devtools available.")
}

# ---- Output dirs ------------------------------------------------------------
out_dir <- "diagnosis/bias/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setting_id <- sprintf("%s_n%d_K%d_p%d_ss%g", tag, n_total, K, p, shift)
csv_path     <- file.path(out_dir, sprintf("%s.csv", setting_id))
summary_path <- file.path(out_dir, sprintf("%s_summary.csv", setting_id))
rdata_path   <- file.path(out_dir, sprintf("%s.RData", setting_id))

# ---- Checkpoint config (mirrors main.R) -------------------------------------
# CHECKPOINT_DIR + CHECKPOINT_JOB_ID are exported by bias_diagnose.cmd; we
# fall back to sensible defaults so the script also runs interactively.
checkpoint_dir <- Sys.getenv("CHECKPOINT_DIR", "diagnosis/bias/checkpoints")
job_id <- Sys.getenv("CHECKPOINT_JOB_ID", "")
if (!nzchar(job_id)) {
  job_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
}
ckpt_config <- init_checkpoint_config(checkpoint_dir, setting_id, job_id)
# Tighten the in-R simulation-level checkpoint interval for the diagnostic
# suite. The package default (20) is fine when a setting fits within one
# walltime; here individual sims are slow and walltime hits before the
# default threshold, which strands progress. Saving every 5 sims keeps
# overhead negligible while ensuring requeue can resume meaningfully.
# (Not an algorithm change — only the persistence cadence of the loop.)
ckpt_config$sim_interval <- 5L
cat(sprintf("[diag/bias] checkpoint file: %s (sim_interval=%d)\n",
            ckpt_config$file, ckpt_config$sim_interval))

# ---- Core count -------------------------------------------------------------
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
n_cores <- if (nzchar(slurm_cpus)) {
  max(1L, as.integer(slurm_cpus))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
n_cores <- min(n_cores, n_sims)
cat(sprintf("[diag/bias] using %d cores (n_sims=%d)\n", n_cores, n_sims))

# ---- Run --------------------------------------------------------------------
t0 <- Sys.time()
results <- run_simulation_study(
  n_sims          = n_sims,
  n_total_vec     = n_total,
  K_vec           = K,
  p_vec           = p,
  configs         = "C1",
  n_cores         = n_cores,
  checkpoint_config = ckpt_config,
  nested_parallel = TRUE,
  estimand_type   = "superpopulation",
  site_allocation = "model",
  transform_type  = "mild",
  outcome_type    = "binary",
  heterogeneity_type = "none",
  shift_strength  = shift,
  n_folds         = 10L,
  use_lambda_cache = TRUE,
  verbose_every   = 25L,
  parallel_strategy = "outer_priority",
  estimate_ate    = FALSE,
  dgp_type        = "facehd"
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[diag/bias] simulation elapsed = %.2f min\n", elapsed))

# ---- Summarize + persist ----------------------------------------------------
summary_df <- summarize_results(results)
summary_df$tag         <- tag
summary_df$shift       <- shift
summary_df$elapsed_min <- round(elapsed, 3)

write.csv(results,    file = csv_path,     row.names = FALSE)
write.csv(summary_df, file = summary_path, row.names = FALSE)
save(results, summary_df, file = rdata_path)

cleanup_checkpoint(ckpt_config$file,
                   ckpt_config$preempt_signal,
                   ckpt_config$saved_signal)

cat(sprintf("[diag/bias] wrote %s\n", summary_path))
cat("[diag/bias] DONE\n")
