#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/c2/run_setting.R
# ----------------------------------------------------------------------------
# Runs one C2 diagnosis setting from the SLURM array.  This is measurement-only:
# it calls run_simulation_study(), writes raw/summary CSVs, and leaves estimator
# code untouched.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 12L) {
  stop(sprintf(
    "run_setting.R: expected 12 args (tag,n,K,p,config,estimand,allocation,transform,heterogeneity,shift,n_folds,n_sims), got %d.",
    length(args)
  ))
}

tag          <- args[[1]]
n_total      <- as.integer(args[[2]])
K            <- as.integer(args[[3]])
p            <- as.integer(args[[4]])
config       <- args[[5]]
estimand     <- args[[6]]
allocation   <- args[[7]]
transform    <- args[[8]]
heterogeneity <- args[[9]]
shift        <- as.numeric(args[[10]])
n_folds      <- as.integer(args[[11]])
n_sims       <- as.integer(args[[12]])

stopifnot(is.finite(n_total), n_total > 0,
          is.finite(K), K >= 1,
          is.finite(p), p >= 1,
          config %in% c("C1", "C2", "C3", "C4"),
          estimand %in% c("superpopulation", "sample"),
          allocation %in% c("model", "uniform", "balanced", "target_heavy",
                            "source_heavy", "independent"),
          transform %in% c("strong", "mild", "none"),
          heterogeneity %in% c("none", "mild", "strong", "partial"),
          is.finite(shift), shift > 0,
          is.finite(n_folds), n_folds >= 3L,
          is.finite(n_sims), n_sims >= 1L)

cat(sprintf(
  "[diag/c2] tag=%s n=%d K=%d p=%d config=%s estimand=%s allocation=%s transform=%s heterogeneity=%s shift=%g K_f=%d sims=%d\n",
  tag, n_total, K, p, config, estimand, allocation, transform, heterogeneity,
  shift, n_folds, n_sims
))

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) {
    library(RoCE)
  } else if (requireNamespace("devtools", quietly = FALSE)) {
    devtools::load_all(".")
  } else {
    stop("Neither installed RoCE nor devtools available.")
  }
})

out_dir <- "diagnosis/c2/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setting_id <- sprintf(
  "%s_n%d_K%d_p%d_%s_%s_%s_tf%s_ht%s_ss%g_kf%d",
  tag, n_total, K, p, config, estimand, allocation, transform,
  heterogeneity, shift, n_folds
)
csv_path     <- file.path(out_dir, sprintf("%s.csv", setting_id))
summary_path <- file.path(out_dir, sprintf("%s_summary.csv", setting_id))
rdata_path   <- file.path(out_dir, sprintf("%s.RData", setting_id))

checkpoint_dir <- Sys.getenv("CHECKPOINT_DIR", "diagnosis/c2/checkpoints")
job_id <- Sys.getenv("CHECKPOINT_JOB_ID", "")
if (!nzchar(job_id)) job_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
ckpt_config <- init_checkpoint_config(checkpoint_dir, setting_id, job_id)
ckpt_config$sim_interval <- 5L
cat(sprintf("[diag/c2] checkpoint file: %s (sim_interval=%d)\n",
            ckpt_config$file, ckpt_config$sim_interval))

slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")
n_cores <- if (nzchar(slurm_cpus)) {
  max(1L, as.integer(slurm_cpus))
} else {
  max(1L, parallel::detectCores(logical = FALSE) - 1L)
}
n_cores <- min(n_cores, n_sims)
cat(sprintf("[diag/c2] using %d cores\n", n_cores))

t0 <- Sys.time()
results <- run_simulation_study(
  n_sims          = n_sims,
  n_total_vec     = n_total,
  K_vec           = K,
  p_vec           = p,
  configs         = config,
  n_cores         = n_cores,
  checkpoint_config = ckpt_config,
  nlambda_init    = as.integer(Sys.getenv("NLAMBDA_INIT", "100")),
  nested_parallel = TRUE,
  estimand_type   = estimand,
  site_allocation = allocation,
  transform_type  = transform,
  outcome_type    = "binary",
  heterogeneity_type = heterogeneity,
  shift_strength  = shift,
  n_folds         = n_folds,
  use_lambda_cache = TRUE,
  verbose_every   = 25L,
  parallel_strategy = "outer_priority",
  estimate_ate    = FALSE,
  dgp_type        = "roce"
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("[diag/c2] simulation elapsed = %.2f min\n", elapsed))

summary_df <- summarize_results(results)
summary_df$tag <- tag
summary_df$allocation <- allocation
summary_df$transform <- transform
summary_df$shift <- shift
summary_df$n_folds <- n_folds
summary_df$elapsed_min <- round(elapsed, 3)

write.csv(results, file = csv_path, row.names = FALSE)
write.csv(summary_df, file = summary_path, row.names = FALSE)
save(results, summary_df, file = rdata_path)

cleanup_checkpoint(ckpt_config$file,
                   ckpt_config$preempt_signal,
                   ckpt_config$saved_signal)

cat(sprintf("[diag/c2] wrote %s\n", summary_path))
cat("\n========== C2 diagnosis summary ==========\n")
print(summary_df[, intersect(
  c("method", "bias_mean", "bias_sd", "se_mean", "coverage", "n_success",
    "estimand_type", "K", "p", "n_folds", "tag"),
  names(summary_df)
)], row.names = FALSE)
cat("[diag/c2] DONE\n")
