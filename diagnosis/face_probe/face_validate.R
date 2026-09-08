#!/usr/bin/env Rscript
# One-cell FACE validation against the installed RoCE package.
# CLI args: kf  n_total  p  config  ate_deviation  n_deviated
# Runs ROCE_VAL_NSIMS replications of the FACE (negative-transfer) DGP for a
# single grid cell and writes a per-cell summary CSV (incl. coverage, bias, and
# the or_degenerate_folds_mean saturation diagnostic). Arrayed over all cells by
# face_validate.cmd so the full grid runs in parallel.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 6L)
kf      <- as.integer(args[[1L]])
n_total <- as.integer(args[[2L]])
p       <- as.integer(args[[3L]])
config  <- args[[4L]]
ate_dev <- as.numeric(args[[5L]])
n_dev   <- as.integer(args[[6L]])

suppressPackageStartupMessages(library(RoCE))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_VAL_NSIMS", "50"))

out_dir <- file.path("diagnosis", "face_probe", "validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("[validate] kf=%d n=%d p=%d %s dev=%g nd=%d | nsims=%d cores=%d | %s\n",
            kf, n_total, p, config, ate_dev, n_dev, n_sims, n_cores, format(Sys.time())))

res <- run_simulation_study(
  n_sims           = n_sims,
  n_total_vec      = n_total,
  K_vec            = 2L,
  p_vec            = p,
  configs          = config,
  n_cores          = n_cores,
  nested_parallel  = TRUE,
  outcome_type     = "binary",
  n_folds          = kf,
  dgp_type         = "face",
  ate_deviation    = ate_dev,
  n_deviated_sites = n_dev
)

summ <- summarize_results(res)
summ$kf            <- kf
summ$p_cell        <- p
summ$ate_deviation <- ate_dev
summ$n_deviated    <- n_dev

out_file <- file.path(out_dir, sprintf("cell_kf%d_n%d_p%d_%s_dev%g_nd%d.csv",
                                       kf, n_total, p, config, ate_dev, n_dev))
write.csv(summ, out_file, row.names = FALSE)
cat(sprintf("[validate] wrote %s (%d rows) | %s\n", out_file, nrow(summ), format(Sys.time())))
