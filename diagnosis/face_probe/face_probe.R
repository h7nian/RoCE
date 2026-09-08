#!/usr/bin/env Rscript
# Diagnosis-only probe: FACE-paper DGP (dgp_type = "face") at calibrated, modest
# per-site sample sizes, to check whether RoCE attains valid coverage in-regime
# and whether the cross-fitting fold count (kf) still matters under the current
# fold-summed calibration.
#
# Grid (K = 2, so n_k = n_total / (K+1) = n_total / 3):
#   n_total in {3000, 6000}  -> n_k in {1000, 2000}
#   p       in {50, 100}
#   config  in {C1, C2, C3}  (FACE Settings 1/2/3)
#   outcome  = binary, ate_deviation = 0 (all sources informative: operating-range check)
# kf is the single CLI argument (run once per fold count, e.g. 5 and 10).
#
# Reuses the production pipeline (run_simulation_study + summarize_results); no
# package-source or main.sh changes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 3L) {
  stop("usage: face_probe.R <n_folds> [ate_deviation] [n_deviated_sites]", call. = FALSE)
}
kf            <- as.integer(args[[1L]])
ate_deviation <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 0.0
n_deviated    <- if (length(args) >= 3L) as.integer(args[[3L]]) else 0L
stopifnot(is.finite(kf), kf >= 3L, is.finite(ate_deviation),
          ate_deviation >= 0, n_deviated >= 0L)

suppressPackageStartupMessages(library(RoCE))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

out_dir <- file.path("diagnosis", "face_probe", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_sims <- 200L
cat(sprintf("[face_probe] kf=%d ate_dev=%g n_dev=%d n_sims=%d n_cores=%d | %s\n",
            kf, ate_deviation, n_deviated, n_sims, n_cores, format(Sys.time())))

results <- run_simulation_study(
  n_sims           = n_sims,
  n_total_vec      = c(3000L, 6000L),   # n_k = 1000, 2000 at K=2
  K_vec            = 2L,
  p_vec            = c(50L, 100L),
  configs          = c("C1", "C2", "C3"),
  n_cores          = n_cores,
  nested_parallel  = TRUE,
  outcome_type     = "binary",
  n_folds          = kf,
  dgp_type         = "face",
  ate_deviation    = ate_deviation,
  n_deviated_sites = n_deviated
)

summary_df <- summarize_results(results)
summary_df$kf            <- kf
summary_df$ate_deviation <- ate_deviation
summary_df$n_deviated    <- n_deviated
summary_df$n_k           <- summary_df$n_total %/% 3L   # K=2 -> 3 sites

out_file <- file.path(out_dir, sprintf("face_probe_kf%d_dev%g_nd%d.csv",
                                       kf, ate_deviation, n_deviated))
write.csv(summary_df, out_file, row.names = FALSE)
cat(sprintf("[face_probe] wrote %s (%d rows) | %s\n",
            out_file, nrow(summary_df), format(Sys.time())))
