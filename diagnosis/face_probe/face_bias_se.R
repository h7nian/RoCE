#!/usr/bin/env Rscript
# Bias / SE / coverage of every estimator on the FIXED binary FACE DGP.
# Runs against the SOURCE tree (load_all) so it reflects the saturation fix.
# p=10 keeps it fast (the p=50 DR-CV speed fix is not applied yet); p=10 is a
# valid FACE setting (the paper's primary dimension). Reports per method/config:
#   bias  = mean(estimate - truth)
#   empSE = sd(estimate)                (empirical / Monte-Carlo SE)
#   estSE = mean(model-reported SE)     (should match empSE if calibrated)
#   cover = CI coverage (target 0.95)
#   RMSE, CI width, and the OR-degeneracy diagnostic (should be ~0 now).
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "100"))

message(sprintf("[bias/se] binary FACE | p=10 K=2 n=3000 kf=5 | C1/C2/C3 | nsims=%d cores=%d | %s",
                n_sims, n_cores, format(Sys.time())))

res <- run_simulation_study(
  n_sims          = n_sims,
  n_total_vec     = 3000L,
  K_vec           = 2L,
  p_vec           = 10L,
  configs         = c("C1", "C2", "C3"),
  n_cores         = n_cores,
  nested_parallel = FALSE,
  outcome_type    = "binary",
  n_folds         = 5L,
  dgp_type        = "face"
)
summ <- summarize_results(res)

# Order methods consistently and print a compact, readable table.
ord <- c("target_only","sample_size","inverse_variance","federated_dr",
         "pooled_dr","tilted_aipw","oracle_dr","one_round_crossfit","two_round_crossfit")
summ$method <- factor(summ$method, levels = union(ord, unique(summ$method)))
summ <- summ[order(summ$config, summ$method), ]

message("\n==================== BIAS / SE / COVERAGE (binary FACE, p=10) ====================")
for (cfg in c("C1","C2","C3")) {
  s <- summ[summ$config == cfg, ]
  if (nrow(s) == 0) next
  message(sprintf("\n--- %s  (n_success up to %d) ---", cfg, n_sims))
  message(sprintf("%-20s %8s %8s %8s %7s %8s %8s %6s",
                  "method","bias","empSE","estSE","cover","rmse","ciW","ordeg"))
  for (i in seq_len(nrow(s))) {
    message(sprintf("%-20s %8.4f %8.4f %8.4f %7.3f %8.4f %8.4f %6.2f",
                    as.character(s$method[i]), s$bias_mean[i], s$bias_sd[i],
                    s$se_mean[i], s$coverage[i], s$rmse[i], s$ci_width_mean[i],
                    if (!is.null(s$or_degenerate_folds_mean)) s$or_degenerate_folds_mean[i] else NA))
  }
}

out_csv <- "diagnosis/face_probe/validation/bias_se_binary_p10.csv"
write.csv(summ, out_csv, row.names = FALSE)
message(sprintf("\n[bias/se] wrote %s | %s", out_csv, format(Sys.time())))
