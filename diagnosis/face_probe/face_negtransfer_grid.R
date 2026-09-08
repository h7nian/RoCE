#!/usr/bin/env Rscript
# FACE.pdf-style negative-transfer experiment (data for the Fig-1 analog).
# Sweeps a 1-D deviation gradient (one source made increasingly non-informative)
# and records, per method, the TATE Bias / RMSE / Coverage / 95%-CI length --
# the four FACE Fig.1 panels. Methods: our cross-fit (one/two_round) + the FACE
# baselines (target_only, SS, IVW, tilted-AIPW). Uses the INSTALLED package, so
# it reflects whatever truncation radius is current; run AFTER the final M_tau is
# locked + reinstalled. Writes a tidy CSV consumed by the plotting script.
suppressPackageStartupMessages(library(RoCE))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "200"))
pp      <- as.integer(Sys.getenv("ROCE_PROBE_P", "10"))
cfg     <- Sys.getenv("ROCE_PROBE_CONFIG", "C1")
KK      <- as.integer(Sys.getenv("ROCE_PROBE_K", "2"))
# Hold per-site sample size fixed as K varies (the fair FACE-style comparison:
# more source sites, same n per site).
n_site  <- as.integer(Sys.getenv("ROCE_PROBE_NSITE", "1000"))
n_total <- n_site * (KK + 1L)
# Deviation gradient: one source deviates by rho on the binary log-odds scale
# (n_dev=1), rho increasing -> stronger negative transfer. rho=0 = all-informative.
rho_grid <- c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5)

message(sprintf("[negT] FACE Fig1 data | p=%d %s K=%d n=%d (%d/site) kf=5 binary TATE | nsims=%d cores=%d | rho in {%s}",
                pp, cfg, KK, n_total, n_site, n_sims, n_cores, paste(rho_grid, collapse=",")))

all_summ <- list()
for (rho in rho_grid) {
  t0 <- Sys.time()
  res <- run_simulation_study(
    n_sims = n_sims, n_total_vec = n_total, K_vec = KK, p_vec = pp,
    configs = cfg, n_cores = n_cores, nested_parallel = FALSE,
    outcome_type = "binary", n_folds = 5L, dgp_type = "face",
    estimate_ate = TRUE,                 # TATE = mu1 - mu0 (FACE Fig.1 estimand)
    ate_deviation = rho, n_deviated_sites = if (rho > 0) 1L else 0L
  )
  summ <- summarize_results(res)
  summ$rho <- rho
  summ$K <- KK
  summ$n_deviated <- if (rho > 0) 1L else 0L
  all_summ[[length(all_summ) + 1]] <- summ
  message(sprintf("  rho=%.1f done in %.1f min", rho,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- do.call(rbind, all_summ)
out_csv <- sprintf("diagnosis/face_probe/validation/negtransfer_%s_p%d_K%d.csv", cfg, pp, KK)
write.csv(out, out_csv, row.names = FALSE)
message(sprintf("[negT] wrote %s (%d rows)", out_csv, nrow(out)))

# Quick console view of the key panels
ord <- c("target_only","sample_size","inverse_variance","tilted_aipw",
         "one_round_crossfit","two_round_crossfit")
for (m in ord) {
  s <- out[out$method == m, ]
  if (nrow(s) == 0) next
  message(sprintf("\n  %-18s  rho:  %s", m,
                  paste(sprintf("%.1f", s$rho), collapse = "    ")))
  message(sprintf("    bias        %s", paste(sprintf("%+.3f", s$bias_mean), collapse = "  ")))
  message(sprintf("    rmse        %s", paste(sprintf("%.3f", s$rmse), collapse = "   ")))
  message(sprintf("    coverage    %s", paste(sprintf("%.2f", s$coverage), collapse = "   ")))
  message(sprintf("    ci_width    %s", paste(sprintf("%.3f", s$ci_width_mean), collapse = "  ")))
}
message("\n[negT] DONE")
