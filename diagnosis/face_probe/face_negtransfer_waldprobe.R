#!/usr/bin/env Rscript
# Wald-truncation ON/OFF probe: identical to face_negtransfer_grid.R, but
#   (1) tags the output csv with ROCE_WALD_TAG so it never clobbers the
#       production negtransfer_*.csv, and
#   (2) reads the rho grid from ROCE_RHO_GRID (comma-separated) so a smoke run
#       can use a 2-point grid.
# The Wald truncation itself is toggled NOT here but by the installed library:
# RoCE built with ROCE_AGG_WALD_LAMBDA=0.5 -> Wald-on; =1e-6 -> Wald-off
# (penalty factor (lambda*t-1)_+ == 0, i.e. pure variance-minimizing aggregation).
suppressPackageStartupMessages(library(RoCE))
cat("RoCE from:", dirname(dirname(getNamespaceInfo("RoCE", "path"))), "\n")
cat("AGG_WALD_LAMBDA =", RoCE:::AGG_WALD_LAMBDA, "\n")

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "200"))
pp      <- as.integer(Sys.getenv("ROCE_PROBE_P", "50"))
cfg     <- Sys.getenv("ROCE_PROBE_CONFIG", "C1")
KK      <- as.integer(Sys.getenv("ROCE_PROBE_K", "4"))
n_site  <- as.integer(Sys.getenv("ROCE_PROBE_NSITE", "1000"))
n_total <- n_site * (KK + 1L)
tag     <- Sys.getenv("ROCE_WALD_TAG", "on")

rho_env <- Sys.getenv("ROCE_RHO_GRID", "")
rho_grid <- if (nzchar(rho_env)) as.numeric(strsplit(rho_env, ",")[[1]]) else
  c(0.0, 0.5, 1.0, 1.5, 2.0, 2.5)

message(sprintf("[waldprobe tag=%s] p=%d %s K=%d n=%d (%d/site) kf=5 binary TATE | nsims=%d cores=%d | rho in {%s}",
                tag, pp, cfg, KK, n_total, n_site, n_sims, n_cores, paste(rho_grid, collapse=",")))

all_summ <- list()
for (rho in rho_grid) {
  t0 <- Sys.time()
  res <- run_simulation_study(
    n_sims = n_sims, n_total_vec = n_total, K_vec = KK, p_vec = pp,
    configs = cfg, n_cores = n_cores, nested_parallel = FALSE,
    outcome_type = "binary", n_folds = 5L, dgp_type = "face",
    estimate_ate = TRUE,
    ate_deviation = rho, n_deviated_sites = if (rho > 0) 1L else 0L
  )
  summ <- summarize_results(res)
  summ$rho <- rho
  summ$K <- KK
  summ$n_deviated <- if (rho > 0) 1L else 0L
  summ$wald_tag <- tag
  all_summ[[length(all_summ) + 1]] <- summ
  message(sprintf("  rho=%.1f done in %.1f min", rho,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- do.call(rbind, all_summ)
out_csv <- sprintf("diagnosis/face_probe/validation/negtransfer_%s_p%d_K%d_wald%s.csv",
                   cfg, pp, KK, tag)
write.csv(out, out_csv, row.names = FALSE)
message(sprintf("[waldprobe] wrote %s (%d rows)", out_csv, nrow(out)))

# Console view of RoCE (one_round_crossfit) + target_only across rho
for (m in c("target_only", "one_round_crossfit")) {
  s <- out[out$method == m, ]
  if (nrow(s) == 0) next
  message(sprintf("\n  %-18s [wald=%s]  rho:  %s", m, tag,
                  paste(sprintf("%.1f", s$rho), collapse = "    ")))
  message(sprintf("    bias        %s", paste(sprintf("%+.3f", s$bias_mean), collapse = "  ")))
  message(sprintf("    coverage    %s", paste(sprintf("%.2f", s$coverage), collapse = "   ")))
  message(sprintf("    ci_width    %s", paste(sprintf("%.3f", s$ci_width_mean), collapse = "  ")))
}
message("\n[waldprobe] DONE")
