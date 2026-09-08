#!/usr/bin/env Rscript
# Why are OUR method's (one/two_round_crossfit) metrics not yet "reasonable"?
# In C1/C2 the empirical SE (~0.040) is ~2x the model SE (~0.022) even though
# coverage holds (~0.95). Diagnose the cause: is the point-estimate distribution
# HEAVY-TAILED (a few unstable cross-fit/density-ratio realizations inflate the
# marginal sd) or is the variance estimator systematically off? Compare:
#   empSE  = sd(estimate)                       (outlier-sensitive)
#   trimSE = sd after dropping 5% each tail      (robust scale)
#   madSE  = 1.4826 * median|est - median|       (robust scale)
#   estSE  = mean(model SE);  z = (est-truth)/se (per-sim standardized error)
# If trimSE ~ madSE ~ estSE << empSE, the gap is heavy tails -> stabilize the
# density-ratio weights. If the whole z-distribution is wide, it's calibration.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "200"))

message(sprintf("[metrics] our method | binary FACE p=10 K=2 n=3000 kf=5 | C1/C3 | nsims=%d cores=%d",
                n_sims, n_cores))

res <- run_simulation_study(
  n_sims = n_sims, n_total_vec = 3000L, K_vec = 2L, p_vec = 10L,
  configs = c("C1", "C3"), n_cores = n_cores, nested_parallel = FALSE,
  outcome_type = "binary", n_folds = 5L, dgp_type = "face"
)

analyze <- function(sub, cfg, m) {
  sub <- sub[is.finite(sub$estimate) & is.finite(sub$se), ]
  est <- sub$estimate; se <- sub$se; bias <- sub$bias
  if (length(est) < 10) { message(sprintf("  %s/%s: too few (%d)", cfg, m, length(est))); return(invisible()) }
  trim <- function(x, p = 0.05) { q <- quantile(x, c(p, 1 - p)); sd(x[x >= q[1] & x <= q[2]]) }
  z <- bias / se
  ord <- order(abs(bias), decreasing = TRUE)
  message(sprintf("\n  [%s / %-18s] n=%d", cfg, m, length(est)))
  message(sprintf("    empSE=%.4f  trimSE(5%%)=%.4f  madSE=%.4f  | estSE(mean)=%.4f median=%.4f",
                  sd(est), trim(est), mad(est), mean(se), median(se)))
  message(sprintf("    z=(bias/se): sd=%.2f  P(|z|>1.96)=%.3f  P(|z|>3)=%.3f  max|z|=%.2f",
                  sd(z), mean(abs(z) > 1.96), mean(abs(z) > 3), max(abs(z))))
  message(sprintf("    worst 3 |bias|: %s",
                  paste(sprintf("(bias=%.3f se=%.3f z=%.1f)", bias[ord[1:3]], se[ord[1:3]], z[ord[1:3]]),
                        collapse = " ")))
}

for (cfg in c("C1", "C3")) {
  for (m in c("target_only", "one_round_crossfit", "two_round_crossfit")) {
    analyze(res[res$config == cfg & res$method == m, ], cfg, m)
  }
}
message("\n[metrics] DONE")
