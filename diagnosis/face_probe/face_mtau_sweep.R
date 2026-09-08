#!/usr/bin/env Rscript
# Decisive M_tau_inference sweep on the failing cell (p=10 C1 binary), to (i)
# MEASURE the real density-ratio tilt scale (max|phi^T gamma|) -- which tells us
# where SMMAL-style truncation should sit -- and (ii) find the M_tau that
# collapses empSE onto estSE WITHOUT introducing bias (bias ~ 0 => truncation
# inactive at the truth, i.e. the bounded-tilt assumption holds for that radius).
# Runs run_crossfit directly via load_all + fork (mclapply), passing M_tau_inference
# explicitly, so it is independent of the installed constant. Each sim's data is
# reused across all M values (paired comparison).
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(parallel))

n_sims  <- as.integer(Sys.getenv("ROCE_NSIMS", "150"))
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
cfg     <- Sys.getenv("ROCE_PROBE_CONFIG", "C1")
pp      <- as.integer(Sys.getenv("ROCE_PROBE_P", "10"))
mode    <- Sys.getenv("ROCE_PROBE_MODE", "one_round")
M_grid  <- c(Inf, 9, 7, 5, 4, 3, 2)

message(sprintf("[mtau] %s | p=%d %s binary n=3000 kf=5 | nsims=%d cores=%d | M_tau_inf in {%s}",
                mode, pp, cfg, n_sims, n_cores, paste(M_grid, collapse=",")))

one_sim <- function(sim) {
  set.seed(sim)
  d  <- generate_face_data(n_total = 3000L, K = 2L, p = pp, config = cfg,
                           outcome_type = "binary", estimand_type = "superpopulation")
  ds <- split_data_by_site(d)
  truth <- d$mu1_true
  rows <- lapply(M_grid, function(M) {
    r <- tryCatch(run_crossfit(ds, n_folds = 5L, communication_mode = mode,
                               lambda_selection = "cv", verbose = FALSE,
                               family = "binomial", M_tau_inference = M),
                  error = function(e) NULL)
    if (is.null(r)) return(NULL)
    cd <- r$correction_clip_diagnostics
    data.frame(sim = sim, M = M, est = r$estimate, se = r$se, truth = truth,
               max_abs_logit = if (!is.null(cd)) cd$max_abs_logit else NA_real_,
               max_raw_weight = if (!is.null(cd)) cd$max_raw_weight else NA_real_)
  })
  do.call(rbind, rows)
}

res <- do.call(rbind, mclapply(seq_len(n_sims), one_sim, mc.cores = n_cores))

message("\n M_tau_inf |   bias    empSE    estSE  emp/est | cover |  median max|logit|  p95 max|logit|")
for (M in M_grid) {
  s <- res[res$M == M & is.finite(res$est) & is.finite(res$se), ]
  if (nrow(s) < 10) { message(sprintf("  %6.0f   (too few: %d)", M, nrow(s))); next }
  bias  <- mean(s$est - s$truth)
  empSE <- sd(s$est)
  estSE <- mean(s$se)
  cover <- mean(s$truth >= s$est - 1.96*s$se & s$truth <= s$est + 1.96*s$se)
  mlog  <- quantile(s$max_abs_logit, c(0.5, 0.95), na.rm = TRUE)
  message(sprintf("  %6s | %+0.4f  %0.4f  %0.4f   %0.2f  | %0.3f |    %6.2f          %6.2f",
                  ifelse(is.finite(M), sprintf("%.0f", M), "Inf"),
                  bias, empSE, estSE, empSE/estSE, cover, mlog[1], mlog[2]))
}
message("\n[mtau] DONE")
