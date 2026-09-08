#!/usr/bin/env Rscript
# Measure the population density-ratio tilt bound M_0 = ess-sup |phi^T gamma*| on
# the binary FACE DGP. As n grows, the fitted gamma-hat -> gamma*, so the largest
# |phi^T gamma-hat| over a large source sample approximates M_0 (the constant in
# the bounded-tilt / positivity assumption that justifies SMMAL-style truncation).
# Reporting the WHOLE distribution (not just the max) separates the population
# bound M_0 (the bulk) from finite-sample excursions (the tail). Uses run_crossfit
# with M_tau_inference=Inf so nothing is clipped; reads max_abs_logit from the
# clip diagnostics. Diagnostic only.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)
suppressPackageStartupMessages(library(parallel))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L
pp  <- as.integer(Sys.getenv("ROCE_PROBE_P", "10"))
cfg <- Sys.getenv("ROCE_PROBE_CONFIG", "C1")
n_grid  <- c(3000L, 9000L, 27000L)
n_reps  <- as.integer(Sys.getenv("ROCE_NSIMS", "24"))

message(sprintf("[tilt] p=%d %s binary | n in {%s} | reps=%d cores=%d | M_tau_inf=Inf (no clip)",
                pp, cfg, paste(n_grid, collapse=","), n_reps, n_cores))

probe <- function(args) {
  n <- args$n; rep <- args$rep
  set.seed(1000L * rep + n)
  d  <- generate_face_data(n_total = n, K = 2L, p = pp, config = cfg,
                           outcome_type = "binary", estimand_type = "superpopulation")
  ds <- split_data_by_site(d)
  r  <- tryCatch(run_crossfit(ds, n_folds = 5L, communication_mode = "one_round",
                              lambda_selection = "cv", verbose = FALSE,
                              family = "binomial", M_tau_inference = Inf),
                 error = function(e) NULL)
  if (is.null(r) || is.null(r$correction_clip_diagnostics)) return(NULL)
  data.frame(n = n, rep = rep,
             max_abs_logit  = r$correction_clip_diagnostics$max_abs_logit,
             max_raw_weight = r$correction_clip_diagnostics$max_raw_weight)
}
grid <- do.call(rbind, lapply(n_grid, function(n) data.frame(n = n, rep = seq_len(n_reps))))
tasks <- split(grid, seq_len(nrow(grid)))
res <- do.call(rbind, mclapply(tasks, probe, mc.cores = n_cores))

message("\n      n |  median max|logit|   p90    max   | median max_weight   p90")
for (n in n_grid) {
  s <- res[res$n == n & is.finite(res$max_abs_logit), ]
  if (nrow(s) < 3) { message(sprintf("  %6d  (too few: %d)", n, nrow(s))); next }
  ql <- quantile(s$max_abs_logit, c(0.5, 0.9), na.rm = TRUE)
  qw <- quantile(s$max_raw_weight, c(0.5, 0.9), na.rm = TRUE)
  message(sprintf("  %6d |    %6.2f          %6.2f  %6.2f |   %8.1f       %8.1f",
                  n, ql[1], ql[2], max(s$max_abs_logit), qw[1], qw[2]))
}
message("\n[tilt] As n grows, max|logit| should stabilize near M_0 (population bound).")
message("[tilt] DONE")
