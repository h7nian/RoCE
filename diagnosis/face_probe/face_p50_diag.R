#!/usr/bin/env Rscript
# Verify the graceful-degradation instrumentation end-to-end: run a couple of
# FACE binary p=50 seeds (the saturating setting that needs the fallback),
# print the per-sim degenerate-OR-fold count, and confirm summarize_results()
# surfaces the per-setting rate column.
suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)

methods <- c("target_only", "federated_dr", "tilted_aipw")
rows <- list()
for (sid in 1:2) {
  r <- run_single_simulation(
    sim_id = sid, n_total = 3000L, K = 2L, p = 50L, config = "C1",
    methods = methods, verbose = FALSE, n_cores_internal = 1L,
    outcome_type = "binary", n_folds = 5L, dgp_type = "face"
  )
  cat(sprintf("[diag] sim %d: or_degenerate_folds = %d | methods = %s\n",
              sid, r$or_degenerate_folds[1L], paste(unique(r$method), collapse = ",")))
  rows[[sid]] <- r
}

res  <- do.call(rbind, rows)
summ <- summarize_results(res)
cat("\n[diag] summary columns: ", paste(names(summ), collapse = ", "), "\n", sep = "")
show_cols <- intersect(c("method", "config", "coverage", "bias_mean",
                         "or_degenerate_folds_mean"), names(summ))
print(summ[, show_cols])
cat("\n[diag] DONE -- graceful-degradation rate surfaced in the summary.\n")
