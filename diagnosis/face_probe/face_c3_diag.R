#!/usr/bin/env Rscript
# Reproduce the C3 failures on the CURRENT source (load_all), single-core so the
# error halts with a real traceback (no clusterApply wrapping). Full default
# method set, at the exact failing cell (face DGP, p=50, C3). Prints the first
# error location so we can fix it, then re-run.
suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)
options(error = function() { traceback(2); quit(status = 1) })

for (sid in 1:6) {
  cat(sprintf("[c3] sim %d: face C3 p=50 n=3000 kf=5 (full methods) | %s\n",
              sid, format(Sys.time())))
  r <- run_single_simulation(
    sim_id = sid, n_total = 3000L, K = 2L, p = 50L, config = "C3",
    methods = c("one_round_crossfit"),
    verbose = FALSE, n_cores_internal = 1L, outcome_type = "binary",
    n_folds = 5L, dgp_type = "face"
  )
  cat(sprintf("[c3] sim %d OK: %d rows | or_degenerate_folds=%s | methods=%s\n",
              sid, nrow(r),
              if ("or_degenerate_folds" %in% names(r)) r$or_degenerate_folds[1L] else "NA",
              paste(unique(r$method), collapse = ",")))
}
cat("[c3] ALL sims passed on current source.\n")
