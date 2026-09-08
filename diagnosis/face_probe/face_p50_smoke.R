#!/usr/bin/env Rscript
# Fast smoke: does the FACE DGP at p=50, config C1 still trigger the
# estimate_complement_fold_aipw "OR model failed: non-conformable arguments"
# bug on the CURRENT source? (The failed validation jobs ran a stale installed
# package.) Runs a few seeds with the baseline + target-only methods that
# exercise the complement-fold AIPW path, and reports the first error (with a
# traceback) or success. Uses devtools::load_all() so it tests the working tree.

suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)

# Baselines that go through estimate_target_only_from_complement /
# estimate_complement_fold_aipw (the failing path), at the exact failing setting.
methods <- c("target_only", "federated_dr", "tilted_aipw")

n_ok <- 0L
for (sid in 1:4) {
  cat(sprintf("[p50-smoke] sim %d: face DGP, n=3000, K=2, p=50, C1, binary ...\n", sid))
  res <- tryCatch(
    run_single_simulation(
      sim_id = sid, n_total = 3000L, K = 2L, p = 50L, config = "C1",
      methods = methods, verbose = FALSE, n_cores_internal = 1L,
      outcome_type = "binary", n_folds = 5L, dgp_type = "face"
    ),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    cat(sprintf("[p50-smoke] FAILED on sim %d: %s\n", sid, conditionMessage(res)))
    cat("---- traceback ----\n"); traceback()
    quit(status = 1)
  }
  cat(sprintf("[p50-smoke] sim %d OK (%d method-rows): %s\n",
              sid, nrow(res), paste(unique(res$method), collapse = ",")))
  n_ok <- n_ok + 1L
}
cat(sprintf("[p50-smoke] ALL %d sims passed at p=50 face DGP -- non-conformable bug NOT reproduced on current source.\n", n_ok))
