#!/usr/bin/env Rscript
# Minimal single-process reproduction of the face-DGP OR-model failure
# ("estimate_complement_fold_aipw: OR model failed: ... 'drop': non-conformable
# arguments") seen in face_probe at K=2, p=50, binary, C1. Runs n_cores=1 with a
# traceback handler so the full call stack surfaces (the parallel workers hid it).

options(error = function() { traceback(2L); quit(save = "no", status = 1L) })
suppressPackageStartupMessages(library(RoCE))

cat("[repro] face DGP | n_total=3000 K=2 p=50 C1 binary kf=5 | reproducing OR failure\n")
res <- run_simulation_study(
  n_sims           = 2L,
  n_total_vec      = 3000L,   # n_k = 1000 at K=2
  K_vec            = 2L,
  p_vec            = 50L,
  configs          = "C1",
  n_cores          = 1L,
  nested_parallel  = FALSE,
  outcome_type     = "binary",
  n_folds          = 5L,
  dgp_type         = "face",
  ate_deviation    = 0.0,
  n_deviated_sites = 0L
)
cat("[repro] NO ERROR — completed; rows:", nrow(res), "\n")
