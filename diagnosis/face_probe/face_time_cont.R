#!/usr/bin/env Rscript
# Per-method timing on C1 p=50 (single-core) to isolate which estimator is slow
# or hangs after the saturation fixes. message() -> stderr so progress is
# visible even if the job is killed (cat to a redirected file is buffered).
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)
methods <- c("target_only", "tilted_aipw", "federated_dr",
             "one_round_crossfit", "two_round_crossfit")
for (m in methods) {
  message(sprintf("[time] START %s | %s", m, format(Sys.time())))
  t0 <- Sys.time()
  r <- tryCatch(
    run_single_simulation(sim_id = 1L, n_total = 3000L, K = 2L, p = 50L,
      config = "C1", methods = m, verbose = FALSE, n_cores_internal = 1L,
      outcome_type = "continuous", n_folds = 5L, dgp_type = "face"),
    error = function(e) e)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(r, "error"))
    message(sprintf("[time] %s ERROR after %.1fs: %s", m, dt, conditionMessage(r)))
  else
    message(sprintf("[time] %s OK in %.1fs | or_degen=%s", m, dt,
                    if ("or_degenerate_folds" %in% names(r)) r$or_degenerate_folds[1] else "NA"))
}
message("[time] DONE")
