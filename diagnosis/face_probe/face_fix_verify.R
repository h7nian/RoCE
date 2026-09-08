#!/usr/bin/env Rscript
# Verify the saturation fixes unblock every config x dimension on the CURRENT
# source (load_all). One full-method simulation per (config, p) cell; tryCatch
# per cell so a failure is reported without stopping the rest. Exits non-zero if
# any cell crashes.
suppressPackageStartupMessages(library(devtools))
load_all(".", quiet = TRUE)

grid <- expand.grid(config = c("C1", "C2", "C3"), p = c(50L, 100L),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
n_fail <- 0L
for (i in seq_len(nrow(grid))) {
  cfg <- grid$config[i]; pp <- grid$p[i]
  cat(sprintf("[verify] %s p=%d (full methods, 4 inner cores) | %s\n",
              cfg, pp, format(Sys.time())))
  r <- tryCatch(
    run_single_simulation(
      sim_id = 1L, n_total = 3000L, K = 2L, p = pp, config = cfg,
      methods = c("one_round_crossfit", "two_round_crossfit",
                  "target_only", "tilted_aipw", "federated_dr"),
      verbose = FALSE, n_cores_internal = 4L, outcome_type = "binary",
      n_folds = 5L, dgp_type = "face"
    ),
    error = function(e) e
  )
  if (inherits(r, "error")) {
    cat(sprintf("[verify] FAIL %s p=%d: %s\n", cfg, pp, conditionMessage(r)))
    n_fail <- n_fail + 1L
    next
  }
  cat(sprintf("[verify] OK   %s p=%d: %d rows | methods=%s\n",
              cfg, pp, nrow(r), paste(unique(r$method), collapse = ",")))
}
cat(sprintf("[verify] DONE: %d/%d cells failed.\n", n_fail, nrow(grid)))
if (n_fail > 0L) quit(status = 1)
