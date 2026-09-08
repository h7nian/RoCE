#!/usr/bin/env Rscript
# Bracket / localize the one_round_crossfit slowdown on the FACE DGP. Runs
# one_round alone (single-core, so the run_crossfit [stage] timers stream to
# stderr) at a configurable dimension and outcome. Parameters via env:
#   ROCE_PROBE_P       (default 50)
#   ROCE_PROBE_OUTCOME (default "continuous")
#   ROCE_PROBE_NFOLDS  (default 5)
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

pp <- as.integer(Sys.getenv("ROCE_PROBE_P", "50"))
oc <- Sys.getenv("ROCE_PROBE_OUTCOME", "continuous")
nf <- as.integer(Sys.getenv("ROCE_PROBE_NFOLDS", "5"))

message(sprintf("[probe] one_round | FACE p=%d C1 %s nfolds=%d single-core | %s",
                pp, oc, nf, format(Sys.time())))
t0 <- Sys.time()
r <- run_single_simulation(sim_id = 1L, n_total = 3000L, K = 2L, p = pp,
       config = "C1", methods = c("one_round_crossfit"), verbose = FALSE,
       n_cores_internal = 1L, outcome_type = oc, n_folds = nf, dgp_type = "face")
message(sprintf("[probe] one_round DONE in %.1fs (%d rows)",
                as.numeric(difftime(Sys.time(), t0, units = "secs")), nrow(r)))
