#!/usr/bin/env Rscript
# One small chunk of the negative-transfer study: runs sim_ids CK_SIM_START..END
# for a single (K, config, p, rho, n_dev) cell and writes the RAW per-sim/per-method
# rows. Many such chunks are submitted as a job array across partitions and then
# combined + summarized by face_negT_aggregate.R. Uses the installed RoCE
# (M_tau=5); forks over sim_ids within the chunk.
suppressPackageStartupMessages(library(RoCE))
suppressPackageStartupMessages(library(parallel))

K     <- as.integer(Sys.getenv("CK_K", "2"))
CFG   <- Sys.getenv("CK_CFG", "C1")
P     <- as.integer(Sys.getenv("CK_P", "10"))
RHO   <- as.numeric(Sys.getenv("CK_RHO", "0"))
NDEV  <- as.integer(Sys.getenv("CK_NDEV", "0"))
NSITE <- as.integer(Sys.getenv("CK_NSITE", "1000"))
NFOLD <- as.integer(Sys.getenv("CK_NFOLDS", "5"))
S0    <- as.integer(Sys.getenv("CK_SIM_START", "1"))
S1    <- as.integer(Sys.getenv("CK_SIM_END", "20"))
EATE  <- toupper(Sys.getenv("CK_ESTIMATE_ATE", "TRUE")) %in% c("TRUE", "T", "1", "YES")
# Methods to compute (default: lean set used by the negative-transfer figure --
# only the 1-round RoCE plus the penalized baselines; 2-round and tilted_aipw
# are dropped to halve per-sim cost, and tilted_aipw is degenerate in high dim).
METHODS <- trimws(strsplit(Sys.getenv("CK_METHODS",
  "one_round_crossfit,target_only,sample_size,inverse_variance,federated_dr,pooled_dr"),
  ",")[[1]])
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(n_cores) || n_cores < 1L) n_cores <- 1L

n_total <- NSITE * (K + 1L)
message(sprintf("[chunk] K=%d %s p=%d rho=%g nd=%d n=%d kf=%d sims=%d-%d cores=%d ate=%s methods=%s | %s",
                K, CFG, P, RHO, NDEV, n_total, NFOLD, S0, S1, n_cores, EATE,
                paste(METHODS, collapse="+"), format(Sys.time())))

one <- function(s) tryCatch(
  run_single_simulation(sim_id = s, n_total = n_total, K = K, p = P, config = CFG,
    outcome_type = "binary", estimand_type = "superpopulation", n_folds = NFOLD,
    dgp_type = "face", ate_deviation = RHO, n_deviated_sites = NDEV, methods = METHODS,
    estimate_ate = EATE, verbose = FALSE, n_cores_internal = 1L),
  error = function(e) { message("  sim ", s, " ERROR: ", conditionMessage(e)); NULL })

res <- do.call(rbind, mclapply(S0:S1, one, mc.cores = n_cores))
if (is.null(res) || nrow(res) == 0L) { message("[chunk] no results"); quit(status = 1L) }
# Tag with the deviation grid coordinates for grouping during aggregation.
res$rho <- RHO
res$n_deviated <- NDEV

out <- sprintf("diagnosis/face_probe/validation/chunks/negT_K%d_%s_p%d_rho%g_nd%d_%05d-%05d.csv",
               K, CFG, P, RHO, NDEV, S0, S1)
write.csv(res, out, row.names = FALSE)
message(sprintf("[chunk] wrote %s (%d rows) | %s", out, nrow(res), format(Sys.time())))
