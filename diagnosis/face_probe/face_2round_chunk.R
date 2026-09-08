#!/usr/bin/env Rscript
# 1-round vs 2-round under source-only effect modification (EM). All sources are
# transportable (ate_deviation=0); the only source-vs-target difference is the
# EM *shape* (effect_mod_strength), the regime that separates the two algorithms.
# Loads the EM-enabled RoCE from ~/Rlibs_em so the main ~/Rlibs (p=100 jobs) is
# untouched. Forks over sim_ids; records mu1 and TATE rows for both algorithms.
.libPaths(c(path.expand("~/Rlibs_em"), .libPaths()))
suppressPackageStartupMessages({library(RoCE); library(parallel)})
K     <- as.integer(Sys.getenv("CK_K", "4"))
CFG   <- Sys.getenv("CK_CFG", "C1")
P     <- as.integer(Sys.getenv("CK_P", "10"))
EM    <- as.numeric(Sys.getenv("CK_EM", "0"))
NSITE <- as.integer(Sys.getenv("CK_NSITE", "1000"))
NFOLD <- as.integer(Sys.getenv("CK_NFOLDS", "5"))
S0    <- as.integer(Sys.getenv("CK_SIM_START", "1"))
S1    <- as.integer(Sys.getenv("CK_SIM_END", "20"))
OUTCOME <- Sys.getenv("CK_OUTCOME", "binary")   # "binary" or "continuous"
RHO   <- as.numeric(Sys.getenv("CK_RHO", "0"))  # source ATE deviation (non-transportability)
NDEV  <- as.integer(Sys.getenv("CK_NDEV", "0")) # number of non-transportable sources
METHODS <- c("one_round_crossfit", "two_round_crossfit", "target_only")
n_cores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
n_total <- NSITE * (K + 1L)
message(sprintf("[2r] K=%d %s p=%d em=%g n=%d kf=%d sims=%d-%d cores=%d | %s",
                K, CFG, P, EM, n_total, NFOLD, S0, S1, n_cores, format(Sys.time())))
# CK_NCORES_INTERNAL>1 parallelizes each sim over its K source-site fits (used
# for the slow high-dim p=50 cells with chunk size 1); otherwise fork over sims.
NCI <- max(1L, as.integer(Sys.getenv("CK_NCORES_INTERNAL", "1")))
chunk_size <- S1 - S0 + 1L
one <- function(s) tryCatch(
  run_single_simulation(sim_id = s, n_total = n_total, K = K, p = P, config = CFG,
    outcome_type = OUTCOME, estimand_type = "superpopulation", n_folds = NFOLD,
    dgp_type = "face", ate_deviation = RHO, n_deviated_sites = NDEV,
    effect_mod_strength = EM, methods = METHODS,
    estimate_ate = TRUE, verbose = FALSE, n_cores_internal = NCI),
  error = function(e) { message("  sim ", s, " ERROR: ", conditionMessage(e)); NULL })
res <- do.call(rbind, mclapply(S0:S1, one, mc.cores = if (chunk_size > 1L) n_cores else 1L))
if (is.null(res) || nrow(res) == 0L) { message("[2r] no results"); quit(status = 1L) }
res$em <- EM; res$rho <- RHO; res$n_deviated <- NDEV
subdir <- if (NDEV > 0L) "chunks2r_nt" else if (OUTCOME == "continuous") "chunks2r_cont" else "chunks2r"
dir.create(file.path("diagnosis/face_probe/validation", subdir), showWarnings = FALSE)
out <- if (NDEV > 0L) {
  sprintf("diagnosis/face_probe/validation/%s/nt_K%d_%s_p%d_rho%g_nd%d_%05d-%05d.csv",
          subdir, K, CFG, P, RHO, NDEV, S0, S1)
} else {
  sprintf("diagnosis/face_probe/validation/%s/em_K%d_%s_p%d_em%g_%05d-%05d.csv",
          subdir, K, CFG, P, EM, S0, S1)
}
write.csv(res, out, row.names = FALSE)
message(sprintf("[2r] wrote %s (%d rows) | %s", out, nrow(res), format(Sys.time())))
