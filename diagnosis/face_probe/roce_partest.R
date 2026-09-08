suppressPackageStartupMessages(library(RoCE))
nc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","8"))
cat("cores available:", nc, "\n")
t0 <- Sys.time()
rs <- RoCE::run_single_simulation(sim_id=1, n_total=3000, K=2, p=50, config="C1",
        outcome_type="binary", estimand_type="superpopulation", n_folds=5, dgp_type="face",
        ate_deviation=1.0, n_deviated_sites=1, estimate_ate=TRUE, verbose=FALSE,
        n_cores_internal=nc)
cat(sprintf("p=50 n_folds=5 ALL baselines+ATE, n_cores_internal=%d: %.1f min, %d rows\n",
    nc, as.numeric(difftime(Sys.time(),t0,units="mins")), nrow(rs)))
