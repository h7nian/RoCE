suppressPackageStartupMessages(library(RoCE))
cat("CV_MAX_ITER in installed pkg:", tryCatch(RoCE:::CV_MAX_ITER, error=function(e) "(C++ side)"), "\n")
set.seed(1)
d  <- RoCE:::generate_face_data(3000L, K=2L, p=50L, config="C1", outcome_type="binary",
        estimand_type="superpopulation", ate_deviation=1.0, n_deviated_sites=1L)
ds <- RoCE:::split_data_by_site(d)
for (nf in c(2L, 3L)) {
  t0 <- Sys.time()
  r <- RoCE:::run_crossfit(ds, n_folds=nf, communication_mode="one_round", verbose=FALSE,
         family="binomial", M_tau_inference=5)
  cat(sprintf("[our method] run_crossfit one_round p=50 n_folds=%d: %.1f min (est=%.3f)\n",
      nf, as.numeric(difftime(Sys.time(),t0,units="mins")), r$estimate))
}
# full panel (all baselines) at n_folds=2
t0 <- Sys.time()
rs <- RoCE::run_single_simulation(sim_id=1, n_total=3000, K=2, p=50, config="C1",
        outcome_type="binary", estimand_type="superpopulation", n_folds=2, dgp_type="face",
        ate_deviation=1.0, n_deviated_sites=1, estimate_ate=TRUE, verbose=FALSE, n_cores_internal=1)
cat(sprintf("[full panel] run_single_simulation p=50 n_folds=2 (all baselines+ATE): %.1f min, %d rows\n",
    as.numeric(difftime(Sys.time(),t0,units="mins")), nrow(rs)))
