#!/usr/bin/env Rscript
.libPaths(c(path.expand("~/Rlibs_em"), .libPaths()))
suppressPackageStartupMessages({library(RoCE); library(parallel)})
runs <- function(s, em, p) tryCatch(run_single_simulation(sim_id=s, n_total=1000*5,
   K=4, p=p, config="C1", outcome_type="binary", estimand_type="superpopulation",
   n_folds=5, dgp_type="face", ate_deviation=0, n_deviated_sites=0L,
   effect_mod_strength=em, methods=c("one_round_crossfit","two_round_crossfit","target_only"),
   estimate_ate=TRUE, verbose=FALSE, n_cores_internal=1L),
   error=function(e){message("err ",s," em=",em,": ",conditionMessage(e));NULL})
nc <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","4")))
cat("\n==== EM sanity: mu1 + TATE, p=10, K=4, 16 reps/cell ====\n")
for (em in c(0,1,2)) {
  res <- do.call(rbind, mclapply(1:16, runs, em=em, p=10, mc.cores=nc))
  if (is.null(res)) { cat("em=",em," ALL FAILED\n"); next }
  for (m in c("one_round_crossfit","two_round_crossfit","target_only",
              "one_round_crossfit_ate","two_round_crossfit_ate")) {
    r <- res[res$method==m,]; if (nrow(r)==0) next
    cat(sprintf("em=%g  %-24s n=%2d  bias=%+.4f  rmse=%.4f  cov=%.3f  ciw=%.4f\n",
      em, m, nrow(r), mean(r$bias,na.rm=T), sqrt(mean(r$bias^2,na.rm=T)),
      mean(r$coverage,na.rm=T), mean(r$ci_width,na.rm=T)))
  }
  cat("----\n")
}
