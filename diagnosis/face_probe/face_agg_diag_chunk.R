#!/usr/bin/env Rscript
# Diagnose the negative-transfer-control failure: per sim, extract the adaptive
# weights eta-hat, the source-assisted estimates mu_ts, the target-only mu_ot,
# the selected aggregation lambda, and the final estimate/truth. Reveals whether
# the aggregation shrinks the deviated source (s1) and why/why not.
suppressPackageStartupMessages(library(RoCE)); suppressPackageStartupMessages(library(parallel))
RHO  <- as.numeric(Sys.getenv("CK_RHO", "1"));  NDEV <- as.integer(Sys.getenv("CK_NDEV", "1"))
K <- as.integer(Sys.getenv("CK_K","2")); P <- as.integer(Sys.getenv("CK_P","10"))
NSITE <- as.integer(Sys.getenv("CK_NSITE","1000")); NFOLD <- as.integer(Sys.getenv("CK_NFOLDS","5"))
S0 <- as.integer(Sys.getenv("CK_SIM_START","1")); S1 <- as.integer(Sys.getenv("CK_SIM_END","15"))
nc <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1")); if(is.na(nc)||nc<1) nc<-1
n_total <- NSITE*(K+1L)
one <- function(s) tryCatch({
  set.seed(s)
  d  <- RoCE:::generate_face_data(n_total, K=K, p=P, config="C1", outcome_type="binary",
          estimand_type="superpopulation", ate_deviation=RHO, n_deviated_sites=NDEV)
  ds <- RoCE:::split_data_by_site(d)
  r  <- RoCE:::run_crossfit(ds, n_folds=NFOLD, communication_mode="one_round", verbose=FALSE,
          family="binomial", M_tau_inference=5)
  w  <- as.numeric(r$weights); se_src <- as.numeric(r$source_estimates)
  data.frame(rho=RHO, sim=s, estimate=r$estimate, se=r$se, truth=d$mu1_true,
             mu_ot=as.numeric(r$target_only$estimate), lambda=mean(as.numeric(r$fold_lambdas), na.rm=TRUE),
             eta_s1=w[1], eta_s2=if(length(w)>1) w[2] else NA_real_,
             mu_ts1=se_src[1], mu_ts2=if(length(se_src)>1) se_src[2] else NA_real_)
}, error=function(e){ message("sim ",s," ERR: ",conditionMessage(e)); NULL })
res <- do.call(rbind, mclapply(S0:S1, one, mc.cores=nc))
if (is.null(res)||nrow(res)==0) { message("no results"); quit(status=1L) }
out <- sprintf("diagnosis/face_probe/validation/chunks/aggdiag_rho%g_%05d-%05d.csv", RHO, S0, S1)
write.csv(res, out, row.names=FALSE); message(sprintf("wrote %s (%d rows)", out, nrow(res)))
