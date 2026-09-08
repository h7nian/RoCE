suppressPackageStartupMessages(library(RoCE))
Sys.setenv(ROCE_AGG_DEBUG="1")
set.seed(7)
d  <- RoCE:::generate_face_data(3000L, K=2L, p=10L, config="C1", outcome_type="binary",
        estimand_type="superpopulation", ate_deviation=2.5, n_deviated_sites=1L)
ds <- RoCE:::split_data_by_site(d)
r  <- RoCE:::run_crossfit(ds, n_folds=3L, communication_mode="one_round",
        verbose=FALSE, family="binomial", M_tau_inference=5)
cat(sprintf("FINAL: eta=[%s] lambda=%.1f est=%.4f truth=%.4f\n",
    paste(round(as.numeric(r$weights),3),collapse=","), mean(as.numeric(r$fold_lambdas)),
    r$estimate, d$mu1_true))
