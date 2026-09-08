#!/usr/bin/env Rscript
# Settle the anomaly: does config actually flow into Z_site/W_outcome and the
# estimate, in the generate_face_data -> split_data_by_site -> run_crossfit path?
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)
for (cfg in c("C1","C2","C3")) {
  set.seed(1)
  d <- generate_face_data(3000L, K=2L, p=10L, config=cfg,
                          outcome_type="binary", estimand_type="superpopulation")
  ds <- split_data_by_site(d)
  cat(sprintf("config=%s | d$Z_site=%dx%d d$W_outcome=%dx%d | split s1 Z=%d W=%d\n",
              cfg, nrow(d$Z_site), ncol(d$Z_site), nrow(d$W_outcome), ncol(d$W_outcome),
              ncol(ds$s1$Z_site), ncol(ds$s1$W_outcome)))
  r <- run_crossfit(ds, n_folds=5L, communication_mode="one_round",
                    verbose=FALSE, family="binomial", M_tau_inference=Inf)
  cat(sprintf("   one_round: estimate=%.5f se=%.5f\n", r$estimate, r$se))
}
cat("DONE\n")
