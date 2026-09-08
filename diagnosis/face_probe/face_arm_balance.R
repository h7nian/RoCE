#!/usr/bin/env Rscript
# Why does fit_initial_outcome still hit a degenerate binomial class on the FIXED
# DGP? The data fix standardizes the logit signal using the TARGET law, but source
# sites have nu=0.2 skewed covariates (a covariate shift), so a source arm's
# prevalence can still be extreme. Measure, per site x arm, the prevalence and
# minority-class count -- at full-site level AND at a cross-fit fold subset
# (n_arm~250, which is what the OR model actually sees and what cv.glmnet's
# internal 10-fold CV needs >=~10 minority for).
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

report <- function(cfg, seed) {
  set.seed(seed)
  d <- generate_face_data(n_total = 3000L, K = 2L, p = 10L, config = cfg,
                          outcome_type = "binary", estimand_type = "superpopulation")
  for (site in c("t", "s1", "s2")) {
    for (a in c(0L, 1L)) {
      idx <- which(d$R == site & d$A == a)
      if (length(idx) == 0) next
      y <- d$Y[idx]; n1 <- sum(y == 1L); n0 <- sum(y == 0L); mn <- min(n1, n0)
      # what a single 5-fold cross-fit fold subset of this arm looks like
      fold_n   <- round(length(idx) * 0.8)          # OR fit uses ~4/5 folds
      fold_min <- round(mn * 0.8)
      flag <- if (fold_min < 10L) "  <-- cv.glmnet 10-fold UNSAFE" else ""
      message(sprintf("  %s | site=%-2s A=%d | n=%4d prev=%.3f minority=%4d | fold(n~%d): minority~%d%s",
                      cfg, site, a, length(idx), mean(y), mn, fold_n, fold_min, flag))
    }
  }
}

for (cfg in c("C1", "C2", "C3")) {
  message(sprintf("=== config %s (3 seeds) ===", cfg))
  for (s in 1:3) report(cfg, s)
}
message("[arm-balance] DONE")
