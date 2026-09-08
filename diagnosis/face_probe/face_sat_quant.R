#!/usr/bin/env Rscript
# Quantify the binary FACE-DGP outcome saturation across dimension. Reports the
# target-arm potential-outcome prevalences and the treated-arm minority-class
# count (what the OR model actually sees, and what drives its degeneracy). This
# isolates the *design* issue (linear-predictor scale) from estimation.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

for (p in c(10L, 50L, 100L)) {
  set.seed(1)
  d <- generate_face_data(n_total = 6000L, K = 2L, p = p, config = "C1",
                          estimand_type = "sample", outcome_type = "binary")
  t_idx <- d$R == "t"
  tr    <- which(d$A == 1L)
  n1 <- sum(d$Y[tr] == 1L); n0 <- sum(d$Y[tr] == 0L)
  message(sprintf(
    "[sat] p=%3d | target E[Y(1)]=%.3f E[Y(0)]=%.3f | treated arm: n=%d, ones=%d zeros=%d (minority=%d)",
    p, mean(d$Y_1[t_idx]), mean(d$Y_0[t_idx]), length(tr), n1, n0, min(n1, n0)))
}
message("[sat] DONE")
