#!/usr/bin/env Rscript
# Verify the binary FACE-DGP saturation fix (standardized logit signal + moderate
# log-odds ATE) AND that the continuous path is untouched. Pure data-generation
# checks (no estimation), so this is fast and single-core-safe.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

ok <- TRUE
chk <- function(cond, msg) {
  ok <<- ok && isTRUE(cond)
  message(sprintf("[%s] %s", if (isTRUE(cond)) "PASS" else "FAIL", msg))
}

message("=== (1) Binary saturation removed across dimension (sample estimand) ===")
for (p in c(10L, 50L, 100L)) {
  set.seed(1)
  d <- generate_face_data(n_total = 8000L, K = 2L, p = p, config = "C1",
                          estimand_type = "sample", outcome_type = "binary")
  t_idx <- d$R == "t"; tr <- which(d$A == 1L)
  n1 <- sum(d$Y[tr] == 1L); n0 <- sum(d$Y[tr] == 0L)
  message(sprintf(
    "  p=%3d | target E[Y(1)]=%.3f E[Y(0)]=%.3f RD=%.3f | treated minority=%d/%d (%.1f%%)",
    p, mean(d$Y_1[t_idx]), mean(d$Y_0[t_idx]),
    mean(d$Y_1[t_idx]) - mean(d$Y_0[t_idx]),
    min(n1, n0), length(tr), 100 * min(n1, n0) / length(tr)))
  chk(mean(d$Y_1[t_idx]) < 0.9 && mean(d$Y_0[t_idx]) > 0.1,
      sprintf("p=%d both arms bounded away from 0/1", p))
  chk(min(n1, n0) / length(tr) > 0.15,
      sprintf("p=%d treated minority class healthy (>15%%)", p))
  # data <-> estimand consistency (sample truth uses realized target X)
  chk(abs(mean(d$Y_1[t_idx]) - d$mu1_true) < 0.05 &&
      abs(mean(d$Y_0[t_idx]) - d$mu0_true) < 0.05,
      sprintf("p=%d sample truth matches generated potential outcomes", p))
}

message("=== (2) Binary superpopulation truth matches generated data (p=50) ===")
set.seed(2)
ds <- generate_face_data(n_total = 20000L, K = 2L, p = 50L, config = "C1",
                         estimand_type = "superpopulation", outcome_type = "binary")
ti <- ds$R == "t"
message(sprintf("  superpop mu1=%.3f mu0=%.3f | realized mu1=%.3f mu0=%.3f",
                ds$mu1_true, ds$mu0_true, mean(ds$Y_1[ti]), mean(ds$Y_0[ti])))
chk(abs(mean(ds$Y_1[ti]) - ds$mu1_true) < 0.05 &&
    abs(mean(ds$Y_0[ti]) - ds$mu0_true) < 0.05,
    "superpop truth matches generated data")

message("=== (3) Continuous path untouched: estimand exactly Δ_T, deterministic ===")
set.seed(3)
c1 <- generate_face_data(2000L, K = 2L, p = 10L, config = "C1",
                         outcome_type = "continuous", estimand_type = "sample")
set.seed(3)
c2 <- generate_face_data(2000L, K = 2L, p = 10L, config = "C1",
                         outcome_type = "continuous", estimand_type = "sample")
chk(abs((c1$mu1_true - c1$mu0_true) - FACE_ATE_TARGET) < 1e-8,
    "continuous sample estimand == FACE_ATE_TARGET (3.0)")
chk(identical(c1$Y, c2$Y), "continuous outcomes deterministic (binary path inert)")
csup <- generate_face_data(2000L, K = 2L, p = 10L, config = "C1",
                           outcome_type = "continuous",
                           estimand_type = "superpopulation")
chk(abs((csup$mu1_true - csup$mu0_true) - FACE_ATE_TARGET) < 1e-8,
    "continuous superpop estimand == FACE_ATE_TARGET (3.0)")

message("=== (4) Deviation grid on the log-odds scale {0.5,1.0,1.5}: source not re-saturated ===")
for (dev in c(0.5, 1.0, 1.5)) {
  set.seed(4)
  dd <- generate_face_data(n_total = 9000L, K = 2L, p = 50L, config = "C1",
                           estimand_type = "sample", outcome_type = "binary",
                           ate_deviation = dev, n_deviated_sites = 2L)
  s1 <- dd$R == "s1"; s1tr <- which(dd$A == 1L & s1)
  if (length(s1tr) > 0) {
    n1 <- sum(dd$Y[s1tr] == 1L); n0 <- sum(dd$Y[s1tr] == 0L)
    message(sprintf("  dev=%.1f | deviated-source s1 treated minority=%d/%d (%.1f%%) mu1=%.3f",
                    dev, min(n1, n0), length(s1tr),
                    100 * min(n1, n0) / length(s1tr), mean(dd$Y_1[s1])))
    chk(min(n1, n0) / length(s1tr) > 0.05,
        sprintf("dev=%.1f deviated-source arm not saturated (>5%%)", dev))
  }
}

message(sprintf("=== OVERALL: %s ===", if (ok) "ALL PASS" else "SOME FAILED"))
if (!ok) quit(status = 1L)
