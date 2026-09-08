#!/usr/bin/env Rscript
# Empirically test the density-ratio conditioning hypothesis: the initial-DR
# lambda-path CV is catastrophically slow at p=50 because the exp-tilting solver
# runs on the UNSTANDARDIZED [X, X^2] design (column scales spanning orders of
# magnitude). glmnet standardizes columns by default; our custom DR solver does
# not. Here we A/B the SAME fit on raw vs. column-scaled Z (scale-only
# reparameterization Z~ = Z/s, mean_phi~ = mean_phi/s, gamma = gamma~/s) and time
# both. A large speedup confirms the fix; finite/sane weights confirm correctness.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

run_one <- function(p) {
  set.seed(1)
  d <- generate_face_data(n_total = 3000L, K = 2L, p = p, config = "C1",
                          outcome_type = "continuous", estimand_type = "sample")
  ti <- d$R == "t"; si <- d$R == "s1"
  Zt <- d$Z_site[ti, , drop = FALSE]
  Zs <- d$Z_site[si, , drop = FALSE]
  As <- d$A[si]
  mean_phi <- colMeans(Zt)                 # target moment the initial DR matches
  s <- apply(Zs, 2, sd); s[s < 1e-8] <- 1
  message(sprintf("p=%d | source Z=%dx%d | col-sd range [%.3f, %.1f]",
                  p, nrow(Zs), ncol(Zs), min(s), max(s)))

  t0 <- Sys.time()
  g_raw <- tryCatch(fit_initial_density_ratio(Zs, As, mean_phi, lambda = NULL, A_val = 1L),
                    error = function(e) { message("  [raw] ERROR: ", conditionMessage(e)); NULL })
  t_raw <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  [raw]    init-DR CV in %.1fs", t_raw))

  Zs_t <- sweep(Zs, 2, s, "/"); mp_t <- mean_phi / s
  t0 <- Sys.time()
  g_til <- tryCatch(fit_initial_density_ratio(Zs_t, As, mp_t, lambda = NULL, A_val = 1L),
                    error = function(e) { message("  [scaled] ERROR: ", conditionMessage(e)); NULL })
  t_scl <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("  [scaled] init-DR CV in %.1fs  (speedup %.0fx)",
                  t_scl, t_raw / max(t_scl, 1e-6)))

  if (!is.null(g_raw) && !is.null(g_til)) {
    g_back <- as.numeric(g_til) / s
    wr <- exp(as.numeric(Zs %*% as.numeric(g_raw)))
    wb <- exp(as.numeric(Zs %*% g_back))
    message(sprintf("  weights: raw range [%.3g, %.3g] | scaled range [%.3g, %.3g] | corr=%.3f",
                    min(wr), max(wr), min(wb), max(wb), suppressWarnings(cor(wr, wb))))
  }
}

for (p in c(10L, 50L)) run_one(p)
message("[dr-cond] DONE")
