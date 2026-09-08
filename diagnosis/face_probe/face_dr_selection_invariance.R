#!/usr/bin/env Rscript
# Selection-invariance safety check for the proposed DR-CV speedup.
# The bottleneck is that the CV-selection pass runs up to CV_MAX_ITER=10000
# iterations even on never-selected small-lambda fits. The intended knob is the
# (currently inert) cv_max_iter = min(max_iter, CV_MAX_ITER) cap. This probe asks
# the only question that matters for correctness: does capping the CV-pass
# iterations (and/or raising lambda_min_ratio) change the SELECTED lambda and the
# resulting FINAL (full-max_iter) gamma? We compare against the current behavior
# (effective cap 10000, ratio 1e-4) on the FACE p=50 cell that is slow.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

set.seed(1)
p <- 50L
d <- generate_face_data(n_total = 3000L, K = 2L, p = p, config = "C1",
                        outcome_type = "continuous", estimand_type = "sample")
ti <- d$R == "t"; si <- d$R == "s1"
Zs <- d$Z_site[si, , drop = FALSE]; As <- d$A[si]
mean_phi <- c(1, colMeans(d$Z_site[ti, , drop = FALSE]))
lmax <- compute_lambda_max_initial_dr(Zs, As, mean_phi, 1L)
nf   <- .nuisance_cv_fold_count(As, 1L, "probe")

# Helper: select lambda via the CV kernel at a given (effective cv cap, ratio),
# then do the FINAL refit at full max_iter=10000, and return both.
select_and_refit <- function(cv_cap, ratio) {
  grid <- build_lambda_grid(lmax, lambda_min_ratio = ratio)
  t0 <- Sys.time()
  cv  <- select_lambda_cv_initial_density_ratio_cpp(Zs, As, mean_phi, grid, nf,
                                                    cv_cap, 1e-6, 1L)
  t_cv <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  lam  <- cv$lambda_min
  g    <- fit_initial_density_ratio_cpp(Zs, As, mean_phi, lam, 10000L, 1e-6, 1L, numeric(0))
  list(t_cv = t_cv, lambda = lam, gamma = as.numeric(g$gamma),
       conv = g$converged, iters = g$iterations)
}

ref <- select_and_refit(10000L, 1e-4)   # current behavior (baseline)
message(sprintf("[baseline] cap=10000 ratio=1e-4 | CV %.1fs | sel lambda=%.5g | refit conv=%s iters=%d",
                ref$t_cv, ref$lambda, ref$conv, ref$iters))

settings <- list(
  list(cap = 1000L, ratio = 1e-4), list(cap = 500L, ratio = 1e-4),
  list(cap = 300L,  ratio = 1e-4), list(cap = 1000L, ratio = 0.05),
  list(cap = 2000L, ratio = 0.02), list(cap = 10000L, ratio = 0.05)
)
for (s in settings) {
  r <- select_and_refit(s$cap, s$ratio)
  d_lam  <- abs(r$lambda - ref$lambda) / max(ref$lambda, 1e-12)
  d_gam  <- max(abs(r$gamma - ref$gamma))
  rel_g  <- d_gam / max(max(abs(ref$gamma)), 1e-12)
  message(sprintf(
    "  cap=%5d ratio=%.0e | CV %6.1fs | sel lambda=%.5g (Δrel=%.1e) | maxΔgamma=%.2e (rel=%.1e) | %s",
    s$cap, s$ratio, r$t_cv, r$lambda, d_lam, d_gam, rel_g,
    if (d_lam < 1e-6 && rel_g < 1e-3) "SELECTION+FIT IDENTICAL"
    else if (rel_g < 0.02) "negligible fit change" else "CHANGED"))
}
message("[dr-selection-invariance] DONE")
