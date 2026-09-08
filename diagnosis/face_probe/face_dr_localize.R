#!/usr/bin/env Rscript
# Decisively localize the p=50 initial-DR CV bottleneck (observed 1259s for the
# first fit, 0.3s when lambda is cached). Hypothesis: the lambda path descends to
# 1e-4*lambda_max (LOW_DIM ratio), and the lightly-regularized fits on the
# collinear [X, X^2] design fail to converge within max_iter=10000, so each grinds
# the full 10000 iterations. We (A) reproduce the full-grid CV, (B) time single
# fixed-lambda fits across the path to find the slow region, and (C) test two
# candidate fixes (higher lambda_min_ratio; capped CV max_iter). Diagnostic only.
suppressPackageStartupMessages(library(devtools)); load_all(".", quiet = TRUE)

run_one <- function(p, n_total = 3000L) {
  set.seed(1)
  d <- generate_face_data(n_total = n_total, K = 2L, p = p, config = "C1",
                          outcome_type = "continuous", estimand_type = "sample")
  ti <- d$R == "t"; si <- d$R == "s1"
  Zt <- d$Z_site[ti, , drop = FALSE]
  Zs <- d$Z_site[si, , drop = FALSE]
  As <- d$A[si]
  mean_phi <- c(1, colMeans(Zt))            # intercept-augmented target moment
  collin <- max(abs(diag(cor(Zs[, 1:p, drop = FALSE],
                             Zs[, (p + 1):(2 * p), drop = FALSE]))))
  message(sprintf("=== p=%d | source Z=%dx%d treated(A=1)=%d | max|cor(X_j,X_j^2)|=%.2f ===",
                  p, nrow(Zs), ncol(Zs), sum(As == 1L), collin))
  lmax <- compute_lambda_max_initial_dr(Zs, As, mean_phi, 1L)
  nf   <- .nuisance_cv_fold_count(As, 1L, "probe")
  message(sprintf("  lambda_max=%.4g | cv_folds=%d | grid=%d | low-dim ratio=1e-4",
                  lmax, nf, LAMBDA_GRID_SIZE_STANDARD))

  # (A) full-grid CV at default settings (reproduce the slow path)
  t0 <- Sys.time()
  g <- fit_initial_density_ratio(Zs, As, mean_phi, lambda = NULL, A_val = 1L)
  message(sprintf("  [A] CV full grid (ratio=1e-4, maxit=10000): %.1fs  selected lambda=%.4g",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  attr(g, "lambda_min")))

  # (B) single fixed-lambda fits across the path -> which region is slow?
  for (r in c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {
    lam <- lmax * r
    t0 <- Sys.time()
    gg <- fit_initial_density_ratio(Zs, As, mean_phi, lambda = lam, A_val = 1L)
    message(sprintf("  [B] single lambda=%.4g (r=%.0e): %.2fs  ||g||1=%.2f nnz=%d",
                    lam, r, as.numeric(difftime(Sys.time(), t0, units = "secs")),
                    sum(abs(gg)), sum(abs(as.numeric(gg)) > 1e-8)))
  }

  # (C1) candidate fix: higher lambda_min_ratio (don't descend to near-zero penalty)
  for (lmr in c(0.05, 0.01)) {
    grid <- build_lambda_grid(lmax, lambda_min_ratio = lmr)
    t0 <- Sys.time()
    invisible(select_lambda_cv_initial_density_ratio_cpp(
      Zs, As, mean_phi, grid, nf, 10000L, 1e-6, 1L))
    message(sprintf("  [C1] CV ratio=%.2f maxit=10000: %.1fs",
                    lmr, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  # (C2) candidate fix: capped CV max_iter (exact convergence not needed for selection)
  for (mi in c(300L, 1000L)) {
    grid <- build_lambda_grid(lmax, lambda_min_ratio = 1e-4)
    t0 <- Sys.time()
    invisible(select_lambda_cv_initial_density_ratio_cpp(
      Zs, As, mean_phi, grid, nf, mi, 1e-6, 1L))
    message(sprintf("  [C2] CV ratio=1e-4 maxit=%d: %.1fs",
                    mi, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
}

run_one(50L)
message("[dr-localize] DONE")
