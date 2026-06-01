#!/usr/bin/env Rscript
# C2 lambda-grid comparison: FACE-HD dense-grid lambda.min  vs  RCAL x0.5+tune.cut.
#
# Diagnosis-only. Quantifies how differently TWO lambda-selection STRATEGIES pick the
# density-ratio penalty lambda_gamma, evaluated on the IDENTICAL out-of-fold calibrated-
# gamma CV entropy curve produced during a REAL C2 run_crossfit run:
#   (A) FACE-HD : lambda.min = argmin over its dense 100-point log grid (down to 1e-4*lmax when n>p).
#   (B) RCAL   : lambda.min on the coarse grid  lmax * 0.5^(0:10)  (nrho=11, tune.fac=0.5),
#                loss="cal" (entropy), tune.cut=TRUE  -- Tan's RCAL::glm.regu.cv default,
#                also what SMMAL used.
# Both minimize the SAME entropy curve; only grid resolution + floor (+ tune.cut) differ.
#
# Mechanism: the coverage-relevant lambda_gamma is chosen inside process_source_site()
# (R/cross_fitting_algorithms.R:208) by a fold-summed block-design call to
# fit_unified_density_ratio(..., calibrated=TRUE), which calls
# select_lambda_cv_calibrated_density_ratio_cpp() and picks lambda.min. We do NOT
# reconstruct that block design; instead we install a runtime namespace patch on that
# C++ wrapper to CAPTURE every (lambda_grid, cv_scores) curve from the live run, then
# post-process each curve two ways. NO package source files are modified.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(sprintf(
    "expected 9 args: tag n_total K p seed n_folds nlambda_init n_cores modes; got %d",
    length(args)
  ), call. = FALSE)
}

tag <- args[[1L]]
n_total <- as.integer(args[[2L]])
K_sites <- as.integer(args[[3L]])
p <- as.integer(args[[4L]])
seed <- as.integer(args[[5L]])
n_folds <- as.integer(args[[6L]])
nlambda_init <- as.integer(args[[7L]])
# n_cores (arg 8) is intentionally FORCED to 1 below: the namespace patch records into a
# parent-process environment, and under mclapply fork the child writes would be lost.
modes <- strsplit(args[[9L]], ",", fixed = TRUE)[[1L]]
modes <- modes[nzchar(modes)]
modes <- modes[modes %in% "one_round"]
if (length(modes) == 0L) modes <- "one_round"

stopifnot(n_total > 0L, K_sites > 0L, p > 0L, seed > 0L,
          n_folds >= 3L, nlambda_init > 0L)

suppressPackageStartupMessages({
  if (requireNamespace("FACEHD", quietly = FALSE)) {
    library(FACEHD)
  } else if (requireNamespace("devtools", quietly = FALSE)) {
    devtools::load_all(".")
  } else {
    stop("Neither installed FACEHD nor devtools available.", call. = FALSE)
  }
})

facehd_constant <- function(name, default) {
  ns <- asNamespace("FACEHD")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  default
}

facehd_function <- function(name) {
  ns <- asNamespace("FACEHD")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("FACEHD namespace does not contain required function '%s'.", name),
         call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}

needed_functions <- c(
  "generate_simulation_data", "split_data_by_site", "build_crossfit_folds",
  "run_crossfit", "fit_unified_density_ratio_cpp",
  "select_lambda_cv_calibrated_density_ratio_cpp"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, facehd_function(fn_name))
  }
}

MAX_ITER_DEFAULT <- facehd_constant("MAX_ITER_DEFAULT", 1000L)
TOL_DEFAULT <- facehd_constant("TOL_DEFAULT", 1e-7)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Sparsity (nnz) refit is OPT-IN: it adds two cold-start fixed-lambda gamma fits per captured
# curve, which at the tiny FACE-HD/RCAL lambda values run to MAX_ITER and dominate wall-clock.
# The core deliverable (lambda_FACEHD vs lambda_RCAL on the captured entropy curve) needs no refit.
NNZ_REFIT <- tolower(Sys.getenv("C2LG_NNZ", "false")) %in% c("1", "true", "yes")

out_root <- Sys.getenv(
  "C2_LAMBDA_GRID_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "lambda_grid_compare")
)
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d_kf%d", tag, n_total, K_sites, p, seed, n_folds
))

# ---------------------------------------------------------------------------
# RCAL emulation on a captured FACE-HD CV curve.
# RCAL::glm.regu.cv: rho.seq = tune.fac^seq(nrho-1,0,-1) * rho.max0 ; loss="cal";
#   sel <- which.min(err.ave) (lambda.min). tune.fac=0.5, nrho=1+10=11 (Tan vignette / SMMAL).
# tune.cut drops non-converged (small) rho; on FACE-HD's curve every grid point is finite
# (the C++ throws otherwise), and RCAL's floor (lmax*0.5^10 ~= 9.8e-4*lmax) sits ABOVE
# FACE-HD's 1e-4*lmax floor, so no rung is dropped -- the divergence is grid resolution + floor.
# ---------------------------------------------------------------------------
RCAL_NRHO <- 11L
RCAL_TUNE_FAC <- 0.5

rcal_emulate <- function(LG, CV) {
  LG <- as.numeric(LG); CV <- as.numeric(CV)
  ok <- is.finite(LG) & is.finite(CV) & LG > 0
  LG <- LG[ok]; CV <- CV[ok]
  if (length(LG) < 2L) return(NULL)
  lambda_max <- max(LG)
  rcal_grid <- lambda_max * RCAL_TUNE_FAC^(0:(RCAL_NRHO - 1L))   # 11 rungs, decreasing
  k_seq <- 0:(RCAL_NRHO - 1L)
  inrange <- rcal_grid >= min(LG) & rcal_grid <= max(LG)
  rg <- rcal_grid[inrange]; kk <- k_seq[inrange]
  if (length(rg) == 0L) return(NULL)
  logLG <- log(LG); ord <- order(logLG)
  ent <- stats::approx(logLG[ord], CV[ord], xout = log(rg))$y   # log-linear interp on entropy curve
  jmin <- which.min(ent)
  list(
    lambda_rcal_min = rg[[jmin]],
    rcal_rung_k = kk[[jmin]],
    entropy_rcal_min = ent[[jmin]]
  )
}

# Refit calibrated gamma at a FIXED lambda and count non-zero SLOPE coefficients
# (the intercept, element 1, is unpenalised and excluded from the sparsity count).
refit_nnz <- function(Z, A, mean_grad_psi, alpha_init, lambda, M_tau, W, A_val, fi, li) {
  tryCatch({
    fit <- fit_unified_density_ratio_cpp(
      Z, A, mean_grad_psi, alpha_init, lambda,
      MAX_ITER_DEFAULT, TOL_DEFAULT, TRUE, M_tau, W, A_val, fi, li, numeric(0)
    )
    g <- as.numeric(fit$gamma)
    if (length(g) > 1L) sum(abs(g[-1L]) > 1e-8) else sum(abs(g) > 1e-8)
  }, error = function(e) NA_integer_)
}

# ---------------------------------------------------------------------------
# Install the namespace patch that captures the calibrated-gamma CV curve.
# ---------------------------------------------------------------------------
CAP <- new.env(parent = emptyenv())
CAP$rec <- list()
CAP$n_calls <- 0L
CAP$t_cv <- 0     # cumulative seconds inside the original CV
CAP$t_post <- 0   # cumulative seconds inside capture post-processing (emulate + optional nnz)

ns <- asNamespace("FACEHD")
orig_cv <- get("select_lambda_cv_calibrated_density_ratio_cpp", envir = ns)

# Formal order MUST match R/RcppExports.R exactly.
patched_cv <- function(Z_site, A_source, mean_grad_psi, alpha_init, lambda_grid, n_folds,
                       max_iter, tol, M_tau, W_outcome, A_val, family_int, link_int) {
  t0 <- proc.time()[["elapsed"]]
  res <- orig_cv(Z_site, A_source, mean_grad_psi, alpha_init, lambda_grid, n_folds,
                 max_iter, tol, M_tau, W_outcome, A_val, family_int, link_int)
  CAP$n_calls <- CAP$n_calls + 1L
  CAP$t_cv <- CAP$t_cv + (proc.time()[["elapsed"]] - t0)
  cat(sprintf("[patch] CV call %d  M_tau=%.3g  n=%d p=%d  cv_secs=%.1f\n",
              CAP$n_calls, M_tau, nrow(Z_site), ncol(Z_site),
              proc.time()[["elapsed"]] - t0))
  # Only finite-M_tau calls are the calibrated-gamma curves that drive coverage.
  # (fit_unified_density_ratio reuses this kernel with M_tau=Inf for the refined branch.)
  if (is.finite(M_tau)) {
    tp <- proc.time()[["elapsed"]]
    LG <- as.numeric(lambda_grid)
    CV <- as.numeric(res$cv_scores)
    idx_min <- as.integer(res$idx_min)
    rec <- list(
      n = nrow(Z_site), p = ncol(Z_site), M_tau = as.numeric(M_tau),
      lambda_max = max(LG),
      lambda_facehd_min = as.numeric(res$lambda_min),
      idx_facehd_min = idx_min,
      lambda_facehd_1se = as.numeric(res$lambda_1se),
      entropy_facehd_min = if (idx_min >= 1L && idx_min <= length(CV)) CV[[idx_min]] else min(CV),
      facehd_at_floor = as.integer(idx_min == length(LG))
    )
    em <- rcal_emulate(LG, CV)
    if (is.null(em)) {
      rec$lambda_rcal_min <- NA_real_; rec$rcal_rung_k <- NA_integer_
      rec$entropy_rcal_min <- NA_real_
      rec$nnz_facehd <- NA_integer_; rec$nnz_rcal <- NA_integer_
    } else {
      rec$lambda_rcal_min <- em$lambda_rcal_min
      rec$rcal_rung_k <- em$rcal_rung_k
      rec$entropy_rcal_min <- em$entropy_rcal_min
      if (NNZ_REFIT) {
        rec$nnz_facehd <- refit_nnz(Z_site, A_source, mean_grad_psi, alpha_init,
                                   rec$lambda_facehd_min, M_tau, W_outcome,
                                   A_val, family_int, link_int)
        rec$nnz_rcal <- refit_nnz(Z_site, A_source, mean_grad_psi, alpha_init,
                                  rec$lambda_rcal_min, M_tau, W_outcome,
                                  A_val, family_int, link_int)
      } else {
        rec$nnz_facehd <- NA_integer_
        rec$nnz_rcal <- NA_integer_
      }
    }
    CAP$t_post <- CAP$t_post + (proc.time()[["elapsed"]] - tp)
    CAP$rec[[length(CAP$rec) + 1L]] <- rec
  }
  res
}

unlockBinding("select_lambda_cv_calibrated_density_ratio_cpp", ns)
assign("select_lambda_cv_calibrated_density_ratio_cpp", patched_cv, envir = ns)
restore_patch <- function() {
  assign("select_lambda_cv_calibrated_density_ratio_cpp", orig_cv, envir = ns)
  lockBinding("select_lambda_cv_calibrated_density_ratio_cpp", ns)
}

# ---------------------------------------------------------------------------
# Generate C2 data and run the real cross-fit (mirrors c2_score_moment_audit.R).
# ---------------------------------------------------------------------------
set.seed(seed)
data <- generate_simulation_data(
  n_total = n_total,
  K = K_sites,
  p = p,
  config = "C2",
  estimand_type = "superpopulation",
  site_allocation = "model",
  transform_type = "mild",
  outcome_type = "binary",
  heterogeneity_type = "none",
  shift_strength = 0.5,
  dgp_type = "facehd",
  warn_ignored = FALSE
)
split <- split_data_by_site(data)
folds <- build_crossfit_folds(split, n_folds)

t_run0 <- proc.time()[["elapsed"]]
run_ok <- tryCatch({
  for (mode in modes) {
    set.seed(seed + 200000L)   # one_round seed offset, matching the audit probe
    target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
    run_crossfit(
      split,
      n_folds = n_folds,
      communication_mode = mode,
      lambda_selection = "cv",
      lambda_rule = "min",
      verbose = FALSE,
      n_cores = 1L,                 # FORCED sequential so the capture is not lost to fork
      nlambda_init = nlambda_init,
      family = "binomial",
      A_val = 1L,
      use_lambda_cache = TRUE,
      precomputed_folds = folds,
      target_only_ps_cache = target_only_ps_cache,
      nuisance_lambda_rule = "min"
    )
  }
  TRUE
}, finally = restore_patch())
cat(sprintf("[timing] run_crossfit=%.1fs  captured_calls=%d  in_CV=%.1fs  in_post=%.1fs  (post incl. nnz=%s)\n",
            proc.time()[["elapsed"]] - t_run0, CAP$n_calls, CAP$t_cv, CAP$t_post,
            if (NNZ_REFIT) "TRUE" else "FALSE"))

records <- CAP$rec
if (length(records) == 0L) {
  stop("c2_lambda_grid_compare: captured 0 calibrated-gamma CV curves; the namespace patch did not intercept the calibrated path.",
       call. = FALSE)
}

# ---- per-call table ----
per_call <- do.call(rbind, lapply(seq_along(records), function(i) {
  r <- records[[i]]
  ratio <- r$lambda_facehd_min / r$lambda_rcal_min
  data.frame(
    tag = tag, seed = seed, n_total = n_total, K = K_sites, p = p,
    call_idx = i, n_obs = r$n, p_dim = r$p, M_tau = r$M_tau,
    lambda_max = r$lambda_max,
    lambda_facehd_min = r$lambda_facehd_min,
    idx_facehd_min = r$idx_facehd_min,
    lambda_facehd_1se = r$lambda_facehd_1se,
    lambda_rcal_min = r$lambda_rcal_min,
    rcal_rung_k = r$rcal_rung_k,
    ratio_facehd_over_rcal = ratio,
    log2_ratio = log2(ratio),
    entropy_facehd_min = r$entropy_facehd_min,
    entropy_rcal_min = r$entropy_rcal_min,
    entropy_gap = r$entropy_rcal_min - r$entropy_facehd_min,
    nnz_facehd = r$nnz_facehd,
    nnz_rcal = r$nnz_rcal,
    facehd_at_floor = r$facehd_at_floor,
    stringsAsFactors = FALSE
  )
}))

per_call_file <- paste0(prefix, "_per_call.csv")
write.csv(per_call, per_call_file, row.names = FALSE)

# ---- per-seed summary ----
fin <- function(x) x[is.finite(x)]
summary_row <- data.frame(
  tag = tag, seed = seed, n_total = n_total, K = K_sites, p = p,
  n_calls = nrow(per_call),
  median_ratio = stats::median(fin(per_call$ratio_facehd_over_rcal)),
  mean_ratio = mean(fin(per_call$ratio_facehd_over_rcal)),
  median_log2_ratio = stats::median(fin(per_call$log2_ratio)),
  mean_log2_ratio = mean(fin(per_call$log2_ratio)),
  mean_entropy_gap = mean(fin(per_call$entropy_gap)),
  mean_nnz_facehd = mean(fin(per_call$nnz_facehd)),
  mean_nnz_rcal = mean(fin(per_call$nnz_rcal)),
  frac_facehd_at_floor = mean(per_call$facehd_at_floor),
  stringsAsFactors = FALSE
)
summary_file <- paste0(prefix, "_summary.csv")
write.csv(summary_row, summary_file, row.names = FALSE)

cat(sprintf("[c2_lambda_grid_compare] seed=%d  n_calls=%d  median(lambda_FACEHD/lambda_RCAL)=%.4f  median(log2)=%.3f  nnz_FACEHD=%.1f  nnz_RCAL=%.1f  at_floor=%.2f\n",
            seed, summary_row$n_calls, summary_row$median_ratio,
            summary_row$median_log2_ratio, summary_row$mean_nnz_facehd,
            summary_row$mean_nnz_rcal, summary_row$frac_facehd_at_floor))
cat(sprintf("[c2_lambda_grid_compare] wrote %s and %s\n", per_call_file, summary_file))
