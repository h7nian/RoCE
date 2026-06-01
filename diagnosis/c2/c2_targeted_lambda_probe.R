#!/usr/bin/env Rscript
# C2 targeted-lambda GATE probe: ENTROPY lambda.min vs BALANCING-MOMENT lambda.min.
#
# Diagnosis-only. On the IDENTICAL out-of-fold calibrated-gamma CV folds produced during a
# REAL C2 run_crossfit run, this probe compares TWO lambda-selection objectives evaluated on
# the SAME captured lambda_grid:
#   (A) lambda_entropy   = current FACE-HD pick = lambda.min of the out-of-fold density-ratio
#                          ENTROPY curve (res$lambda_min) -- what the estimator actually uses.
#   (B) lambda_balancing = argmin over the grid of the out-of-fold BALANCING-MOMENT MSE
#                          ( mean over our OWN CV folds of ||U_gamma||^2 ), where U_gamma is the
#                          first-order bias driver (target_mean - source_mean) authoritatively
#                          defined by fold_gamma_score() in c2_score_moment_audit.R:129-147.
#
# GATE LOGIC: if lambda_balancing ~= lambda_entropy (as FACE-HD-vs-RCAL turned out), a targeted
# lambda cannot move coverage and we stop here. If they differ materially (|log2 ratio| >= one
# full grid octave on a non-trivial fraction of calls), a separate end-to-end coverage run
# follows. This probe does NOT change the estimator's lambda; it only COMPUTES and COMPARES the
# two selections plus their out-of-fold balancing.
#
# Mechanism: same runtime namespace patch on select_lambda_cv_calibrated_density_ratio_cpp() as
# c2_lambda_grid_compare.R. We CAPTURE every (lambda_grid, cv_scores, lambda_min, idx_min) curve
# from the live run, then -- only for finite-M_tau (calibrated-gamma) calls -- run our OWN
# stratified K-fold balancing CV on the same source rows over the same grid. NO package source
# files are modified.

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
  "run_crossfit", "fit_unified_density_ratio_cpp", "predict_glm_cpp",
  "select_lambda_cv_calibrated_density_ratio_cpp"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, facehd_function(fn_name))
  }
}

MAX_ITER_DEFAULT <- facehd_constant("MAX_ITER_DEFAULT", 1000L)
TOL_DEFAULT <- facehd_constant("TOL_DEFAULT", 1e-7)
LOGISTIC_CLIP <- facehd_constant("LOGISTIC_CLIP", 50)
FAMILY_BINOMIAL <- facehd_constant("FAMILY_BINOMIAL", 1L)
LINK_LOGIT <- facehd_constant("LINK_LOGIT", 1L)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Number of OWN balancing CV folds (stratified by A). Kept small (3) so the warm-started path
# over the captured grid stays cheap relative to the entropy CV the live run already paid for.
K_BAL <- 3L
# Local RNG offset: the balancing-fold assignment uses a DEDICATED stream so run_crossfit's
# own .Random.seed is never perturbed (we save/restore around every captured call).
BAL_SEED_OFFSET <- 777000L

out_root <- Sys.getenv(
  "C2_TARGETED_LAMBDA_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "targeted_lambda_probe")
)
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d_kf%d", tag, n_total, K_sites, p, seed, n_folds
))

# ---------------------------------------------------------------------------
# Balancing-moment helpers (authoritative = fold_gamma_score, c2_score_moment_audit.R:129-147).
# psi' depends ONLY on alpha (NOT on gamma/lambda); only exp(-eta_gamma) changes per lambda.
# ---------------------------------------------------------------------------
clip_to <- function(eta, bound) pmax(pmin(as.numeric(eta), bound), -bound)

# source_mean over a set of rows: colMeans( w * Z_int ), w = I(A==A_val) exp(-eta_gamma) psi'.
# U_gamma = mean_grad_psi (the fixed aggregate target-side gradient passed to the patch) -
#           source_mean. We do NOT recompute the target side.
balancing_source_mean <- function(Z_int, A_rows, A_val, gamma, psi_prime) {
  eta_gamma <- clip_to(as.numeric(Z_int %*% gamma), LOGISTIC_CLIP)
  w <- as.numeric(A_rows == A_val) * exp(-eta_gamma) * psi_prime
  colMeans(sweep(Z_int, 1L, w, "*"))
}

# ---------------------------------------------------------------------------
# Install the namespace patch that captures the calibrated-gamma CV curve AND, for finite-M_tau
# calls, runs the OWN stratified balancing CV over the captured grid.
# ---------------------------------------------------------------------------
CAP <- new.env(parent = emptyenv())
CAP$rec <- list()
CAP$n_calls <- 0L
CAP$t_cv <- 0     # cumulative seconds inside the original CV
CAP$t_post <- 0   # cumulative seconds inside capture post-processing (balancing CV)

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
    idx_entropy <- as.integer(res$idx_min)
    n_grid <- length(LG)

    rec <- list(
      n = nrow(Z_site), p = ncol(Z_site), M_tau = as.numeric(M_tau),
      lambda_max = max(LG),
      lambda_entropy = as.numeric(res$lambda_min),
      idx_entropy = idx_entropy,
      lambda_balancing = NA_real_,
      idx_balancing = NA_integer_,
      bal_mse_at_entropy = NA_real_,
      bal_mse_at_balancing = NA_real_
    )

    # ---- OWN stratified balancing CV (wrapped so an error never aborts the live run) ----
    bal <- tryCatch({
      Z <- as.matrix(Z_site)
      A <- as.numeric(A_source)
      alpha <- as.numeric(alpha_init)
      target_mean <- as.numeric(mean_grad_psi)   # fixed target side; DO NOT recompute.
      n_rows <- nrow(Z)
      Z_int <- cbind(1, Z)                        # gamma length = ncol(Z)+1; dims must match.
      # W_outcome here is the plugin BLOCK design (ncol = n_folds*(p+1)); the C++ and
      # predict_glm_cpp prepend a GLOBAL intercept, and alpha = c(0, per-fold...) has length
      # 1 + n_folds*(p+1). So cbind(1, W_outcome) %*% alpha matches the C++ contract and the
      # authoritative fold_gamma_score (c2_score_moment_audit.R:140) / .mean_glm_gradient_site_basis.
      W_int <- cbind(1, W_outcome)

      # psi' computed ONCE for ALL rows (depends only on alpha, not gamma/lambda).
      eta_alpha <- as.numeric(W_int %*% alpha)
      # M_tau truncation on eta_alpha is INERT at C2 p50 (proven), but clip for fidelity.
      eta_alpha <- clip_to(eta_alpha, M_tau)
      mu <- 1 / (1 + exp(-clip_to(eta_alpha, LOGISTIC_CLIP)))
      psi_prime <- mu * (1 - mu)

      # Sanity: block matrix mult / intercept handling matches predict_glm_cpp.
      mu_chk <- as.numeric(predict_glm_cpp(W_outcome, alpha, family_int, link_int))
      mu_ref <- 1 / (1 + exp(-clip_to(as.numeric(W_int %*% alpha), LOGISTIC_CLIP)))
      if (max(abs(mu_chk - mu_ref)) > 1e-6) {
        stop(sprintf("predict_glm_cpp mismatch (max abs diff=%.3e): block matrix mult / intercept handling failed.",
                     max(abs(mu_chk - mu_ref))), call. = FALSE)
      }

      # Stratified-by-A fold assignment via a DEDICATED local RNG stream (save/restore the
      # global .Random.seed so run_crossfit's stream is untouched).
      had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      saved_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv, inherits = FALSE) else NULL
      on.exit({
        if (had_seed) {
          assign(".Random.seed", saved_seed, envir = .GlobalEnv)
        } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed + BAL_SEED_OFFSET + CAP$n_calls)

      assign_strat <- function(idx) {
        ni <- length(idx)
        if (ni == 0L) return(integer(0))
        # cyclic assignment after a shuffle: balances counts per fold within each A stratum,
        # guaranteeing every fold receives A==A_val rows when enough exist.
        ((sample.int(ni) - 1L) %% K_BAL) + 1L
      }
      fold_id <- integer(n_rows)
      is_val <- (A == A_val)
      idx_pos <- which(is_val)
      idx_neg <- which(!is_val)
      fold_id[idx_pos] <- assign_strat(idx_pos)
      fold_id[idx_neg] <- assign_strat(idx_neg)

      bal_sq <- numeric(n_grid)
      for (v in seq_len(K_BAL)) {
        val <- which(fold_id == v)
        train <- which(fold_id != v)
        if (length(train) == 0L || length(val) == 0L) next
        Z_int_val <- Z_int[val, , drop = FALSE]
        A_val_rows <- A[val]
        psi_val <- psi_prime[val]
        warm <- numeric(0)   # cold at i=1, then warm-started down the decreasing grid.
        for (i in seq_len(n_grid)) {
          fit <- fit_unified_density_ratio_cpp(
            Z[train, , drop = FALSE], A[train], target_mean, alpha, LG[[i]],
            MAX_ITER_DEFAULT, TOL_DEFAULT, TRUE, M_tau,
            W_outcome[train, , drop = FALSE], A_val, family_int, link_int, warm
          )
          warm <- as.numeric(fit$gamma)
          source_mean_v <- balancing_source_mean(Z_int_val, A_val_rows, A_val, warm, psi_val)
          U <- target_mean - source_mean_v
          bal_sq[[i]] <- bal_sq[[i]] + sum(U * U)
        }
      }
      bal_mse <- bal_sq / K_BAL
      idx_balancing <- which.min(bal_mse)
      list(
        bal_mse = bal_mse,
        idx_balancing = idx_balancing,
        lambda_balancing = LG[[idx_balancing]],
        bal_mse_at_balancing = bal_mse[[idx_balancing]],
        bal_mse_at_entropy = if (idx_entropy >= 1L && idx_entropy <= n_grid)
          bal_mse[[idx_entropy]] else NA_real_
      )
    }, error = function(e) {
      cat(sprintf("[patch] balancing CV failed on call %d: %s\n", CAP$n_calls,
                  conditionMessage(e)))
      NULL
    })

    if (!is.null(bal)) {
      rec$lambda_balancing <- bal$lambda_balancing
      rec$idx_balancing <- as.integer(bal$idx_balancing)
      rec$bal_mse_at_entropy <- bal$bal_mse_at_entropy
      rec$bal_mse_at_balancing <- bal$bal_mse_at_balancing
    }
    rec$n_grid <- n_grid
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
# Generate C2 data and run the real cross-fit (identical to c2_lambda_grid_compare.R).
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
cat(sprintf("[timing] run_crossfit=%.1fs  captured_calls=%d  in_CV=%.1fs  in_post(balancing)=%.1fs\n",
            proc.time()[["elapsed"]] - t_run0, CAP$n_calls, CAP$t_cv, CAP$t_post))

records <- CAP$rec
if (length(records) == 0L) {
  stop("c2_targeted_lambda_probe: captured 0 calibrated-gamma CV curves; the namespace patch did not intercept the calibrated path.",
       call. = FALSE)
}

# ---- per-call table ----
per_call <- do.call(rbind, lapply(seq_along(records), function(i) {
  r <- records[[i]]
  ratio <- r$lambda_entropy / r$lambda_balancing
  reduction <- if (is.finite(r$bal_mse_at_entropy) && is.finite(r$bal_mse_at_balancing) &&
                   r$bal_mse_at_entropy > 0) {
    1 - r$bal_mse_at_balancing / r$bal_mse_at_entropy
  } else NA_real_
  data.frame(
    tag = tag, seed = seed, n_total = n_total, K = K_sites, p = p,
    call_idx = i, n_obs = r$n, p_dim = r$p, M_tau = r$M_tau,
    lambda_max = r$lambda_max,
    lambda_entropy = r$lambda_entropy,
    idx_entropy = r$idx_entropy,
    lambda_balancing = r$lambda_balancing,
    idx_balancing = r$idx_balancing,
    ratio_ent_over_bal = ratio,
    log2_ratio = log2(ratio),
    bal_mse_at_entropy = r$bal_mse_at_entropy,
    bal_mse_at_balancing = r$bal_mse_at_balancing,
    bal_mse_reduction = reduction,
    entropy_at_floor = as.integer(r$idx_entropy == r$n_grid),
    balancing_at_floor = as.integer(!is.na(r$idx_balancing) && r$idx_balancing == r$n_grid),
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
  median_ratio_ent_over_bal = stats::median(fin(per_call$ratio_ent_over_bal)),
  mean_ratio_ent_over_bal = mean(fin(per_call$ratio_ent_over_bal)),
  median_log2_ratio = stats::median(fin(per_call$log2_ratio)),
  mean_log2_ratio = mean(fin(per_call$log2_ratio)),
  median_bal_mse_reduction = stats::median(fin(per_call$bal_mse_reduction)),
  frac_log2_ge_octave = mean(abs(fin(per_call$log2_ratio)) >= 1),
  frac_entropy_at_floor = mean(per_call$entropy_at_floor),
  frac_balancing_at_floor = mean(per_call$balancing_at_floor),
  stringsAsFactors = FALSE
)
summary_file <- paste0(prefix, "_summary.csv")
write.csv(summary_row, summary_file, row.names = FALSE)

cat(sprintf("[c2_targeted_lambda_probe] seed=%d  n_calls=%d  median(lambda_ent/lambda_bal)=%.4f  median(log2)=%.3f  median_bal_reduction=%.3f  frac|log2|>=1=%.2f  ent_floor=%.2f  bal_floor=%.2f\n",
            seed, summary_row$n_calls, summary_row$median_ratio_ent_over_bal,
            summary_row$median_log2_ratio, summary_row$median_bal_mse_reduction,
            summary_row$frac_log2_ge_octave, summary_row$frac_entropy_at_floor,
            summary_row$frac_balancing_at_floor))
cat(sprintf("[c2_targeted_lambda_probe] wrote %s and %s\n", per_call_file, summary_file))
