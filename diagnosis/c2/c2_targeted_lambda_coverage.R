#!/usr/bin/env Rscript
# C2 targeted-lambda COVERAGE: baseline (entropy lambda.min) vs targeted (balancing-MSE argmin).
#
# Diagnosis-only END-TO-END coverage measurement. Companion to (and downstream of) the GATE probe
# diagnosis/c2/c2_targeted_lambda_probe.R, which established that on the IDENTICAL out-of-fold
# calibrated-gamma CV folds the entropy lambda.min and the out-of-fold BALANCING-MOMENT MSE argmin
# can differ. THIS script measures whether forcing the estimator to use lambda_balancing changes
# the FACE-HD estimator's coverage, PAIRED on the same data / folds / seed.
#
#   ARM baseline : the estimator as-is. The namespace patch on
#                  select_lambda_cv_calibrated_density_ratio_cpp() is CAPTURE-ONLY -- it returns the
#                  original entropy-lambda.min CV result UNCHANGED and does NOT run any balancing CV.
#   ARM targeted : the SAME estimator, but the patch -- for each finite-M_tau (calibrated-gamma)
#                  CV call -- runs the OWN out-of-fold stratified balancing CV (the VERBATIM block
#                  from c2_targeted_lambda_probe.R) over a COARSER subgrid of the captured
#                  lambda_grid, finds lambda_balancing = grid[which.min(bal_mse)], and OVERRIDES
#                  the CV result's lambda_min / lambda_1se (+ best_lambda / idx_min / best_idx) to
#                  lambda_balancing so the downstream .select_nuisance_cv_lambda(res,"min") in
#                  fit_unified_density_ratio fits gamma at lambda_balancing.
#
# OVERRIDE SUFFICIENCY (verified against R/model_fitting.R:102-126): for lambda_rule="min",
# .select_nuisance_cv_lambda reads cv_result$lambda_min (line 120-121) and REQUIRES that BOTH
# lambda_min AND lambda_1se be present, length-1, and finite (line 108-119). No other field is
# consulted for the selected lambda value. We therefore set lambda_min <- lambda_1se <-
# best_lambda <- lambda_balancing and idx_min <- best_idx <- idx_balancing for full consistency.
#
# The balancing-CV math (psi' precompute from eta_alpha = cbind(1,W_outcome) %*% alpha, the
# stratified-by-A K_bal folds with a dedicated save/restore RNG stream, the warm-started gamma path
# via fit_unified_density_ratio_cpp passing RAW W_outcome[train,], the balancing_source_mean
# helper, U = target_mean - source_mean, bal_mse) is REUSED VERBATIM from the validated probe and
# is NOT re-derived here. NO package source files are modified.

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
mode <- modes[[1L]]   # coverage is measured for a single communication_mode (one_round).

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

# Number of OWN balancing CV folds (stratified by A). Kept small (default 3) so the warm-started
# path over the (subsampled) captured grid stays cheap relative to the entropy CV the live run
# already paid for. Only the targeted arm pays this cost.
K_BAL <- as.integer(Sys.getenv("C2TLC_K_BAL", "3"))
if (!is.finite(K_BAL) || K_BAL < 2L) {
  stop(sprintf("C2TLC_K_BAL must be an integer >= 2; got '%s'.",
               Sys.getenv("C2TLC_K_BAL")), call. = FALSE)
}
# Coarser balancing search grid: subsample the captured lambda_grid with this stride to cut the
# ~9000-fit full-grid cost. idx_entropy is ALWAYS forced into the subgrid (see below) so the
# targeted arm can never do worse than baseline on the entropy point.
BAL_STRIDE <- as.integer(Sys.getenv("C2TLC_BAL_STRIDE", "4"))
if (!is.finite(BAL_STRIDE) || BAL_STRIDE < 1L) {
  stop(sprintf("C2TLC_BAL_STRIDE must be an integer >= 1; got '%s'.",
               Sys.getenv("C2TLC_BAL_STRIDE")), call. = FALSE)
}
# Local RNG offset: the balancing-fold assignment uses a DEDICATED stream so run_crossfit's
# own .Random.seed is never perturbed (we save/restore around every captured call).
BAL_SEED_OFFSET <- 777000L

out_root <- Sys.getenv(
  "C2TLC_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "targeted_lambda_coverage")
)
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d_kf%d_arms.csv", tag, n_total, K_sites, p, seed, n_folds
))

# ---------------------------------------------------------------------------
# Balancing-moment helpers (authoritative = fold_gamma_score, c2_score_moment_audit.R:129-147).
# psi' depends ONLY on alpha (NOT on gamma/lambda); only exp(-eta_gamma) changes per lambda.
# REUSED VERBATIM from c2_targeted_lambda_probe.R.
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

ns <- asNamespace("FACEHD")
orig_cv <- get("select_lambda_cv_calibrated_density_ratio_cpp", envir = ns)

# ---------------------------------------------------------------------------
# Per-arm patch factory. OVERRIDE=FALSE -> capture-only baseline (NO balancing CV; returns the
# original entropy CV result UNCHANGED). OVERRIDE=TRUE -> compute the out-of-fold balancing CV on
# the coarser subgrid and OVERRIDE the CV result so gamma is fit at lambda_balancing.
# CAP records per-call lambda_entropy / lambda_balancing so we can report the ratio.
# ---------------------------------------------------------------------------
make_patched_cv <- function(OVERRIDE, CAP) {
  function(Z_site, A_source, mean_grad_psi, alpha_init, lambda_grid, n_folds,
           max_iter, tol, M_tau, W_outcome, A_val, family_int, link_int) {
    t0 <- proc.time()[["elapsed"]]
    res <- orig_cv(Z_site, A_source, mean_grad_psi, alpha_init, lambda_grid, n_folds,
                   max_iter, tol, M_tau, W_outcome, A_val, family_int, link_int)
    CAP$n_calls <- CAP$n_calls + 1L
    CAP$t_cv <- CAP$t_cv + (proc.time()[["elapsed"]] - t0)

    # Only finite-M_tau calls are the calibrated-gamma curves that drive coverage.
    # (fit_unified_density_ratio reuses this kernel with M_tau=Inf for the refined branch; we leave
    #  those UNCHANGED in both arms so only the calibrated-gamma lambda differs.)
    if (!is.finite(M_tau)) {
      return(res)
    }

    LG <- as.numeric(lambda_grid)
    idx_entropy <- as.integer(res$idx_min)
    lambda_entropy <- as.numeric(res$lambda_min)
    n_grid <- length(LG)

    # In BASELINE mode SKIP the balancing CV entirely (keeps baseline fast); record entropy only.
    if (!OVERRIDE) {
      cat(sprintf("[patch:base] CV call %d  M_tau=%.3g  n=%d p=%d  cv_secs=%.1f  lambda_entropy=%.4g\n",
                  CAP$n_calls, M_tau, nrow(Z_site), ncol(Z_site),
                  proc.time()[["elapsed"]] - t0, lambda_entropy))
      CAP$rec[[length(CAP$rec) + 1L]] <- list(
        lambda_entropy = lambda_entropy, lambda_balancing = NA_real_,
        overridden = FALSE
      )
      return(res)   # UNCHANGED entropy lambda.min.
    }

    # ---- TARGETED mode: OWN stratified balancing CV over a COARSER subgrid. ----
    tp <- proc.time()[["elapsed"]]
    cat(sprintf("[patch:tgt]  CV call %d  M_tau=%.3g  n=%d p=%d  cv_secs=%.1f  lambda_entropy=%.4g\n",
                CAP$n_calls, M_tau, nrow(Z_site), ncol(Z_site),
                proc.time()[["elapsed"]] - t0, lambda_entropy))

    # COARSER balancing search subgrid: stride BAL_STRIDE over the full captured grid, ALWAYS
    # including idx_entropy so the targeted arm can never do worse than baseline on the entropy
    # point. sub_idx maps positions in the subgrid back to the full grid.
    sub_idx <- sort(unique(c(seq.int(1L, n_grid, by = BAL_STRIDE), idx_entropy)))
    LG_sub <- LG[sub_idx]
    n_sub <- length(LG_sub)
    # position of the entropy lambda within the subgrid (guaranteed present).
    pos_entropy_sub <- match(idx_entropy, sub_idx)

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

      # bal_sq accumulated over the SUBGRID only (warm-started down the decreasing subgrid).
      bal_sq <- numeric(n_sub)
      for (v in seq_len(K_BAL)) {
        val <- which(fold_id == v)
        train <- which(fold_id != v)
        if (length(train) == 0L || length(val) == 0L) next
        Z_int_val <- Z_int[val, , drop = FALSE]
        A_val_rows <- A[val]
        psi_val <- psi_prime[val]
        warm <- numeric(0)   # cold at i=1, then warm-started down the decreasing subgrid.
        for (i in seq_len(n_sub)) {
          fit <- fit_unified_density_ratio_cpp(
            Z[train, , drop = FALSE], A[train], target_mean, alpha, LG_sub[[i]],
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
      pos_balancing_sub <- which.min(bal_mse)
      idx_balancing <- sub_idx[[pos_balancing_sub]]   # map back to the FULL grid index.
      list(
        idx_balancing = idx_balancing,
        lambda_balancing = LG[[idx_balancing]],
        bal_mse_at_balancing = bal_mse[[pos_balancing_sub]],
        bal_mse_at_entropy = if (!is.na(pos_entropy_sub)) bal_mse[[pos_entropy_sub]] else NA_real_
      )
    }, error = function(e) {
      cat(sprintf("[patch:tgt]  balancing CV failed on call %d: %s\n", CAP$n_calls,
                  conditionMessage(e)))
      NULL
    })

    CAP$t_post <- CAP$t_post + (proc.time()[["elapsed"]] - tp)

    if (is.null(bal)) {
      # Balancing CV failed: leave the entropy CV result UNCHANGED (degrade to baseline behaviour
      # on this call) rather than aborting the live run.
      CAP$rec[[length(CAP$rec) + 1L]] <- list(
        lambda_entropy = lambda_entropy, lambda_balancing = NA_real_, overridden = FALSE
      )
      return(res)
    }

    lambda_balancing <- as.numeric(bal$lambda_balancing)
    idx_balancing <- as.integer(bal$idx_balancing)
    cat(sprintf("[patch:tgt]  call %d  lambda_entropy=%.4g -> lambda_balancing=%.4g  ratio=%.3f  bal_secs=%.1f\n",
                CAP$n_calls, lambda_entropy, lambda_balancing,
                lambda_entropy / lambda_balancing,
                proc.time()[["elapsed"]] - tp))
    CAP$rec[[length(CAP$rec) + 1L]] <- list(
      lambda_entropy = lambda_entropy, lambda_balancing = lambda_balancing,
      overridden = TRUE
    )

    # ---- OVERRIDE the CV result so .select_nuisance_cv_lambda(res,"min") picks lambda_balancing.
    # Verified against R/model_fitting.R:102-126: rule="min" reads res$lambda_min and REQUIRES
    # finite length-1 lambda_min AND lambda_1se. Set all selection fields consistently.
    res$lambda_min <- lambda_balancing
    res$lambda_1se <- lambda_balancing
    res$best_lambda <- lambda_balancing
    res$idx_min <- idx_balancing
    res$best_idx <- idx_balancing
    res
  }
}

install_patch <- function(patched) {
  unlockBinding("select_lambda_cv_calibrated_density_ratio_cpp", ns)
  assign("select_lambda_cv_calibrated_density_ratio_cpp", patched, envir = ns)
}
restore_patch <- function() {
  if (bindingIsLocked("select_lambda_cv_calibrated_density_ratio_cpp", ns)) {
    unlockBinding("select_lambda_cv_calibrated_density_ratio_cpp", ns)
  }
  assign("select_lambda_cv_calibrated_density_ratio_cpp", orig_cv, envir = ns)
  lockBinding("select_lambda_cv_calibrated_density_ratio_cpp", ns)
}

# ---------------------------------------------------------------------------
# Generate C2 data and build folds ONCE (identical call to the probe). Both arms share these.
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
truth <- as.numeric(data$mu1_true)

# ---------------------------------------------------------------------------
# run_arm: install the per-OVERRIDE patch, run the SAME run_crossfit on identical data/folds/seed,
# restore the patch on exit. The pre-run seed is IDENTICAL across arms so ONLY the calibrated-gamma
# lambda differs between baseline and targeted.
# ---------------------------------------------------------------------------
run_arm <- function(OVERRIDE) {
  CAP <- new.env(parent = emptyenv())
  CAP$rec <- list()
  CAP$n_calls <- 0L
  CAP$t_cv <- 0
  CAP$t_post <- 0

  patched <- make_patched_cv(OVERRIDE, CAP)
  install_patch(patched)
  arm_label <- if (OVERRIDE) "targeted" else "baseline"
  t_run0 <- proc.time()[["elapsed"]]
  res <- tryCatch({
    set.seed(seed + 200000L)   # SAME pre-run seed for BOTH arms; only lambda differs.
    target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
    run_crossfit(
      split,
      n_folds = n_folds,
      communication_mode = mode,
      lambda_selection = "cv",
      lambda_rule = "min",
      verbose = FALSE,
      n_cores = 1L,                 # FORCED sequential so the capture/override is not lost to fork
      nlambda_init = nlambda_init,
      family = "binomial",
      A_val = 1L,
      use_lambda_cache = TRUE,
      precomputed_folds = folds,
      target_only_ps_cache = target_only_ps_cache,
      nuisance_lambda_rule = "min"
    )
  }, finally = restore_patch())

  ratios <- vapply(CAP$rec, function(r) {
    if (isTRUE(r$overridden) && is.finite(r$lambda_balancing) && r$lambda_balancing > 0) {
      r$lambda_entropy / r$lambda_balancing
    } else NA_real_
  }, numeric(1L))
  cat(sprintf("[arm:%s] run_crossfit=%.1fs  calls=%d  in_CV=%.1fs  in_post(balancing)=%.1fs  estimate=%.5f  se=%.5f\n",
              arm_label, proc.time()[["elapsed"]] - t_run0, CAP$n_calls,
              CAP$t_cv, CAP$t_post,
              as.numeric(res$estimate), as.numeric(res$se)))

  list(
    estimate = as.numeric(res$estimate),
    se = as.numeric(res$se),
    n_calls = CAP$n_calls,
    ratios = ratios
  )
}

# ---------------------------------------------------------------------------
# Paired arms on identical data/folds/seed. Baseline first (capture-only, fast), then targeted
# (pays the balancing-CV cost). The patch is installed AND restored around EACH arm separately so
# baseline truly has no override.
# ---------------------------------------------------------------------------
base <- run_arm(OVERRIDE = FALSE)
tgt  <- run_arm(OVERRIDE = TRUE)

bias_base <- base$estimate - truth
bias_tgt  <- tgt$estimate - truth
covered_base <- as.integer(abs(bias_base) <= 1.96 * base$se)
covered_tgt  <- as.integer(abs(bias_tgt)  <= 1.96 * tgt$se)

fin_ratios <- tgt$ratios[is.finite(tgt$ratios)]
median_lambda_ratio_ent_over_bal <- if (length(fin_ratios) > 0L)
  stats::median(fin_ratios) else NA_real_

row <- data.frame(
  seed = seed,
  n_total = n_total,
  K = K_sites,
  p = p,
  truth = truth,
  est_base = base$estimate,
  se_base = base$se,
  bias_base = bias_base,
  covered_base = covered_base,
  est_tgt = tgt$estimate,
  se_tgt = tgt$se,
  bias_tgt = bias_tgt,
  covered_tgt = covered_tgt,
  median_lambda_ratio_ent_over_bal = median_lambda_ratio_ent_over_bal,
  n_calls = tgt$n_calls,
  stringsAsFactors = FALSE
)
write.csv(row, out_file, row.names = FALSE)

cat(sprintf("[c2_targeted_lambda_coverage] seed=%d  truth=%.5f  | BASE est=%.5f se=%.5f bias=%+.5f cov=%d  | TGT est=%.5f se=%.5f bias=%+.5f cov=%d  | median(lambda_ent/lambda_bal)=%.4f  n_calls=%d\n",
            seed, truth,
            base$estimate, base$se, bias_base, covered_base,
            tgt$estimate, tgt$se, bias_tgt, covered_tgt,
            median_lambda_ratio_ent_over_bal, tgt$n_calls))
cat(sprintf("[c2_targeted_lambda_coverage] wrote %s\n", out_file))
