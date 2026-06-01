#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/probe.R
# ----------------------------------------------------------------------------
# Fast diagnostic probe for the one_round / two_round cross-fit blow-up.
# Walks the pipeline function-by-function and records magnitudes at every
# step so we can pinpoint WHERE the nuisance values explode (rather than
# only seeing the bias_mean=±100 aggregate result of run_simulation_study).
#
# THE ALGORITHM IS NOT MODIFIED — this script only calls existing exported
# / internal helpers in sequence and records their outputs.
#
# Three phases, each independently runnable:
#   A  data_sanity     One generate_simulation_data() call + structural
#                      checks (per-site sample sizes, per-arm sample sizes,
#                      treatment / outcome balance). Instant.
#   B  per_function    For ONE source and ONE (k1, k2) pair, call each
#                      nuisance fit (initial outcome, initial DR,
#                      calibrated DR, calibrated outcome) IN ISOLATION
#                      and record output magnitudes. Loops over secondary
#                      folds k2 = secondary so we see fold-to-fold
#                      variability. ~30 seconds.
#   C  end_to_end      For each of N seeds, run full run_crossfit() in
#                      both communication modes and record per-seed
#                      estimate + per-source mu_ts so we can see which
#                      seeds blow up and which source caused it.
#                      ~2 min × n_seeds.
#
# Usage:
#   Rscript diagnosis/bias/probe.R <phase> [n_total] [K] [p] [n_seeds]
# Defaults: phase=ALL  n_total=1000  K=3  p=20  n_seeds=10
# phase ∈ {A, B, C, ALL}
#
# All output goes to diagnosis/bias/probe_out/<phase>_n<N>_K<K>_p<P>.*
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
phase   <- if (length(args) >= 1L) toupper(args[[1]]) else "ALL"
n_total <- if (length(args) >= 2L) as.integer(args[[2]]) else 1000L
K       <- if (length(args) >= 3L) as.integer(args[[3]]) else 3L
p       <- if (length(args) >= 4L) as.integer(args[[4]]) else 20L
n_seeds <- if (length(args) >= 5L) as.integer(args[[5]]) else 10L

stopifnot(phase %in% c("A", "B", "C", "ALL"))
stopifnot(is.finite(n_total), n_total > 0,
          is.finite(K), K >= 2,
          is.finite(p), p >= 1,
          is.finite(n_seeds), n_seeds >= 1)

cat(sprintf("[probe] phase=%s n=%d K=%d p=%d n_seeds=%d\n",
            phase, n_total, K, p, n_seeds))

# ---- Load package via devtools::load_all so internal helpers are visible -
suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else if (requireNamespace("FACEHD", quietly = TRUE)) {
    library(FACEHD)
  } else {
    stop("Neither devtools nor FACEHD available.")
  }
})

out_dir <- "diagnosis/bias/probe_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setting_id <- sprintf("n%d_K%d_p%d", n_total, K, p)

# ---- Compact formatter for magnitude inspection --------------------------
mag <- function(x, label = "", limit_print = 8) {
  if (is.null(x)) {
    cat(sprintf("  %-32s NULL\n", label))
    return(invisible())
  }
  x <- as.numeric(x)
  finite_x <- x[is.finite(x)]
  n_na <- sum(!is.finite(x))
  if (length(finite_x) == 0L) {
    cat(sprintf("  %-32s ALL NON-FINITE (n=%d)\n", label, length(x)))
    return(invisible())
  }
  abs_max <- max(abs(finite_x))
  alert <- if (abs_max > 50) "  <<<--- LARGE" else ""
  cat(sprintf("  %-32s len=%-4d  min=%+.3e  max=%+.3e  |max|=%+.3e  non_finite=%d%s\n",
              label, length(x),
              min(finite_x), max(finite_x), abs_max, n_na, alert))
}

# ============================================================================
# PHASE A — data sanity
# ============================================================================
run_phase_A <- function() {
  cat("\n=========================  PHASE A: DATA SANITY  =========================\n")
  set.seed(101)
  data <- generate_simulation_data(
    n_total = n_total, K = K, p = p, config = "C1",
    estimand_type = "superpopulation", site_allocation = "model",
    transform_type = "mild", outcome_type = "binary",
    heterogeneity_type = "none", shift_strength = 0.5,
    dgp_type = "facehd"
  )
  data_split <- split_data_by_site(data)

  cat(sprintf("  site members: %s\n", paste(names(data_split), collapse = ", ")))
  for (s in names(data_split)) {
    d <- data_split[[s]]
    n_treated <- sum(d$A == 1L)
    n_control <- d$n - n_treated
    cat(sprintf("  site=%-2s  n=%-5d  n_A=1: %-4d  n_A=0: %-4d  E[Y]=%.3f\n",
                s, d$n, n_treated, n_control, mean(d$Y)))
  }
  cat(sprintf("  mu1_true=%.4f  mu0_true=%.4f  estimand_type=%s\n",
              data$mu1_true, data$mu0_true, data$estimand_type))
  cat("\n  >>> Effective sample sizes per (k1,k2) training fold (K_f=10):\n")
  source_sites <- setdiff(names(data_split), "t")
  for (s in source_sites) {
    n_site <- data_split[[s]]$n
    n_arm <- sum(data_split[[s]]$A == 1L)
    n_per_fold <- n_site / 10
    n_train_fold <- n_site * 8 / 10  # K_f - 2
    n_train_arm <- n_arm * 8 / 10
    cat(sprintf("  site=%-2s : per-fold ~%.0f obs, per (k1,k2)-train ~%.0f obs, A=1 in train ~%.0f obs (with p=%d)\n",
                s, n_per_fold, n_train_fold, n_train_arm, p))
  }
  cat("\n")
  invisible(list(data = data, data_split = data_split))
}

# ============================================================================
# PHASE B — per-function probe on ONE source, all (k1=1, k2 ∈ 2..K_f) pairs
# ============================================================================
run_phase_B <- function(data_split = NULL) {
  cat("\n=====================  PHASE B: PER-FUNCTION PROBE  ======================\n")
  if (is.null(data_split)) {
    set.seed(101)
    data <- generate_simulation_data(
      n_total = n_total, K = K, p = p, config = "C1",
      estimand_type = "superpopulation", site_allocation = "model",
      transform_type = "mild", outcome_type = "binary",
      heterogeneity_type = "none", shift_strength = 0.5,
      dgp_type = "facehd"
    )
    data_split <- split_data_by_site(data)
  }

  n_folds <- 10L
  source_sites <- setdiff(names(data_split), "t")
  target_data <- data_split[["t"]]

  fold_seed_base <- target_data$n + sum(sapply(source_sites, function(s) data_split[[s]]$n))
  target_folds <- partition_into_folds(target_data, n_folds, seed = fold_seed_base)
  source_folds <- setNames(
    lapply(seq_along(source_sites), function(i) {
      partition_into_folds(data_split[[source_sites[i]]], n_folds,
                           seed = fold_seed_base + i)
    }),
    source_sites
  )

  # Pick the first source, k1=1, sweep k2 over secondary folds
  s <- source_sites[[1]]
  k1 <- 1L
  secondary <- setdiff(1:n_folds, k1)
  glm_spec <- resolve_glm_family("binomial")
  family_int <- glm_spec$family_int
  link_int   <- glm_spec$link_int
  A_val <- 1L
  M_tau <- 10

  probe_records <- list()

  for (k2 in secondary) {
    cat(sprintf("\n-------- source=%s  k1=%d  k2=%d  --------\n", s, k1, k2))
    training_folds <- setdiff(1:n_folds, c(k1, k2))
    source_train <- combine_folds(source_folds[[s]], training_folds)
    source_calib <- materialize_fold(source_folds[[s]], k2)
    target_calib <- materialize_fold(target_folds, k2)

    cat(sprintf("  sizes: train=%d (A=1: %d)  calib=%d (A=1: %d)  target_calib=%d\n",
                source_train$n, sum(source_train$A == 1L),
                source_calib$n, sum(source_calib$A == 1L),
                target_calib$n))

    # ---- Step 1: fit_initial_outcome -----------------------------------
    alpha_init <- tryCatch(
      fit_initial_outcome(source_train$W_outcome, source_train$Y,
                          source_train$A, A_val = A_val, family = "binomial"),
      error = function(e) { cat("  [ERR] fit_initial_outcome: ", conditionMessage(e), "\n"); NULL })
    mag(alpha_init, "alpha_init (incl intercept)")

    # Prediction range on source_train (treated arm)
    if (!is.null(alpha_init)) {
      arm <- which(source_train$A == A_val)
      Xa <- cbind(1, source_train$W_outcome[arm, , drop = FALSE])
      eta_arm <- as.numeric(Xa %*% alpha_init)
      mu_arm  <- 1 / (1 + exp(-eta_arm))
      mag(eta_arm, "eta_arm (treated)")
      mag(mu_arm,  "logistic(eta_arm)")
    }

    # ---- Step 2: target summaries (mean_phi, mean_grad_psi_init) -------
    mean_phi <- c(1, colMeans(target_calib$Z_site))
    mag(mean_phi, "mean_phi (target k2)")

    mean_grad_psi <- tryCatch(
      .mean_glm_gradient_site_basis(target_calib$W_outcome, target_calib$Z_site,
                                    alpha_init, family_int, link_int),
      error = function(e) { cat("  [ERR] mean_grad_psi: ", conditionMessage(e), "\n"); NULL })
    mag(mean_grad_psi, "mean_grad_psi_init")

    # ---- Step 3: fit_initial_density_ratio -----------------------------
    gamma_init <- tryCatch(
      fit_initial_density_ratio(source_train$Z_site, source_train$A,
                                mean_phi, A_val = A_val),
      error = function(e) { cat("  [ERR] fit_initial_density_ratio: ", conditionMessage(e), "\n"); NULL })
    mag(gamma_init, "gamma_init")

    # Weight range on source_calib for the calibrated DR objective
    if (!is.null(gamma_init)) {
      Zc <- cbind(1, source_calib$Z_site)
      w_lin <- as.numeric(Zc %*% gamma_init)
      mag(w_lin, "Z_calib %*% gamma_init")
      mag(exp(w_lin), "exp(Z_calib %*% gamma_init)")
    }

    # ---- Step 4: fit_unified_density_ratio (calibrated) ----------------
    gamma_cal <- tryCatch(
      fit_unified_density_ratio(
        Z_site = source_calib$Z_site, A = source_calib$A,
        mean_grad_psi = mean_grad_psi, alpha_init = alpha_init,
        calibrated = TRUE, M_tau = M_tau,
        W_outcome = source_calib$W_outcome, A_val = A_val,
        family_int = family_int, link_int = link_int),
      error = function(e) { cat("  [ERR] fit_unified_density_ratio: ", conditionMessage(e), "\n"); NULL })
    mag(gamma_cal, "gamma_cal (calibrated)")

    # ---- Step 5: fit_unified_outcome (calibrated) ----------------------
    alpha_cal <- tryCatch(
      fit_unified_outcome(
        W_outcome = source_calib$W_outcome, Y = source_calib$Y,
        A = source_calib$A, A_val = A_val, gamma_s = gamma_init,
        family_int = family_int, link_int = link_int,
        calibrated = TRUE, M_tau = M_tau, Z_site = source_calib$Z_site),
      error = function(e) { cat("  [ERR] fit_unified_outcome: ", conditionMessage(e), "\n"); NULL })
    mag(alpha_cal, "alpha_cal (calibrated)")

    # ---- Step 6: prediction on TARGET main fold (k1) -------------------
    target_main <- materialize_fold(target_folds, k1)
    if (!is.null(alpha_cal)) {
      mu_pred_ts <- mean(predict_glm_cpp(target_main$W_outcome, alpha_cal,
                                         family_int, link_int))
      cat(sprintf("  mu_pred_ts (target main, k1=%d) = %+.4e\n", k1, mu_pred_ts))
    }

    probe_records[[paste0("k2_", k2)]] <- list(
      alpha_init = alpha_init, gamma_init = gamma_init,
      gamma_cal = gamma_cal, alpha_cal = alpha_cal
    )
  }

  rec_path <- file.path(out_dir, sprintf("phaseB_%s.RData", setting_id))
  save(probe_records, file = rec_path)
  cat(sprintf("\n[probe] phase B records saved: %s\n", rec_path))
  invisible(probe_records)
}

# ============================================================================
# PHASE C — end-to-end seed sweep: per-seed estimate + per-source mu_ts
# ============================================================================
run_phase_C <- function() {
  cat("\n=====================  PHASE C: END-TO-END SEED SWEEP  ===================\n")
  cat(sprintf("  Running %d seeds at (n=%d, K=%d, p=%d, C1, ss=0.5)\n",
              n_seeds, n_total, K, p))
  cat(sprintf("  Looking for blow-ups in one_round / two_round estimates.\n\n"))

  records <- vector("list", n_seeds)
  for (i in seq_len(n_seeds)) {
    seed <- 1000L + i
    set.seed(seed)
    t0 <- Sys.time()
    data <- generate_simulation_data(
      n_total = n_total, K = K, p = p, config = "C1",
      estimand_type = "superpopulation", site_allocation = "model",
      transform_type = "mild", outcome_type = "binary",
      heterogeneity_type = "none", shift_strength = 0.5,
      dgp_type = "facehd"
    )
    data_split <- split_data_by_site(data)

    res_or <- tryCatch(
      run_crossfit(data_split, n_folds = 10L,
                   communication_mode = "one_round",
                   verbose = FALSE),
      error = function(e) list(estimate = NA, error = conditionMessage(e)))
    res_tr <- tryCatch(
      run_crossfit(data_split, n_folds = 10L,
                   communication_mode = "two_round",
                   verbose = FALSE),
      error = function(e) list(estimate = NA, error = conditionMessage(e)))

    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    mu1 <- data$mu1_true

    # Per-source mu_ts (if components are present)
    source_mu_ts_or <- res_or$components$source_estimates %||% rep(NA, K)
    source_mu_ts_tr <- res_tr$components$source_estimates %||% rep(NA, K)

    or_alert <- if (is.finite(res_or$estimate) && abs(res_or$estimate - mu1) > 5) " *** BLOW UP" else ""
    tr_alert <- if (is.finite(res_tr$estimate) && abs(res_tr$estimate - mu1) > 5) " *** BLOW UP" else ""

    cat(sprintf("seed=%4d  elapsed=%5.1fs  mu1_true=%.4f\n", seed, elapsed, mu1))
    cat(sprintf("  one_round: est=%+10.4e  se=%+9.2e  source_mu_ts=[%s]%s\n",
                res_or$estimate %||% NA, res_or$se %||% NA,
                paste(sprintf("%+.2e", source_mu_ts_or), collapse = ", "),
                or_alert))
    cat(sprintf("  two_round: est=%+10.4e  se=%+9.2e  source_mu_ts=[%s]%s\n",
                res_tr$estimate %||% NA, res_tr$se %||% NA,
                paste(sprintf("%+.2e", source_mu_ts_tr), collapse = ", "),
                tr_alert))

    records[[i]] <- list(seed = seed, mu1_true = mu1,
                         one_round = res_or, two_round = res_tr,
                         elapsed_sec = elapsed)
  }

  rec_path <- file.path(out_dir, sprintf("phaseC_%s.RData", setting_id))
  save(records, file = rec_path)
  cat(sprintf("\n[probe] phase C records saved: %s\n", rec_path))
  invisible(records)
}

# ============================================================================
# Dispatch
# ============================================================================
t0_all <- Sys.time()
phase_a_out <- NULL
if (phase %in% c("A", "ALL")) phase_a_out <- run_phase_A()
if (phase %in% c("B", "ALL")) run_phase_B(phase_a_out$data_split)
if (phase %in% c("C", "ALL")) run_phase_C()
cat(sprintf("\n[probe] total elapsed: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0_all, units = "mins"))))
