# test-c2-alignment-violation.R
#
# DIAGNOSTIC TEST (does NOT modify the main algorithm).
#
# Background
# ----------
# Production simulations show two_round_crossfit / one_round_crossfit under-cover
# under config C2 (outcome misspecified, weighting correct) at SOME (n, K, p, K_f)
# combinations -- e.g. n=10000, K=3, p=10, K_f=10 gives coverage 0.826 with
# bias = -0.0074, while pooled_dr / tilted_aipw at the same setting cover at
# 0.946 / 0.956 with near-zero bias.  Other (n, p, K_f) settings have normal C2
# coverage (verify_C2_n5000_K3_p10_kf3 ~ 0.95).  So C2 is not *intrinsically*
# broken; the question is what combination of regime + nuisance estimation
# triggers the bias.
#
# Three-layer diagnostic (per Sinian's framing)
# ---------------------------------------------
# Layer 1 (Proof scope check)
#   Compare alpha_init vs alpha_cal under C1/C2.  Only answers: does the
#   supplemental alignment assumption (supplemental.tex) hold in
#   population under each config?  alignment failure  =>  the current proof
#   does not cover this case.  It does NOT prove the estimator is biased,
#   because the DR identity
#     E_t[m(X)] + E_s[w*(X) (Y - m(X)) | A=A_val] = mu^1
#   holds for *any* working outcome m as long as w* is the correct density
#   ratio.
#
# Layer 2 (DR identity check with TRUE gamma)
#   For each (alpha_method in {alpha_init, alpha_cal, alpha_true}), compute
#   the DR bias using the conditional-source normalized TRUE density-ratio
#   coefficients (from DGP):
#     B_true_gamma(alpha_method) =
#       mean_t(psi(W; alpha_method))
#       + (1/n_s) sum_{i in source} I(A_i=A_val) w(Z_i; gamma_true) (Y_i - psi(W_i; alpha_method))
#       - mu^1_true.
#   If these are ~0 across alpha_method choices, the DR identity is intact
#   in the implemented estimand/site-weighting/treatment-arm setup.  Then
#   outcome misspecification alone is NOT the bias source.
#
# Layer 3 (DR identity check with CALIBRATED gamma)
#   Same as Layer 2 but with the *estimated* calibrated gamma (single-shot
#   version, not the production fold-summed block-design).  If
#   B_true_gamma ~ 0 and B_cal_gamma is materially non-zero, the bias source
#   is gamma estimation -- penalty shrinkage, calibration moment alignment,
#   high-dim regularization bias, truncation, etc.
#
# Dimension axis (per Sinian's p=10 observation)
# ----------------------------------------------
# Vary p across {10, 50}.  Production data shows C2 bias scales with p; the
# 3-layer diagnostic at each p tells us whether the high-p degradation is
# alignment, DR identity, or gamma estimation.
#
# What this file does NOT do
# --------------------------
# * Does NOT modify any algorithm code.
# * Does NOT use cross-fitting; everything is single-shot on full site data.
#   The point of this file is population-level identification of the bias
#   mechanism, not finite-sample fold accounting.
# * Does NOT claim the proof is wrong, only checks whether its sufficient
#   conditions are met under the configs being probed.
#
# Outputs
# -------
# tests/testthat/c2_alignment_output/<timestamp>_alignment_probe.csv  (long form)
# Console summary in log/test_*.out.
#
# Submit via:  ./test.sh --filter c2-alignment-violation

library(testthat)

ALIGNMENT_OUTPUT_DIR <- file.path("c2_alignment_output")
dir.create(ALIGNMENT_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

.alignment_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.alignment_numeric_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed) || parsed <= 0) default else parsed
}

.alignment_int_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  parsed <- parsed[is.finite(parsed) & parsed > 0L]
  if (length(parsed) == 0L) default else parsed
}

.alignment_config_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  bad <- setdiff(parsed, c("C1", "C2", "C3", "C4"))
  if (length(bad) > 0L) {
    stop(sprintf("%s contains invalid config(s): %s",
                 name, paste(bad, collapse = ", ")), call. = FALSE)
  }
  if (length(parsed) == 0L) default else parsed
}

.alignment_source_env <- function(default) {
  value <- Sys.getenv("ROCE_C2_ALIGNMENT_SOURCES", unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  if (length(parsed) == 0L) default else parsed
}


# ---------------------------------------------------------------------------
# .psi_outcome
# ---------------------------------------------------------------------------
# Outcome-model prediction psi(W; alpha) under binary logistic family.
# Includes intercept handling: alpha[1] is intercept, alpha[-1] are W coefs.
# ---------------------------------------------------------------------------
.psi_outcome <- function(W, alpha) {
  W_int <- cbind(1, as.matrix(W))
  alpha <- as.numeric(alpha)
  if (length(alpha) != ncol(W_int)) {
    stop(sprintf(".psi_outcome: alpha length %d but design has %d columns (need %d).",
                 length(alpha), ncol(W_int), ncol(W_int)))
  }
  eta <- as.numeric(W_int %*% alpha)
  eta_c <- pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP)
  1 / (1 + exp(-eta_c))
}


# ---------------------------------------------------------------------------
# .dr_weights
# ---------------------------------------------------------------------------
# Density-ratio weights w(Z; gamma) = exp(-Z_int %*% gamma) for one source.
# Returned raw (no normalization) -- the DR formula uses raw weights.
# ---------------------------------------------------------------------------
.dr_weights <- function(Z, gamma) {
  Z_int <- cbind(1, as.matrix(Z))
  gamma <- as.numeric(gamma)
  if (length(gamma) != ncol(Z_int)) {
    stop(sprintf(".dr_weights: gamma length %d but Z design has %d columns (need %d).",
                 length(gamma), ncol(Z_int), ncol(Z_int)))
  }
  eta <- as.numeric(Z_int %*% gamma)
  exp(-eta)
}

# ---------------------------------------------------------------------------
# .dr_bias
# ---------------------------------------------------------------------------
# Source-assisted DR estimator with arbitrary (alpha_method, gamma_method),
# minus the population truth mu1_true.  This is what the user calls B in the
# 3-layer diagnostic.
# ---------------------------------------------------------------------------
.dr_bias <- function(W_target, W_source, Z_source,
                     Y_source, A_source, A_val,
                     alpha, gamma, mu1_true) {
  mu_t   <- mean(.psi_outcome(W_target, alpha))
  m_src  <- .psi_outcome(W_source, alpha)
  w_src  <- .dr_weights(Z_source, gamma)
  treated <- as.numeric(A_source == A_val)
  delta  <- mean(treated * w_src * (Y_source - m_src))
  mu_hat <- mu_t + delta
  list(
    mu_t   = mu_t,
    delta  = delta,
    mu_hat = mu_hat,
    bias   = mu_hat - mu1_true,
    max_w  = max(w_src),
    mean_w_treated = if (sum(treated) > 0) mean(w_src[treated == 1]) else NA_real_
  )
}


# ---------------------------------------------------------------------------
# .probe_one_realization
# ---------------------------------------------------------------------------
# For one DGP realization, fit alpha_init, gamma_init, alpha_cal, gamma_cal
# on the FULL source data (no folds), then compute the 3-layer diagnostic
# for the requested source sites.  Returns a long data.frame.
# ---------------------------------------------------------------------------
.probe_one_realization <- function(data, source_names, A_val = 1L) {
  data_split <- split_data_by_site(data)
  if (!"t" %in% names(data_split)) {
    stop(".probe_one_realization: target site 't' is missing from split_data_by_site(data).",
         call. = FALSE)
  }
  missing_sources <- setdiff(source_names, names(data_split))
  if (length(missing_sources) > 0L) {
    stop(sprintf(
      ".probe_one_realization: requested source site(s) not found: %s. Available sites: %s.",
      paste(missing_sources, collapse = ", "),
      paste(names(data_split), collapse = ", ")
    ), call. = FALSE)
  }

  target_data <- data_split[["t"]]
  # Fitted bases (what alpha_init/alpha_cal/gamma_init/gamma_cal were fit on)
  W_target      <- as.matrix(target_data$W_outcome)
  Z_target      <- as.matrix(target_data$Z_site)
  # TRUE DGP bases (what alpha_true/gamma_true are calibrated to).  Under C2/C4
  # these differ from the fitted bases; under C1/C3 they coincide.  Using the
  # fitted basis with true-coef rows inflates "oracle" bias artificially --
  # that was the failure mode of the previous run.
  W_target_true <- as.matrix(target_data$W_outcome_true %||% target_data$W_outcome)
  Z_target_true <- as.matrix(target_data$Z_site_true   %||% target_data$Z_site)

  mu1_true <- as.numeric(data$mu1_true)
  alpha_true_global <- as.numeric(data$alpha1_true)
  gamma_params <- data$gamma_params

  per_source_rows <- list()
  for (s in source_names) {
    source_data <- data_split[[s]]
    W_source      <- as.matrix(source_data$W_outcome)
    Z_source      <- as.matrix(source_data$Z_site)
    W_source_true <- as.matrix(source_data$W_outcome_true %||% source_data$W_outcome)
    Z_source_true <- as.matrix(source_data$Z_site_true   %||% source_data$Z_site)
    Y_source <- source_data$Y
    A_source <- source_data$A

    site_num <- as.integer(gsub("s", "", s))
    gamma_key <- paste0("s", site_num, "_", A_val)
    gamma_true_raw <- as.numeric(gamma_params[[gamma_key]])
    if (is.null(gamma_true_raw) || length(gamma_true_raw) == 0L) {
      stop(sprintf(".probe_one_realization: gamma_true missing for key '%s'.", gamma_key))
    }
    gamma_true <- c2_true_source_calibration_gamma(
      data, s, A_val = A_val, gamma_raw = gamma_true_raw
    )

    # ---- LAYER 1: alpha_init vs alpha_cal (population alignment probe) ----
    alpha_init <- as.numeric(fit_initial_outcome(
      W_outcome = W_source, Y = Y_source, A = A_source,
      A_val = A_val, family = "binomial"
    ))

    mean_phi_target <- c(1, colMeans(Z_target))
    probe_lambda <- 0.01
    gamma_init <- as.numeric(fit_initial_density_ratio(
      Z_site = Z_source, A = A_source, mean_phi = mean_phi_target,
      lambda = probe_lambda, A_val = A_val
    ))

    alpha_cal <- as.numeric(fit_unified_outcome(
      W_outcome = W_source, Y = Y_source, A = A_source, A_val = A_val,
      gamma_s = gamma_init,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT,
      Z_site = Z_source
    ))

    # mean_grad_psi computed from TARGET side with alpha_init plug-in, matching
    # the calibrated-gamma loss formulation.
    mean_grad_psi <- .mean_glm_gradient_site_basis(
      W_outcome = W_target, Z_site = Z_target, alpha = alpha_init,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    )
    gamma_cal_fit <- fit_unified_density_ratio(
      Z_site = Z_source, A = A_source, mean_grad_psi = mean_grad_psi,
      alpha_init = alpha_init,
      lambda = probe_lambda,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT,
      W_outcome = W_source, A_val = A_val,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    )
    gamma_cal <- as.numeric(gamma_cal_fit)
    lambda_gamma_cal <- as.numeric(attr(gamma_cal_fit, "lambda_used") %||% NA_real_)

    # ----- gamma-side CV / regularization diagnostics ---------------------
    # (Sinian's nuisance-CV concern: is the chosen lambda_gamma trading off
    # the right thing for the final DR estimator?)
    gamma_cal_support <- sum(abs(gamma_cal[-1L]) > 1e-8)  # exclude intercept
    w_cal   <- .dr_weights(Z_source, gamma_cal)
    w_true  <- .dr_weights(Z_source_true, gamma_true)
    treated_mask <- as.numeric(A_source == A_val)

    .ess <- function(w, mask) {
      ww <- w * mask
      s1 <- sum(ww)
      s2 <- sum(ww^2)
      if (s2 <= 0) NA_real_ else (s1^2) / s2
    }
    ess_cal_treated  <- .ess(w_cal,  treated_mask)
    ess_true_treated <- .ess(w_true, treated_mask)

    # Balance gap on Z_site: target feature mean vs source weighted mean
    # restricted to A=A_val.  If gamma_cal balances correctly, this is ~0.
    Z_target_mean <- colMeans(Z_target)
    .weighted_treated_mean <- function(Z, w, mask) {
      ww <- w * mask
      if (sum(ww) <= 0) return(rep(NA_real_, ncol(Z)))
      colSums(Z * ww) / sum(ww)
    }
    Z_src_w_cal  <- .weighted_treated_mean(Z_source, w_cal,  treated_mask)
    Z_src_w_true <- .weighted_treated_mean(Z_source, w_true, treated_mask)
    balance_gap_cal_l2  <- sqrt(mean((Z_target_mean - Z_src_w_cal)^2))
    balance_gap_true_l2 <- sqrt(mean((Z_target_mean - Z_src_w_true)^2))

    coef_gap_alpha <- sqrt(mean((alpha_init - alpha_cal)^2))
    pred_init <- mean(.psi_outcome(W_target, alpha_init))
    pred_cal  <- mean(.psi_outcome(W_target, alpha_cal))
    pred_gap_alpha <- abs(pred_init - pred_cal)
    pred_signed_alpha <- pred_init - pred_cal

    # ---- LAYERS 2 + 3: DR identity under (alpha_method, gamma_method) ----
    # All 4 informative combinations of {true, cal} x {true, cal}.
    # We also include (alpha_init, gamma_cal) as a sanity reference -- that
    # one is closest to what an uncalibrated outcome + calibrated gamma DR
    # would produce.
    # Each (alpha_method, gamma_method) requires its OWN basis pair:
    #  - alpha_true coefficients must apply to the true outcome basis (X_dagger).
    #  - alpha_init / alpha_cal coefficients must apply to the fitted basis.
    #  - gamma_true similarly to the true site basis; gamma_init / gamma_cal to
    #    the fitted site basis.  Mixing breaks the DR identity *by construction*.
    combos <- list(
      list(am = "true", gm = "true", a = alpha_true_global, g = gamma_true,
           W_t = W_target_true, W_s = W_source_true, Z_s = Z_source_true),
      list(am = "true", gm = "raw",  a = alpha_true_global, g = gamma_true_raw,
           W_t = W_target_true, W_s = W_source_true, Z_s = Z_source_true),
      list(am = "true", gm = "cal",  a = alpha_true_global, g = gamma_cal,
           W_t = W_target_true, W_s = W_source_true, Z_s = Z_source),
      list(am = "cal",  gm = "true", a = alpha_cal,         g = gamma_true,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source_true),
      list(am = "cal",  gm = "raw",  a = alpha_cal,         g = gamma_true_raw,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source_true),
      list(am = "cal",  gm = "cal",  a = alpha_cal,         g = gamma_cal,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source),
      list(am = "init", gm = "true", a = alpha_init,        g = gamma_true,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source_true),
      list(am = "init", gm = "raw",  a = alpha_init,        g = gamma_true_raw,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source_true),
      list(am = "init", gm = "cal",  a = alpha_init,        g = gamma_cal,
           W_t = W_target,      W_s = W_source,      Z_s = Z_source)
    )

    dr_rows <- list()
    for (cc in combos) {
      dr <- tryCatch(
        .dr_bias(cc$W_t, cc$W_s, cc$Z_s,
                 Y_source, A_source, A_val,
                 alpha = cc$a, gamma = cc$g, mu1_true = mu1_true),
        error = function(e) {
          stop(sprintf(
            ".probe_one_realization: .dr_bias failed for source=%s, alpha_method=%s, gamma_method=%s: %s",
            s, cc$am, cc$gm, conditionMessage(e)
          ), call. = FALSE)
        }
      )
      dr_rows[[length(dr_rows) + 1L]] <- data.frame(
        source = s,
        alpha_method = cc$am,
        gamma_method = cc$gm,
        mu_t        = dr$mu_t,
        delta       = dr$delta,
        mu_hat      = dr$mu_hat,
        bias        = dr$bias,
        max_w       = dr$max_w,
        mean_w_trtd = dr$mean_w_treated,
        stringsAsFactors = FALSE
      )
    }
    dr_df <- do.call(rbind, dr_rows)

    # gamma quality diagnostics
    gamma_cal_abs_max <- max(abs(gamma_cal))
    gamma_diff_max <- max(abs(gamma_cal - gamma_true))
    gamma_true_abs_max <- max(abs(gamma_true))
    gamma_true_raw_abs_max <- max(abs(gamma_true_raw))

    # Annotate the DR rows with the alpha-layer + gamma-quality summaries so
    # everything for this (rep, source) lives in one long-form row set.
    dr_df$coef_gap_alpha     <- coef_gap_alpha
    dr_df$pred_gap_alpha     <- pred_gap_alpha
    dr_df$pred_signed_alpha  <- pred_signed_alpha
    dr_df$mu1_true           <- mu1_true
    dr_df$gamma_cal_abs_max  <- gamma_cal_abs_max
    dr_df$gamma_true_abs_max <- gamma_true_abs_max
    dr_df$gamma_true_raw_abs_max <- gamma_true_raw_abs_max
    dr_df$gamma_diff_max     <- gamma_diff_max
    dr_df$lambda_gamma_cal   <- lambda_gamma_cal
    dr_df$gamma_cal_support  <- gamma_cal_support
    dr_df$ess_cal_treated    <- ess_cal_treated
    dr_df$ess_true_treated   <- ess_true_treated
    dr_df$balance_gap_cal_l2  <- balance_gap_cal_l2
    dr_df$balance_gap_true_l2 <- balance_gap_true_l2
    dr_df$n_source           <- nrow(W_source)
    dr_df$n_source_treated   <- sum(A_source == A_val)

    per_source_rows[[length(per_source_rows) + 1L]] <- dr_df
  }

  if (length(per_source_rows) == 0L) {
    stop(".probe_one_realization: no source-site diagnostic rows were produced.",
         call. = FALSE)
  }
  do.call(rbind, per_source_rows)
}


# ---------------------------------------------------------------------------
# .run_probe_grid
# ---------------------------------------------------------------------------
# Sweep across n_total, p, config, replications.  Records the cartesian
# product as a single long data.frame.
# ---------------------------------------------------------------------------
.run_probe_grid <- function(grid, n_reps, K_sites, source_subset,
                            seed_base = 0L) {
  rows <- list()
  for (cell_idx in seq_len(nrow(grid))) {
    n_total <- grid$n_total[cell_idx]
    p       <- grid$p[cell_idx]
    config  <- grid$config[cell_idx]

    cat(sprintf("[c2-alignment]   running n=%d, p=%d, config=%s (%d reps)\n",
                n_total, p, config, n_reps))

    cell_seed_base <- seed_base +
      1009L * which(unique(grid$n_total) == n_total) +
      2027L * which(unique(grid$p)       == p)       +
      3041L * match(config, unique(grid$config))

    for (rep_idx in seq_len(n_reps)) {
      seed <- cell_seed_base + rep_idx * 7L
      set.seed(seed)
      data <- tryCatch({
        generate_simulation_data(
          n_total = n_total, K = K_sites, p = p,
          config = config,
          estimand_type = "superpopulation",
          site_allocation = "model",
          transform_type  = "mild",
          outcome_type    = "binary",
          heterogeneity_type = "none",
          shift_strength  = 0.5,
          dgp_type        = "roce",
          warn_ignored    = FALSE
        )
      }, error = function(e) {
        stop(sprintf(
          ".run_probe_grid: generate_simulation_data failed for n_total=%d, p=%d, config=%s, rep=%d: %s",
          n_total, p, config, rep_idx, conditionMessage(e)
        ), call. = FALSE)
      })

      avail_sources <- setdiff(names(split_data_by_site(data)), "t")
      chosen <- intersect(source_subset, avail_sources)
      missing_sources <- setdiff(source_subset, avail_sources)
      if (length(missing_sources) > 0L) {
        stop(sprintf(
          ".run_probe_grid: requested source_subset site(s) not found for n_total=%d, p=%d, config=%s, rep=%d: %s. Available sources: %s.",
          n_total, p, config, rep_idx,
          paste(missing_sources, collapse = ", "),
          paste(avail_sources, collapse = ", ")
        ), call. = FALSE)
      }

      probe <- tryCatch(
        .probe_one_realization(data, chosen, A_val = 1L),
        error = function(e) {
          stop(sprintf(
            ".run_probe_grid: .probe_one_realization failed for n_total=%d, p=%d, config=%s, rep=%d: %s",
            n_total, p, config, rep_idx, conditionMessage(e)
          ), call. = FALSE)
        }
      )

      probe$n_total <- n_total
      probe$p       <- p
      probe$config  <- config
      probe$rep     <- rep_idx
      probe$seed    <- seed
      rows[[length(rows) + 1L]] <- probe
    }
  }
  if (length(rows) == 0L) {
    stop(".run_probe_grid: no diagnostic rows were produced.", call. = FALSE)
  }
  do.call(rbind, rows)
}


# ---------------------------------------------------------------------------
# .summarise
# ---------------------------------------------------------------------------
.summarise <- function(df) {
  if (is.null(df) || nrow(df) == 0L) {
    stop(".summarise: diagnostic data frame is empty.", call. = FALSE)
  }
  key <- with(df, paste(config, n_total, p, alpha_method, gamma_method, sep = "|"))
  by(df, key, function(sub) {
    data.frame(
      config       = sub$config[1L],
      n_total      = sub$n_total[1L],
      p            = sub$p[1L],
      alpha_method = sub$alpha_method[1L],
      gamma_method = sub$gamma_method[1L],
      bias_mean    = mean(sub$bias),
      bias_sd      = sd(sub$bias),
      bias_median  = median(sub$bias),
      mu_hat_mean  = mean(sub$mu_hat),
      n_rows       = nrow(sub),
      gamma_diff_max_mean = mean(sub$gamma_diff_max),
      gamma_cal_abs_max_mean = mean(sub$gamma_cal_abs_max),
      max_w_p95    = quantile(sub$max_w, 0.95, names = FALSE),
      stringsAsFactors = FALSE
    )
  }) -> by_list
  do.call(rbind, lapply(by_list, identity))
}


# ---------------------------------------------------------------------------
# MAIN PROBE
# ---------------------------------------------------------------------------
.c2_alignment_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_ALIGNMENT", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-alignment", filter, fixed = TRUE)
}

test_that("three-layer C1/C2 DR-identity probe across (n_total, p)", {
  skip_if_not(.c2_alignment_enabled(),
              message = "set ROCE_RUN_C2_ALIGNMENT=1 or run with --filter c2-alignment")

  # Grid: vary p (Sinian's observation: p=10 is fine, p=50 shows bias) and
  # n_total around the production focal points.
  K_SITES       <- .alignment_int_env("ROCE_C2_ALIGNMENT_K", 3L)
  N_REPS        <- .alignment_int_env("ROCE_C2_ALIGNMENT_REPS", 20L)
  SOURCE_SUBSET <- .alignment_source_env(c("s1", "s2"))  # 2 of K=3 sources by default

  grid <- expand.grid(
    n_total = .alignment_int_vector_env("ROCE_C2_ALIGNMENT_N", c(5000L, 10000L)),
    p       = .alignment_int_vector_env("ROCE_C2_ALIGNMENT_P", c(10L, 50L)),
    config  = .alignment_config_env("ROCE_C2_ALIGNMENT_CONFIGS", c("C1", "C2")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  cat(sprintf("\n[c2-alignment] Grid: %d cells, %d reps each, %d sources/cell, K=%d\n",
              nrow(grid), N_REPS, length(SOURCE_SUBSET), K_SITES))
  cat("[c2-alignment] Layer 1: alpha alignment (info, NOT a bias claim)\n")
  cat("[c2-alignment] Layer 2: B with TRUE conditional-source gamma (probes DR identity)\n")
  cat("[c2-alignment] Layer 3: B with CAL  gamma (probes gamma estimation)\n\n")

  t0 <- Sys.time()
  df <- .run_probe_grid(grid, N_REPS, K_SITES, SOURCE_SUBSET, seed_base = 10000L)
  elapsed <- difftime(Sys.time(), t0, units = "secs")
  cat(sprintf("[c2-alignment] Total wall time: %.1f sec\n", as.numeric(elapsed)))

  expect_true(!is.null(df) && nrow(df) > 0L,
              info = "Probe produced no rows -- see skip messages above.")

  # ---- Persist long-form CSV ----
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_path <- file.path(ALIGNMENT_OUTPUT_DIR,
                        paste0(ts, "_alignment_probe.csv"))
  write.csv(df, out_path, row.names = FALSE)
  cat(sprintf("\n[c2-alignment] Long-form CSV: %s  (%d rows)\n", out_path, nrow(df)))

  # ---- Summary table ----
  agg <- .summarise(df)
  agg_sorted <- agg[order(agg$config, agg$p, agg$n_total,
                          agg$alpha_method, agg$gamma_method), ]
  rownames(agg_sorted) <- NULL
  cat("\n[c2-alignment] Mean bias by (config, p, n_total, alpha_method, gamma_method):\n")
  print(agg_sorted, row.names = FALSE)

  agg_path <- file.path(ALIGNMENT_OUTPUT_DIR,
                        paste0(ts, "_alignment_probe_summary.csv"))
  write.csv(agg_sorted, agg_path, row.names = FALSE)
  cat(sprintf("\n[c2-alignment] Summary CSV: %s\n", agg_path))

  # ---- Headline diagnostic readout ----
  # Compare {true, true} bias to {cal, cal} bias, by (config, p, n_total).
  cat("\n[c2-alignment] HEADLINE -- B(alpha_method, gamma_method) at n=10000:\n")
  for (cfg in c("C1", "C2")) {
    for (pp in c(10L, 50L)) {
      sub <- agg_sorted[agg_sorted$config == cfg &
                          agg_sorted$p == pp &
                          agg_sorted$n_total == 10000L, , drop = FALSE]
      if (nrow(sub) == 0L) next
      cat(sprintf("  %s p=%d:\n", cfg, pp))
      for (i in seq_len(nrow(sub))) {
        cat(sprintf("    alpha=%-4s gamma=%-4s  bias_mean=%+.5f  bias_sd=%.5f\n",
                    sub$alpha_method[i], sub$gamma_method[i],
                    sub$bias_mean[i], sub$bias_sd[i]))
      }
    }
  }

  cat("\n[c2-alignment] INTERPRETATION GUIDE\n")
  cat("  * gamma=true is the conditional-source normalized DGP gamma; gamma=raw is printed only as a convention check.\n")
  cat("  * If B(true, true) ~ 0 across the board: DGP and DR formula are aligned (sanity).\n")
  cat("  * If B(cal, true) ~ B(true, true) ~ 0: alpha misspec under C2 alone is NOT the issue.\n")
  cat("  * If B(cal, true) ~ 0 but B(cal, cal) is materially non-zero in C2/highP:\n")
  cat("      => gamma estimation (penalty/calibration/regularization) is the bias source.\n")
  cat("  * If B(cal, true) is materially non-zero already in C2:\n")
  cat("      => DR identity does NOT hold under the implemented site weighting/estimand;\n")
  cat("         go back to the DR derivation under C2's misspec, NOT to gamma estimation.\n")
  cat("  * Compare p=10 vs p=50 columns to see if any bias scales with dimension.\n\n")

  # ---- Sanity-only assertions (loose; the file is diagnostic, not gating) ----
  # 1) B(true, true) should be small in magnitude in the default large-n,
  # multi-rep diagnostic.  Tiny one-rep smoke jobs can be noisy because this
  # subtracts the superpopulation truth, not the realized target-sample mean.
  bias_true_true <- agg_sorted$bias_mean[
    agg_sorted$alpha_method == "true" & agg_sorted$gamma_method == "true"
  ]
  expect_true(all(is.finite(bias_true_true)),
              info = "B(true, true) produced non-finite values.")
  strict_true_check <- N_REPS >= 5L && min(grid$n_total) >= 5000L
  if (strict_true_check) {
    true_tol <- .alignment_numeric_env("ROCE_C2_ALIGNMENT_TRUE_TOL", 0.04)
    expect_true(
      all(abs(bias_true_true) < true_tol),
      info = sprintf("B(true, true) should be near zero across large diagnostic cells; got max |bias|=%.5f with tol=%.5f",
                     max(abs(bias_true_true)), true_tol)
    )
  } else {
    cat(sprintf("[c2-alignment] Skipping strict B(true,true) gate for small smoke grid; max |bias|=%.5f\n",
                max(abs(bias_true_true))))
  }
  # 2) Probe must have produced at least one row per (config, p, n_total, alpha_method, gamma_method).
  expect_true(
    all(table(df$config, df$p, df$n_total, df$alpha_method, df$gamma_method) > 0L),
    info = "Some (config, p, n_total, alpha_method, gamma_method) cells produced zero probe rows."
  )
})
