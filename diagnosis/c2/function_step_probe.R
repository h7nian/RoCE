#!/usr/bin/env Rscript
# Function-by-function C2 diagnostic probe.
#
# This is intentionally diagnosis-only. It does not modify estimator code.  For
# one C2 setting and one seed, it records where centering fails: data/truth,
# target nuisance, source nuisance/transport, source estimates, and aggregation.

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
n_cores <- as.integer(args[[8L]])
modes <- strsplit(args[[9L]], ",", fixed = TRUE)[[1L]]
modes <- modes[nzchar(modes)]

stopifnot(n_total > 0L, K_sites > 0L, p > 0L, seed > 0L,
          n_folds >= 3L, nlambda_init > 0L, n_cores >= 1L)
if (!all(modes %in% c("one_round", "two_round"))) {
  stop("modes must be comma-separated one_round/two_round values.", call. = FALSE)
}

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) {
    library(RoCE)
  } else if (requireNamespace("devtools", quietly = FALSE)) {
    devtools::load_all(".")
  } else {
    stop("Neither installed RoCE nor devtools available.", call. = FALSE)
  }
})

roce_constant <- function(name, default) {
  ns <- asNamespace("RoCE")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  default
}

roce_function <- function(name) {
  ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("RoCE namespace does not contain required function '%s'.", name),
         call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}

LOGISTIC_CLIP <- roce_constant("LOGISTIC_CLIP", 50)
DIVISION_FLOOR <- roce_constant("DIVISION_FLOOR", 1e-10)
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("diagnosis", "c2", "c2_true_gamma_utils.R"))

for (fn_name in c("generate_simulation_data", "split_data_by_site",
                  "estimate_target_only_crossfit", "build_crossfit_folds",
                  "run_crossfit", "materialize_fold")) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, roce_function(fn_name))
  }
}

out_dir <- file.path("diagnosis", "c2", "function_probe", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d_kf%d", tag, n_total, K_sites, p, seed, n_folds
))

inv_logit <- function(eta) {
  eta <- pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP)
  1 / (1 + exp(-eta))
}

predict_alpha <- function(W, alpha) {
  inv_logit(as.numeric(cbind(1, as.matrix(W)) %*% as.numeric(alpha)))
}

rmse <- function(x) sqrt(mean(as.numeric(x)^2))

safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  as.numeric(stats::cor(x, y))
}

ess <- function(w) {
  w <- as.numeric(w)
  if (length(w) == 0L || sum(w^2) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

source_dr_components <- function(target_W, source_W, target_Z, source_Z,
                                 source_A, source_Y, target_truth_fold,
                                 gamma, alpha, label) {
  Wt <- as.matrix(target_W)
  Ws <- as.matrix(source_W)
  Zt <- as.matrix(target_Z)
  Zs <- as.matrix(source_Z)
  if (length(alpha) != ncol(Wt) + 1L || length(alpha) != ncol(Ws) + 1L) {
    stop(sprintf(
      "%s: alpha length %d incompatible with target/source W columns %d/%d.",
      label, length(alpha), ncol(Wt), ncol(Ws)
    ), call. = FALSE)
  }
  if (length(gamma) != ncol(Zs) + 1L) {
    stop(sprintf(
      "%s: gamma length %d incompatible with source Z columns %d.",
      label, length(gamma), ncol(Zs)
    ), call. = FALSE)
  }
  w <- exp(-as.numeric(cbind(1, Zs) %*% as.numeric(gamma)))
  treated <- as.numeric(source_A == 1L)
  mt <- predict_alpha(Wt, alpha)
  ms <- predict_alpha(Ws, alpha)
  raw_correction <- mean(treated * w * (source_Y - ms))
  denom <- mean(treated * w)
  normalized_correction <- raw_correction / max(denom, DIVISION_FLOOR)
  raw_est <- mean(mt) + raw_correction
  normalized_est <- mean(mt) + normalized_correction

  raw_weighted_z <- colMeans(Zs * as.numeric(treated * w))
  normalized_weighted_z <- colSums(Zs * as.numeric(treated * w)) /
    max(sum(treated * w), DIVISION_FLOOR)
  target_z <- colMeans(Zt)

  data.frame(
    variant = label,
    mu_pred = mean(mt),
    raw_correction = raw_correction,
    raw_estimate = raw_est,
    raw_bias_vs_fold_truth = raw_est - target_truth_fold,
    normalized_correction = normalized_correction,
    normalized_estimate = normalized_est,
    normalized_bias_vs_fold_truth = normalized_est - target_truth_fold,
    moment_intercept = denom,
    raw_z_gap = rmse(raw_weighted_z - target_z),
    normalized_z_gap = rmse(normalized_weighted_z - target_z),
    weight_mean = mean(w),
    weight_treated_mean = if (any(treated == 1)) mean(w[treated == 1]) else NA_real_,
    weight_treated_ess = ess(w[treated == 1]),
    weight_max = max(w),
    stringsAsFactors = FALSE
  )
}

write_csv <- function(df, suffix) {
  path <- paste0(prefix, "_", suffix, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("[function-probe] wrote %s (%d rows)\n", path, nrow(df)))
}

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
  dgp_type = "roce",
  warn_ignored = FALSE
)
split <- split_data_by_site(data)
source_sites <- setdiff(names(split), "t")
truth <- as.numeric(data$mu1_true)
target <- split[["t"]]

target_true_m <- predict_alpha(target$W_outcome_true, data$alpha1_true)
target_idx_global <- which(data$R == "t")
target_true_ps <- data$p_treat_true[target_idx_global]

data_rows <- data.frame(
  tag = tag,
  seed = seed,
  n_total = n_total,
  K = K_sites,
  p = p,
  n_folds = n_folds,
  truth_superpopulation = truth,
  truth_target_realized = as.numeric(data$mu1_realized),
  truth_gap_realized_minus_super = as.numeric(data$mu1_realized) - truth,
  n_target = target$n,
  n_source_min = min(vapply(source_sites, function(s) split[[s]]$n, integer(1L))),
  n_source_max = max(vapply(source_sites, function(s) split[[s]]$n, integer(1L))),
  z_cols = ncol(target$Z_site),
  w_cols = ncol(target$W_outcome),
  z_true_cols = ncol(target$Z_site_true),
  w_true_cols = ncol(target$W_outcome_true),
  target_treated_fraction = mean(target$A == 1L),
  target_true_ps_mean = mean(target_true_ps),
  stringsAsFactors = FALSE
)
write_csv(data_rows, "data_truth")

# Target-only decomposition: estimated m / estimated pi versus oracle pieces.
set.seed(seed + 100000L)
target_cf <- estimate_target_only_crossfit(target, n_folds = n_folds,
                                           family = "binomial", A_val = 1L)
target_est_m <- target_cf$m_pred
target_est_ps <- target_cf$prop_scores
target_y <- as.numeric(target$Y)
target_a <- as.numeric(target$A)
target_phi <- function(m, ps) {
  ps <- pmax(pmin(as.numeric(ps), 1 - DIVISION_FLOOR), DIVISION_FLOOR)
  as.numeric(m) + as.numeric(target_a == 1L) / ps * (target_y - as.numeric(m))
}
target_rows <- data.frame(
  tag = tag,
  seed = seed,
  estimator = c(
    "true_m_true_ps",
    "estimated_m_true_ps",
    "true_m_estimated_ps",
    "estimated_m_estimated_ps"
  ),
  estimate = c(
    mean(target_phi(target_true_m, target_true_ps)),
    mean(target_phi(target_est_m, target_true_ps)),
    mean(target_phi(target_true_m, target_est_ps)),
    mean(target_phi(target_est_m, target_est_ps))
  ),
  truth_superpopulation = truth,
  truth_target_realized = as.numeric(data$mu1_realized),
  m_mean_bias = c(
    0,
    mean(target_est_m - target_true_m),
    0,
    mean(target_est_m - target_true_m)
  ),
  m_rmse = c(0, rmse(target_est_m - target_true_m), 0,
             rmse(target_est_m - target_true_m)),
  ps_mean_bias = c(0, 0, mean(target_est_ps - target_true_ps),
                   mean(target_est_ps - target_true_ps)),
  ps_rmse = c(0, 0, rmse(target_est_ps - target_true_ps),
              rmse(target_est_ps - target_true_ps)),
  ps_cor = c(NA_real_, NA_real_, safe_cor(target_est_ps, target_true_ps),
             safe_cor(target_est_ps, target_true_ps)),
  stringsAsFactors = FALSE
)
target_rows$bias_vs_superpopulation <- target_rows$estimate - truth
target_rows$bias_vs_realized_target <- target_rows$estimate - as.numeric(data$mu1_realized)
write_csv(target_rows, "target_decomposition")

folds <- build_crossfit_folds(split, n_folds)
run_rows <- list()
weight_rows <- list()
source_rows <- list()
dr_rows <- list()

for (mode in modes) {
  set.seed(seed + if (identical(mode, "one_round")) 200000L else 300000L)
  res <- run_crossfit(
    split,
    n_folds = n_folds,
    communication_mode = mode,
    lambda_selection = "cv",
    lambda_rule = "min",
    verbose = FALSE,
    n_cores = n_cores,
    nlambda_init = nlambda_init,
    family = "binomial",
    A_val = 1L,
    use_lambda_cache = TRUE,
    precomputed_folds = folds,
    target_only_ps_cache = new.env(hash = TRUE, parent = emptyenv()),
    nuisance_lambda_rule = "min"
  )

  run_rows[[length(run_rows) + 1L]] <- data.frame(
    tag = tag,
    seed = seed,
    mode = mode,
    estimate = res$estimate,
    bias = res$estimate - truth,
    se = res$se,
    ci_lower = res$ci_lower,
    ci_upper = res$ci_upper,
    covered = truth >= res$ci_lower && truth <= res$ci_upper,
    target_estimate = res$target_only$estimate,
    target_bias = res$target_only$estimate - truth,
    source_estimate_mean = mean(as.numeric(res$source_estimates)),
    source_estimate_min = min(as.numeric(res$source_estimates)),
    source_estimate_max = max(as.numeric(res$source_estimates)),
    weight_sum = sum(as.numeric(res$weights)),
    weight_l1 = sum(abs(as.numeric(res$weights))),
    weight_min = min(as.numeric(res$weights)),
    weight_max = max(as.numeric(res$weights)),
    aggregation_lambda_mean = mean(as.numeric(res$fold_lambdas)),
    all_phi_mean = mean(as.numeric(res$all_phi_agg)),
    all_phi_site_centered_se = res$se,
    stringsAsFactors = FALSE
  )

  fw <- as.matrix(res$fold_weights)
  colnames(fw) <- source_sites
  for (k1 in seq_len(n_folds)) {
    for (source_idx in seq_along(source_sites)) {
      s <- source_sites[source_idx]
      weight_rows[[length(weight_rows) + 1L]] <- data.frame(
        tag = tag,
        seed = seed,
        mode = mode,
        fold = k1,
        source = s,
        weight = fw[k1, s],
        fold_lambda = res$fold_lambdas[k1],
        fold_target_estimate = res$intermediates$phase1$fold_target_estimate[k1],
        fold_source_estimate =
          res$intermediates$phase1$fold_source_estimates[k1, source_idx],
        fold_mu_pred = res$intermediates$phase1$mu_pred_ts[k1, source_idx],
        fold_delta = res$intermediates$phase1$delta_ts[k1, source_idx],
        stringsAsFactors = FALSE
      )
    }
  }

  for (k1 in seq_len(n_folds)) {
    target_fold <- materialize_fold(folds$target_folds, k1)
    target_true_fold <- predict_alpha(
      target$W_outcome_true[target_fold$original_idx, , drop = FALSE],
      data$alpha1_true
    )
    fold_truth <- mean(target_true_fold)

    for (s in source_sites) {
      src <- res$fold_results[[k1]]$source_results[[s]]
      source_fold <- materialize_fold(folds$source_folds[[s]], k1)
      target_true_W <- target$W_outcome_true[target_fold$original_idx, , drop = FALSE]
      target_true_Z <- target$Z_site_true[target_fold$original_idx, , drop = FALSE]
      source_true_W <- split[[s]]$W_outcome_true[source_fold$original_idx, , drop = FALSE]
      source_true_Z <- split[[s]]$Z_site_true[source_fold$original_idx, , drop = FALSE]
      source_true_fold <- predict_alpha(
        source_true_W,
        data$alpha1_true
      )
      target_est_m_fold <- predict_alpha(target_fold$W_outcome, src$alpha_ts)
      source_est_m_fold <- predict_alpha(source_fold$W_outcome, src$alpha_ts)

      source_rows[[length(source_rows) + 1L]] <- data.frame(
        tag = tag,
        seed = seed,
        mode = mode,
        fold = k1,
        source = s,
        fold_truth_realized = fold_truth,
        target_fold_estimate = res$fold_results[[k1]]$target_only$estimate,
        target_fold_bias = res$fold_results[[k1]]$target_only$estimate - fold_truth,
        source_estimate = src$mu_ts,
        source_bias = src$mu_ts - fold_truth,
        mu_pred = src$mu_pred_ts,
        mu_pred_bias = src$mu_pred_ts - fold_truth,
        delta = src$delta_ts,
        target_m_est_minus_true_mean = mean(target_est_m_fold - target_true_fold),
        target_m_est_rmse = rmse(target_est_m_fold - target_true_fold),
        source_m_est_minus_true_mean = mean(source_est_m_fold - source_true_fold),
        source_m_est_rmse = rmse(source_est_m_fold - source_true_fold),
        gamma_lambda = as.numeric(attr(src$gamma_s, "lambda_used") %||% NA_real_),
        alpha_lambda = as.numeric(attr(src$alpha_ts, "lambda_used") %||% NA_real_),
        gamma_support = sum(abs(as.numeric(src$gamma_s)[-1L]) > 1e-8),
        alpha_support = sum(abs(as.numeric(src$alpha_ts)[-1L]) > 1e-8),
        clip_any = isTRUE(src$correction_clip_diagnostics$any_clipped),
        stringsAsFactors = FALSE
      )

      true_gamma <- c2_true_source_calibration_gamma(data, s, A_val = 1L)
      variants <- rbind(
        source_dr_components(
          target_fold$W_outcome, source_fold$W_outcome,
          target_fold$Z_site, source_fold$Z_site,
          source_fold$A, source_fold$Y, fold_truth,
          src$gamma_s, src$alpha_ts, "estimated_gamma_estimated_alpha"
        ),
        source_dr_components(
          target_true_W, source_true_W,
          target_fold$Z_site, source_fold$Z_site,
          source_fold$A, source_fold$Y, fold_truth,
          src$gamma_s, data$alpha1_true, "estimated_gamma_true_alpha"
        ),
        source_dr_components(
          target_fold$W_outcome, source_fold$W_outcome,
          target_true_Z, source_true_Z,
          source_fold$A, source_fold$Y, fold_truth,
          true_gamma, src$alpha_ts, "true_gamma_estimated_alpha"
        ),
        source_dr_components(
          target_true_W, source_true_W,
          target_true_Z, source_true_Z,
          source_fold$A, source_fold$Y, fold_truth,
          true_gamma, data$alpha1_true, "true_gamma_true_alpha"
        )
      )
      variants$tag <- tag
      variants$seed <- seed
      variants$mode <- mode
      variants$fold <- k1
      variants$source <- s
      variants$fold_truth_realized <- fold_truth
      dr_rows[[length(dr_rows) + 1L]] <- variants
    }
  }
}

write_csv(do.call(rbind, run_rows), "crossfit_runs")
write_csv(do.call(rbind, weight_rows), "aggregation_weights")
write_csv(do.call(rbind, source_rows), "source_fold_diagnostics")
write_csv(do.call(rbind, dr_rows), "source_dr_oracle_variants")

cat("[function-probe] DONE\n")
