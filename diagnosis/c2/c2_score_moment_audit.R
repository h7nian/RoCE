#!/usr/bin/env Rscript
# C2 score/moment audit for source-assisted RoCE components.
#
# This probe is diagnosis-only. It checks fold/source-level nuisance scores and
# DR correction components for the problematic C2 p=10 settings without
# modifying package source code.

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
use_lambda_cache <- tolower(Sys.getenv("C2_SCORE_AUDIT_USE_LAMBDA_CACHE", "true")) %in%
  c("1", "true", "yes")

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

needed_functions <- c(
  "generate_simulation_data", "split_data_by_site", "build_crossfit_folds",
  "run_crossfit", "materialize_fold", "predict_glm_cpp",
  ".mean_glm_gradient_site_basis"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, roce_function(fn_name))
  }
}

LOGISTIC_CLIP <- roce_constant("LOGISTIC_CLIP", 50)
M_TAU_DEFAULT <- roce_constant("M_TAU_DEFAULT", 10)
M_TAU_INFERENCE_DEFAULT <- roce_constant("M_TAU_INFERENCE_DEFAULT", Inf)
# Diagnostic: override the calibrated-loss truncation M_tau via env C2_M_TAU
# (default = package M_TAU_DEFAULT). Used to study truncation sensitivity vs SMMAL.
c2_m_tau_train <- suppressWarnings(as.numeric(Sys.getenv("C2_M_TAU", "")))
if (!is.finite(c2_m_tau_train) || c2_m_tau_train <= 0) c2_m_tau_train <- M_TAU_DEFAULT
FAMILY_BINOMIAL <- roce_constant("FAMILY_BINOMIAL", 1L)
LINK_LOGIT <- roce_constant("LINK_LOGIT", 1L)
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("diagnosis", "c2", "c2_true_gamma_utils.R"))

out_root <- Sys.getenv(
  "C2_SCORE_AUDIT_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "score_moment_audit")
)
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d_kf%d", tag, n_total, K_sites, p, seed, n_folds
))

clip_linear <- function(eta, bound) {
  eta <- as.numeric(eta)
  if (!is.finite(bound)) return(pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP))
  pmax(pmin(eta, bound), -bound)
}

inv_logit <- function(eta) {
  eta <- pmax(pmin(as.numeric(eta), LOGISTIC_CLIP), -LOGISTIC_CLIP)
  1 / (1 + exp(-eta))
}

logit_derivative <- function(eta) {
  mu <- inv_logit(eta)
  mu * (1 - mu)
}

predict_alpha <- function(W, alpha) {
  inv_logit(as.numeric(cbind(1, as.matrix(W)) %*% as.numeric(alpha)))
}

score_norm_row <- function(prefix, score) {
  score <- as.numeric(score)
  tail_score <- if (length(score) > 1L) score[-1L] else numeric(0)
  out <- c(
    intercept = score[[1L]],
    linf = max(abs(score)),
    l2 = sqrt(mean(score^2)),
    tail_linf = if (length(tail_score)) max(abs(tail_score)) else 0,
    tail_l2 = if (length(tail_score)) sqrt(mean(tail_score^2)) else 0
  )
  names(out) <- paste(prefix, names(out), sep = "_")
  as.list(out)
}

fold_gamma_score <- function(target_fold, source_fold, gamma_final, alpha_init) {
  target_mean <- .mean_glm_gradient_site_basis(
    target_fold$W_outcome, target_fold$Z_site,
    alpha_init, FAMILY_BINOMIAL, LINK_LOGIT
  )

  Z_source_int <- cbind(1, as.matrix(source_fold$Z_site))
  W_source_int <- cbind(1, as.matrix(source_fold$W_outcome))
  alpha_init <- as.numeric(alpha_init)
  gamma_final <- as.numeric(gamma_final)

  eta_alpha <- as.numeric(W_source_int %*% alpha_init)
  psi_prime <- logit_derivative(clip_linear(eta_alpha, M_TAU_DEFAULT))
  eta_gamma <- clip_linear(as.numeric(Z_source_int %*% gamma_final), LOGISTIC_CLIP)
  weights <- as.numeric(source_fold$A == 1L) * exp(-eta_gamma) * psi_prime
  source_mean <- colMeans(sweep(Z_source_int, 1L, weights, "*"))

  target_mean - source_mean
}

fold_alpha_score <- function(source_fold, alpha_final, gamma_for_weights,
                             gamma_truncation = M_TAU_DEFAULT) {
  W_int <- cbind(1, as.matrix(source_fold$W_outcome))
  Z_int <- cbind(1, as.matrix(source_fold$Z_site))
  alpha_final <- as.numeric(alpha_final)
  gamma_for_weights <- as.numeric(gamma_for_weights)

  mu_hat <- inv_logit(as.numeric(W_int %*% alpha_final))
  eta_gamma <- clip_linear(as.numeric(Z_int %*% gamma_for_weights), gamma_truncation)
  weights <- as.numeric(source_fold$A == 1L) * exp(-eta_gamma)
  residual <- mu_hat - as.numeric(source_fold$Y)

  colMeans(sweep(W_int, 1L, weights * residual, "*"))
}

mean_score <- function(scores) {
  Reduce(`+`, scores) / length(scores)
}

fold_source_pi <- function(data, split, source_site, source_fold) {
  source_all_idx <- which(data$R == source_site)
  source_pi_all <- c2_true_source_treatment_propensity(data, source_site, A_val = 1L)
  source_pi_all[source_all_idx[source_fold$original_idx]]
}

correction_moments <- function(target_fold, source_fold, data, split,
                               source_site, gamma_final, alpha_final) {
  target_true_W <- split$t$W_outcome_true[target_fold$original_idx, , drop = FALSE]
  source_true_W <- split[[source_site]]$W_outcome_true[source_fold$original_idx, , drop = FALSE]
  target_true_m <- predict_alpha(target_true_W, data$alpha1_true)
  source_true_m <- predict_alpha(source_true_W, data$alpha1_true)

  target_m_hat <- predict_alpha(target_fold$W_outcome, alpha_final)
  source_m_hat <- predict_alpha(source_fold$W_outcome, alpha_final)

  Z_int <- cbind(1, as.matrix(source_fold$Z_site))
  eta_gamma <- clip_linear(as.numeric(Z_int %*% as.numeric(gamma_final)),
                           M_TAU_INFERENCE_DEFAULT)
  weights <- exp(-eta_gamma)
  observed_weights <- as.numeric(source_fold$A == 1L) * weights
  true_pi_weights <- fold_source_pi(data, split, source_site, source_fold) * weights

  source_residual_true <- source_true_m - source_m_hat
  source_residual_observed <- as.numeric(source_fold$Y) - source_m_hat
  target_residual_true <- target_true_m - target_m_hat

  required_delta <- mean(target_residual_true)
  observed_delta <- mean(observed_weights * source_residual_observed)
  model_delta_true_pi <- mean(true_pi_weights * source_residual_true)
  model_delta_observed_A <- mean(observed_weights * source_residual_true)

  data.frame(
    target_true_mean = mean(target_true_m),
    target_hat_mean = mean(target_m_hat),
    source_true_mean = mean(source_true_m),
    source_hat_mean = mean(source_m_hat),
    required_delta = required_delta,
    observed_delta = observed_delta,
    model_delta_true_pi = model_delta_true_pi,
    model_delta_observed_A = model_delta_observed_A,
    source_bias_observed = observed_delta - required_delta,
    source_bias_model_true_pi = model_delta_true_pi - required_delta,
    source_bias_model_observed_A = model_delta_observed_A - required_delta,
    outcome_noise_delta = observed_delta - model_delta_observed_A,
    intercept_gap_true_pi = mean(true_pi_weights) - 1,
    intercept_gap_observed_A = mean(observed_weights) - 1,
    weight_mean_treated = if (any(source_fold$A == 1L)) {
      mean(weights[source_fold$A == 1L])
    } else {
      NA_real_
    },
    weight_max = max(weights),
    stringsAsFactors = FALSE
  )
}

append_csv <- function(df, suffix) {
  path <- paste0(prefix, "_", suffix, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("[score-audit] wrote %s (%d rows)\n", path, nrow(df)))
}

summarize_score_rows <- function(rows) {
  groups <- split(rows, interaction(rows$tag, rows$mode, drop = TRUE, sep = "::"))
  do.call(rbind, lapply(groups, function(df) {
    data.frame(
      tag = df$tag[[1L]],
      mode = df$mode[[1L]],
      n_fold_source = nrow(df),
      mean_gamma_intercept_abs = mean(abs(df$gamma_score_intercept)),
      mean_gamma_linf = mean(df$gamma_score_linf),
      mean_alpha_init_intercept_abs = mean(abs(df$alpha_init_score_intercept)),
      mean_alpha_init_linf = mean(df$alpha_init_score_linf),
      mean_alpha_final_gamma_intercept_abs = mean(abs(df$alpha_final_gamma_score_intercept)),
      mean_alpha_final_gamma_linf = mean(df$alpha_final_gamma_score_linf),
      mean_required_delta = mean(df$required_delta),
      mean_observed_delta = mean(df$observed_delta),
      mean_model_delta_true_pi = mean(df$model_delta_true_pi),
      mean_source_bias_observed = mean(df$source_bias_observed),
      mean_source_bias_model_true_pi = mean(df$source_bias_model_true_pi),
      mean_intercept_gap_true_pi = mean(df$intercept_gap_true_pi),
      stringsAsFactors = FALSE
    )
  }))
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
folds <- build_crossfit_folds(split, n_folds)

run_rows <- list()
score_rows <- list()

for (mode in modes) {
  set.seed(seed + if (identical(mode, "one_round")) 200000L else 300000L)
  target_only_ps_cache <- new.env(hash = TRUE, parent = emptyenv())
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
    use_lambda_cache = use_lambda_cache,
    precomputed_folds = folds,
    target_only_ps_cache = target_only_ps_cache,
    nuisance_lambda_rule = "min",
    M_tau = c2_m_tau_train
  )

  run_rows[[length(run_rows) + 1L]] <- data.frame(
    tag = tag,
    mode = mode,
    estimate = as.numeric(res$estimate),
    truth_superpopulation = as.numeric(data$mu1_true),
    bias_vs_superpopulation = as.numeric(res$estimate) - as.numeric(data$mu1_true),
    se = as.numeric(res$se),
    target_estimate = as.numeric(res$target_only$estimate),
    target_bias_vs_superpopulation =
      as.numeric(res$target_only$estimate) - as.numeric(data$mu1_true),
    weight_sum = sum(as.numeric(res$weights)),
    stringsAsFactors = FALSE
  )

  for (k1 in seq_len(n_folds)) {
    target_main <- materialize_fold(folds$target_folds, k1)

    for (source_site in source_sites) {
      src <- res$fold_results[[k1]]$source_results[[source_site]]
      source_main <- materialize_fold(folds$source_folds[[source_site]], k1)

      k2_keys <- names(src$per_k2_alpha)
      k2_values <- as.integer(sub("^k2_", "", k2_keys))

      gamma_scores <- list()
      alpha_init_scores <- list()
      alpha_final_gamma_scores <- list()
      for (k2 in k2_values) {
        k2_key <- paste0("k2_", k2)
        target_k2 <- materialize_fold(folds$target_folds, k2)
        source_k2 <- materialize_fold(folds$source_folds[[source_site]], k2)

        gamma_scores[[k2_key]] <- fold_gamma_score(
          target_k2, source_k2,
          gamma_final = src$gamma_s,
          alpha_init = src$per_k2_alpha[[k2_key]]
        )
        alpha_init_scores[[k2_key]] <- fold_alpha_score(
          source_k2,
          alpha_final = src$alpha_ts,
          gamma_for_weights = src$per_k2_gamma[[k2_key]],
          gamma_truncation = M_TAU_DEFAULT
        )
        alpha_final_gamma_scores[[k2_key]] <- fold_alpha_score(
          source_k2,
          alpha_final = src$alpha_ts,
          gamma_for_weights = src$gamma_s,
          gamma_truncation = M_TAU_DEFAULT
        )
      }

      moment_row <- correction_moments(
        target_fold = target_main,
        source_fold = source_main,
        data = data,
        split = split,
        source_site = source_site,
        gamma_final = src$gamma_s,
        alpha_final = src$alpha_ts
      )

      score_row <- c(
        list(
          tag = tag,
          mode = mode,
          fold = k1,
          source = source_site,
          gamma_lambda = as.numeric(attr(src$gamma_s, "lambda_used") %||% NA_real_),
          alpha_lambda = as.numeric(attr(src$alpha_ts, "lambda_used") %||% NA_real_)
        ),
        score_norm_row("gamma_score", mean_score(gamma_scores)),
        score_norm_row("alpha_init_score", mean_score(alpha_init_scores)),
        score_norm_row("alpha_final_gamma_score", mean_score(alpha_final_gamma_scores)),
        as.list(moment_row[1L, , drop = TRUE])
      )
      score_rows[[length(score_rows) + 1L]] <- as.data.frame(
        score_row,
        stringsAsFactors = FALSE
      )
    }
  }
}

run_df <- do.call(rbind, run_rows)
score_df <- do.call(rbind, score_rows)
summary_df <- summarize_score_rows(score_df)

append_csv(run_df, "run_summary")
append_csv(score_df, "fold_source_scores")
append_csv(summary_df, "score_summary")

cat("[score-audit] summary:\n")
print(summary_df, row.names = FALSE, digits = 5)
cat("[score-audit] DONE\n")
