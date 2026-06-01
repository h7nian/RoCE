#!/usr/bin/env Rscript
# Diagnosis-only C2 residual-balance probe.
#
# This probe decomposes each source-assisted fold estimate into:
#   required_delta = E_t[m_true(X) - m_hat(X)]
#   observed_delta = E_s[I(A=1) w_hat(X) {Y - m_hat(X)}]
#   model_delta    = E_s[I(A=1) w_hat(X) {m_true(X) - m_hat(X)}]
#   pi_delta       = E_s[pi_true(X) w_hat(X) {m_true(X) - m_hat(X)}]
#
# It is intentionally read-only with respect to estimator source code.

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
  "run_crossfit", "materialize_fold", "fit_unified_density_ratio",
  "fit_unified_outcome", "predict_glm_cpp", ".make_plugin_block_design",
  ".mean_glm_gradient_site_basis"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, facehd_function(fn_name))
  }
}

dr_cv_scale_patch <- Sys.getenv("C2_DR_CV_SCALE_PATCH", "none")
dr_cv_train_scale <- Sys.getenv("C2_DR_CV_TRAIN_SCALE", "production")
if (!dr_cv_scale_patch %in% c("none", "validation_source_scale")) {
  stop("C2_DR_CV_SCALE_PATCH must be 'none' or 'validation_source_scale'.",
       call. = FALSE)
}
if (dr_cv_scale_patch != "none") {
  source(file.path("diagnosis", "c2", "c2_dr_cv_scale_patch.R"))
  patch_info <- install_c2_dr_cv_scale_patch(train_scale = dr_cv_train_scale)
  fit_initial_density_ratio <- facehd_function("fit_initial_density_ratio")
  fit_unified_density_ratio <- facehd_function("fit_unified_density_ratio")
  cat(sprintf(
    "[residual-balance] installed diagnosis-only density-ratio CV scale patch: validation=%s train_scale=%s cpp=%s\n",
    dr_cv_scale_patch, patch_info$train_scale, patch_info$cpp_path
  ))
}

LOGISTIC_CLIP <- facehd_constant("LOGISTIC_CLIP", 50)
M_TAU_DEFAULT <- facehd_constant("M_TAU_DEFAULT", 10)
M_TAU_INFERENCE_DEFAULT <- facehd_constant("M_TAU_INFERENCE_DEFAULT", 10)
FAMILY_BINOMIAL <- facehd_constant("FAMILY_BINOMIAL", 1L)
LINK_LOGIT <- facehd_constant("LINK_LOGIT", 1L)
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("diagnosis", "c2", "c2_true_gamma_utils.R"))

out_root <- Sys.getenv(
  "C2_RBAL_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "residual_balance_probe")
)
out_dir <- file.path(out_root, "results")
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

stack_field <- function(fold_data_list, field) {
  values <- lapply(fold_data_list, `[[`, field)
  if (is.matrix(values[[1L]])) {
    do.call(rbind, values)
  } else {
    unlist(values, use.names = FALSE)
  }
}

mean_list <- function(values) Reduce(`+`, values) / length(values)

rmse <- function(x) sqrt(mean(as.numeric(x)^2))

weighted_col_means <- function(X, weights) {
  weights <- as.numeric(weights)
  colMeans(sweep(as.matrix(X), 1L, weights, "*"))
}

normalized_weighted_col_means <- function(X, weights) {
  weights <- as.numeric(weights)
  denom <- sum(weights)
  if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, ncol(as.matrix(X))))
  colSums(sweep(as.matrix(X), 1L, weights, "*")) / denom
}

safe_cor <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 2L || stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
  as.numeric(stats::cor(x, y))
}

make_lambda_specs <- function(lambda_used, lambda_1se, labels) {
  base <- as.numeric(lambda_used)
  if (!is.finite(base) || base < 0) base <- 0
  lambda_1se <- as.numeric(lambda_1se)
  values <- c(
    current = base,
    zero = 0,
    x0.1 = base * 0.1,
    x10 = base * 10,
    lambda_1se = lambda_1se
  )
  keep <- intersect(labels, names(values))
  values <- values[keep]
  values <- values[is.finite(values) & values >= 0]
  values[!duplicated(paste(names(values), signif(values, 12)))]
}

fit_or_error <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = eval.parent(substitute(expr)), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

fit_refit_nuisances <- function(src, source_calib_list, mean_grad_avg) {
  Z_cal_stack <- stack_field(source_calib_list, "Z_site")
  W_cal_stack <- stack_field(source_calib_list, "W_outcome")
  A_cal_stack <- stack_field(source_calib_list, "A")
  Y_cal_stack <- stack_field(source_calib_list, "Y")

  W_plugin_block <- .make_plugin_block_design(
    lapply(source_calib_list, `[[`, "W_outcome"),
    caller = "c2_residual_balance_probe(gamma)"
  )
  alpha_plugin_block <- c(0, unlist(src$per_k2_alpha, use.names = FALSE))

  Z_plugin_block <- .make_plugin_block_design(
    lapply(source_calib_list, `[[`, "Z_site"),
    caller = "c2_residual_balance_probe(alpha)"
  )
  gamma_plugin_block <- c(0, unlist(src$per_k2_gamma, use.names = FALSE))

  gamma_fits <- list(current = src$gamma_s)
  gamma_errors <- list()
  gamma_lambdas <- make_lambda_specs(
    attr(src$gamma_s, "lambda_used"),
    attr(src$gamma_s, "lambda_1se"),
    labels = c("zero", "x0.1", "x10")
  )
  for (label in names(gamma_lambdas)) {
    if (identical(label, "current")) next
    fit <- fit_or_error(fit_unified_density_ratio(
      Z_site = Z_cal_stack, A = A_cal_stack,
      mean_grad_psi = mean_grad_avg, alpha_init = alpha_plugin_block,
      lambda = gamma_lambdas[[label]],
      calibrated = TRUE, M_tau = M_TAU_DEFAULT,
      W_outcome = W_plugin_block, A_val = 1L,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    ))
    if (fit$ok) {
      gamma_fits[[label]] <- fit$value
    } else {
      gamma_errors[[label]] <- fit$error
    }
  }

  alpha_fits <- list(current = src$alpha_ts)
  alpha_errors <- list()
  alpha_lambdas <- make_lambda_specs(
    attr(src$alpha_ts, "lambda_used"),
    attr(src$alpha_ts, "lambda_1se"),
    labels = c("zero", "lambda_1se")
  )
  for (label in names(alpha_lambdas)) {
    if (identical(label, "current")) next
    fit <- fit_or_error(fit_unified_outcome(
      W_outcome = W_cal_stack, Y = Y_cal_stack, A = A_cal_stack,
      A_val = 1L, gamma_s = gamma_plugin_block,
      lambda = alpha_lambdas[[label]],
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT, Z_site = Z_plugin_block
    ))
    if (fit$ok) {
      alpha_fits[[label]] <- fit$value
    } else {
      alpha_errors[[label]] <- fit$error
    }
  }

  list(
    gamma_fits = gamma_fits,
    alpha_fits = alpha_fits,
    gamma_errors = gamma_errors,
    alpha_errors = alpha_errors
  )
}

evaluate_residual_balance <- function(target_W, source_W, target_Z, source_Z,
                                      target_true_m, source_true_m,
                                      source_A, source_Y, source_pi,
                                      gamma, alpha,
                                      gamma_label, alpha_label,
                                      gamma_basis, alpha_basis) {
  target_m_hat <- predict_alpha(target_W, alpha)
  source_m_hat <- predict_alpha(source_W, alpha)

  logits <- as.numeric(cbind(1, as.matrix(source_Z)) %*% as.numeric(gamma))
  weights <- exp(-pmax(pmin(logits, M_TAU_INFERENCE_DEFAULT), -M_TAU_INFERENCE_DEFAULT))
  treated <- as.numeric(source_A == 1L)
  observed_weights <- treated * weights
  pi_weights <- as.numeric(source_pi) * weights

  source_residual_y <- as.numeric(source_Y) - source_m_hat
  source_residual_true <- as.numeric(source_true_m) - source_m_hat
  target_residual_true <- as.numeric(target_true_m) - target_m_hat

  required_delta <- mean(target_residual_true)
  observed_delta <- mean(observed_weights * source_residual_y)
  model_delta_observed_A <- mean(observed_weights * source_residual_true)
  model_delta_true_pi <- mean(pi_weights * source_residual_true)

  target_z_mean <- colMeans(as.matrix(target_Z))
  raw_z_observed_A <- weighted_col_means(source_Z, observed_weights)
  raw_z_true_pi <- weighted_col_means(source_Z, pi_weights)
  norm_z_observed_A <- normalized_weighted_col_means(source_Z, observed_weights)
  norm_z_true_pi <- normalized_weighted_col_means(source_Z, pi_weights)

  data.frame(
    gamma_label = gamma_label,
    alpha_label = alpha_label,
    gamma_basis = gamma_basis,
    alpha_basis = alpha_basis,
    required_delta = required_delta,
    observed_delta = observed_delta,
    model_delta_observed_A = model_delta_observed_A,
    model_delta_true_pi = model_delta_true_pi,
    source_bias_observed = observed_delta - required_delta,
    source_bias_model_observed_A = model_delta_observed_A - required_delta,
    source_bias_model_true_pi = model_delta_true_pi - required_delta,
    outcome_noise_delta = observed_delta - model_delta_observed_A,
    target_m_hat_mean = mean(target_m_hat),
    target_m_true_mean = mean(target_true_m),
    source_m_hat_mean = mean(source_m_hat),
    source_m_true_mean = mean(source_true_m),
    target_residual_rmse = rmse(target_residual_true),
    source_residual_rmse = rmse(source_residual_true),
    intercept_gap_observed_A = mean(observed_weights) - 1,
    intercept_gap_true_pi = mean(pi_weights) - 1,
    raw_z_gap_observed_A = rmse(raw_z_observed_A - target_z_mean),
    raw_z_gap_true_pi = rmse(raw_z_true_pi - target_z_mean),
    norm_z_gap_observed_A = rmse(norm_z_observed_A - target_z_mean),
    norm_z_gap_true_pi = rmse(norm_z_true_pi - target_z_mean),
    residual_weight_cor_observed_A = safe_cor(source_residual_true, observed_weights),
    residual_weight_cor_true_pi = safe_cor(source_residual_true, pi_weights),
    weight_max = max(weights),
    treated_weight_mean = if (any(treated == 1L)) mean(weights[treated == 1L]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

write_csv <- function(df, suffix) {
  path <- paste0(prefix, "_", suffix, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("[residual-balance] wrote %s (%d rows)\n", path, nrow(df)))
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
  dgp_type = "facehd",
  warn_ignored = FALSE
)
split <- split_data_by_site(data)
source_sites <- setdiff(names(split), "t")
target <- split[["t"]]
folds <- build_crossfit_folds(split, n_folds)

target_global_idx <- which(data$R == "t")
source_global_idx <- lapply(source_sites, function(s) which(data$R == s))
names(source_global_idx) <- source_sites
source_true_pi <- lapply(source_sites, function(s) {
  c2_true_source_treatment_propensity(data, s, A_val = 1L)
})
names(source_true_pi) <- source_sites

run_rows <- list()
balance_rows <- list()
error_rows <- list()

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
    mode = mode,
    estimate = res$estimate,
    truth_superpopulation = as.numeric(data$mu1_true),
    bias_vs_superpopulation = res$estimate - as.numeric(data$mu1_true),
    se = res$se,
    target_estimate = res$target_only$estimate,
    target_bias_vs_superpopulation =
      res$target_only$estimate - as.numeric(data$mu1_true),
    weight_sum = sum(as.numeric(res$weights)),
    stringsAsFactors = FALSE
  )

  for (k1 in seq_len(n_folds)) {
    target_main <- materialize_fold(folds$target_folds, k1)
    target_true_W <- target$W_outcome_true[target_main$original_idx, , drop = FALSE]
    target_true_Z <- target$Z_site_true[target_main$original_idx, , drop = FALSE]
    target_true_m <- predict_alpha(target_true_W, data$alpha1_true)

    for (s in source_sites) {
      src <- res$fold_results[[k1]]$source_results[[s]]
      source_main <- materialize_fold(folds$source_folds[[s]], k1)
      source_true_W <- split[[s]]$W_outcome_true[source_main$original_idx, , drop = FALSE]
      source_true_Z <- split[[s]]$Z_site_true[source_main$original_idx, , drop = FALSE]
      source_true_m <- predict_alpha(source_true_W, data$alpha1_true)
      source_pi <- source_true_pi[[s]][source_global_idx[[s]][source_main$original_idx]]

      k2_keys <- names(src$per_k2_alpha)
      k2_values <- as.integer(sub("^k2_", "", k2_keys))
      source_calib_list <- lapply(k2_values, function(k2) {
        materialize_fold(folds$source_folds[[s]], k2)
      })
      names(source_calib_list) <- k2_keys

      mean_grad_list <- lapply(k2_values, function(k2) {
        key <- paste0("k2_", k2)
        target_k2 <- materialize_fold(folds$target_folds, k2)
        .mean_glm_gradient_site_basis(
          target_k2$W_outcome, target_k2$Z_site,
          src$per_k2_alpha[[key]], FAMILY_BINOMIAL, LINK_LOGIT
        )
      })
      refits <- fit_refit_nuisances(src, source_calib_list, mean_list(mean_grad_list))

      if (length(refits$gamma_errors) > 0L) {
        for (label in names(refits$gamma_errors)) {
          error_rows[[length(error_rows) + 1L]] <- data.frame(
            tag = tag, mode = mode, fold = k1, source = s,
            nuisance = "gamma", label = label, error = refits$gamma_errors[[label]],
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(refits$alpha_errors) > 0L) {
        for (label in names(refits$alpha_errors)) {
          error_rows[[length(error_rows) + 1L]] <- data.frame(
            tag = tag, mode = mode, fold = k1, source = s,
            nuisance = "alpha", label = label, error = refits$alpha_errors[[label]],
            stringsAsFactors = FALSE
          )
        }
      }

      gamma_fits <- refits$gamma_fits
      alpha_fits <- refits$alpha_fits
      combo_labels <- expand.grid(
        gamma_label = names(gamma_fits),
        alpha_label = names(alpha_fits),
        stringsAsFactors = FALSE
      )
      for (idx in seq_len(nrow(combo_labels))) {
        gamma_label <- combo_labels$gamma_label[[idx]]
        alpha_label <- combo_labels$alpha_label[[idx]]
        row <- evaluate_residual_balance(
          target_W = target_main$W_outcome,
          source_W = source_main$W_outcome,
          target_Z = target_main$Z_site,
          source_Z = source_main$Z_site,
          target_true_m = target_true_m,
          source_true_m = source_true_m,
          source_A = source_main$A,
          source_Y = source_main$Y,
          source_pi = source_pi,
          gamma = gamma_fits[[gamma_label]],
          alpha = alpha_fits[[alpha_label]],
          gamma_label = gamma_label,
          alpha_label = alpha_label,
          gamma_basis = "fitted",
          alpha_basis = "fitted"
        )
        row$tag <- tag
        row$mode <- mode
        row$fold <- k1
        row$source <- s
        row$variant <- "refit"
        balance_rows[[length(balance_rows) + 1L]] <- row
      }

      raw_true_gamma <- c2_true_source_joint_gamma(data, s, A_val = 1L)
      true_gamma <- c2_true_source_calibration_gamma(data, s, A_val = 1L)
      if (!is.null(true_gamma)) {
        oracle_specs <- list(
          raw_gamma_current_alpha = list(
            gamma = raw_true_gamma, alpha = src$alpha_ts,
            target_W = target_main$W_outcome, source_W = source_main$W_outcome,
            target_Z = target_true_Z, source_Z = source_true_Z,
            gamma_label = "raw_true", alpha_label = "current",
            gamma_basis = "true", alpha_basis = "fitted"
          ),
          raw_gamma_true_alpha = list(
            gamma = raw_true_gamma, alpha = data$alpha1_true,
            target_W = target_true_W, source_W = source_true_W,
            target_Z = target_true_Z, source_Z = source_true_Z,
            gamma_label = "raw_true", alpha_label = "true",
            gamma_basis = "true", alpha_basis = "true"
          ),
          true_gamma_current_alpha = list(
            gamma = true_gamma, alpha = src$alpha_ts,
            target_W = target_main$W_outcome, source_W = source_main$W_outcome,
            target_Z = target_true_Z, source_Z = source_true_Z,
            gamma_label = "true", alpha_label = "current",
            gamma_basis = "true", alpha_basis = "fitted"
          ),
          current_gamma_true_alpha = list(
            gamma = src$gamma_s, alpha = data$alpha1_true,
            target_W = target_true_W, source_W = source_true_W,
            target_Z = target_main$Z_site, source_Z = source_main$Z_site,
            gamma_label = "current", alpha_label = "true",
            gamma_basis = "fitted", alpha_basis = "true"
          ),
          true_gamma_true_alpha = list(
            gamma = true_gamma, alpha = data$alpha1_true,
            target_W = target_true_W, source_W = source_true_W,
            target_Z = target_true_Z, source_Z = source_true_Z,
            gamma_label = "true", alpha_label = "true",
            gamma_basis = "true", alpha_basis = "true"
          )
        )
        for (spec_name in names(oracle_specs)) {
          spec <- oracle_specs[[spec_name]]
          row <- evaluate_residual_balance(
            target_W = spec$target_W,
            source_W = spec$source_W,
            target_Z = spec$target_Z,
            source_Z = spec$source_Z,
            target_true_m = target_true_m,
            source_true_m = source_true_m,
            source_A = source_main$A,
            source_Y = source_main$Y,
            source_pi = source_pi,
            gamma = spec$gamma,
            alpha = spec$alpha,
            gamma_label = spec$gamma_label,
            alpha_label = spec$alpha_label,
            gamma_basis = spec$gamma_basis,
            alpha_basis = spec$alpha_basis
          )
          row$tag <- tag
          row$mode <- mode
          row$fold <- k1
          row$source <- s
          row$variant <- spec_name
          balance_rows[[length(balance_rows) + 1L]] <- row
        }
      }
    }
  }
}

write_csv(do.call(rbind, run_rows), "run_summary")
write_csv(do.call(rbind, balance_rows), "residual_balance")
if (length(error_rows) > 0L) {
  write_csv(do.call(rbind, error_rows), "refit_errors")
} else {
  write_csv(data.frame(
    tag = character(0), mode = character(0), fold = integer(0),
    source = character(0), nuisance = character(0), label = character(0),
    error = character(0)
  ), "refit_errors")
}

cat("[residual-balance] DONE\n")
