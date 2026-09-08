#!/usr/bin/env Rscript
# Diagnosis-only C2 nuisance refit probe.
#
# For the exact problematic C2 settings, this reruns the production cross-fit
# once, then refits the final source gamma and final source alpha on the same
# folds under fixed lambda multipliers.  The estimator source is untouched.

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

needed_functions <- c(
  "generate_simulation_data", "split_data_by_site", "build_crossfit_folds",
  "run_crossfit", "materialize_fold", "fit_unified_density_ratio",
  "fit_unified_outcome", "calculate_correction_term_cpp", "predict_glm_cpp",
  ".make_plugin_block_design", ".mean_glm_gradient_site_basis"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, roce_function(fn_name))
  }
}

LOGISTIC_CLIP <- roce_constant("LOGISTIC_CLIP", 50)
DIVISION_FLOOR <- roce_constant("DIVISION_FLOOR", 1e-10)
M_TAU_DEFAULT <- roce_constant("M_TAU_DEFAULT", 10)
M_TAU_INFERENCE_DEFAULT <- roce_constant("M_TAU_INFERENCE_DEFAULT", 10)
MAX_ITER_DEFAULT <- roce_constant("MAX_ITER_DEFAULT", 1000)
TOL_DEFAULT <- roce_constant("TOL_DEFAULT", 1e-6)
FAMILY_BINOMIAL <- roce_constant("FAMILY_BINOMIAL", 1L)
LINK_LOGIT <- roce_constant("LINK_LOGIT", 1L)
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("diagnosis", "c2", "c2_true_gamma_utils.R"))

out_dir <- file.path("diagnosis", "c2", "nuisance_refit_probe", "results")
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

ess <- function(w) {
  w <- as.numeric(w)
  if (length(w) == 0L || sum(w^2) <= 0) return(NA_real_)
  sum(w)^2 / sum(w^2)
}

safe_coef_support <- function(x) sum(abs(as.numeric(x)[-1L]) > 1e-8)

make_lambda_specs <- function(lambda_used, lambda_min, lambda_1se) {
  lambda_used <- as.numeric(lambda_used)
  lambda_min <- as.numeric(lambda_min)
  lambda_1se <- as.numeric(lambda_1se)
  base <- lambda_used
  if (!is.finite(base) || base < 0) base <- lambda_min
  if (!is.finite(base) || base < 0) base <- 0
  specs <- data.frame(
    label = c("current", "zero", "x0.01", "x0.1", "x10", "lambda_1se"),
    lambda = c(base, 0, base * 0.01, base * 0.1, base * 10, lambda_1se),
    stringsAsFactors = FALSE
  )
  specs <- specs[is.finite(specs$lambda) & specs$lambda >= 0, , drop = FALSE]
  specs[!duplicated(paste(specs$label, signif(specs$lambda, 12))), , drop = FALSE]
}

source_eval <- function(target_W, source_W, source_Z, source_A, source_Y,
                        gamma, alpha, fold_truth, M_tau_inference,
                        gamma_basis = "fitted") {
  correction <- calculate_correction_term_cpp(
    source_Z, source_A, source_Y, gamma, alpha, source_W, M_tau_inference,
    FAMILY_BINOMIAL, LINK_LOGIT, 1L
  )
  mu_pred <- mean(predict_glm_cpp(target_W, alpha, FAMILY_BINOMIAL, LINK_LOGIT))
  delta <- as.numeric(correction$delta_ts)
  treated <- as.numeric(source_A == 1L)
  logits <- as.numeric(cbind(1, as.matrix(source_Z)) %*% as.numeric(gamma))
  weights <- exp(-pmax(pmin(logits, M_tau_inference), -M_tau_inference))
  active_weights <- weights[treated == 1L]
  data.frame(
    gamma_basis = gamma_basis,
    estimate = mu_pred + delta,
    bias_vs_fold_truth = mu_pred + delta - fold_truth,
    mu_pred = mu_pred,
    mu_pred_bias = mu_pred - fold_truth,
    delta = delta,
    treated_weight_mean = if (length(active_weights)) mean(active_weights) else NA_real_,
    treated_weight_ess = ess(active_weights),
    treated_weight_max = if (length(active_weights)) max(active_weights) else NA_real_,
    clip_any = isTRUE(correction$clip_diagnostics$any_clipped),
    stringsAsFactors = FALSE
  )
}

fit_or_error <- function(expr) {
  tryCatch(
    list(ok = TRUE, value = eval.parent(substitute(expr)), error = NA_character_),
    error = function(e) list(ok = FALSE, value = NULL, error = conditionMessage(e))
  )
}

write_csv <- function(df, suffix) {
  path <- paste0(prefix, "_", suffix, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("[nuisance-refit] wrote %s (%d rows)\n", path, nrow(df)))
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
target <- split[["t"]]
truth <- as.numeric(data$mu1_true)
folds <- build_crossfit_folds(split, n_folds)

run_rows <- list()
refit_rows <- list()
eval_rows <- list()
oracle_rows <- list()

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
    bias_vs_superpopulation = res$estimate - truth,
    se = res$se,
    covered = truth >= res$ci_lower && truth <= res$ci_upper,
    target_estimate = res$target_only$estimate,
    target_bias_vs_superpopulation = res$target_only$estimate - truth,
    weight_sum = sum(as.numeric(res$weights)),
    weight_min = min(as.numeric(res$weights)),
    weight_max = max(as.numeric(res$weights)),
    stringsAsFactors = FALSE
  )

  for (k1 in seq_len(n_folds)) {
    target_main <- materialize_fold(folds$target_folds, k1)
    target_true_W <- target$W_outcome_true[target_main$original_idx, , drop = FALSE]
    fold_truth <- mean(predict_alpha(target_true_W, data$alpha1_true))

    for (s in source_sites) {
      src <- res$fold_results[[k1]]$source_results[[s]]
      source_main <- materialize_fold(folds$source_folds[[s]], k1)
      source_true_W <- split[[s]]$W_outcome_true[source_main$original_idx, , drop = FALSE]
      source_true_Z <- split[[s]]$Z_site_true[source_main$original_idx, , drop = FALSE]

      k2_keys <- names(src$per_k2_alpha)
      k2_values <- as.integer(sub("^k2_", "", k2_keys))
      source_calib_list <- lapply(k2_values, function(k2) {
        materialize_fold(folds$source_folds[[s]], k2)
      })
      names(source_calib_list) <- k2_keys

      Z_cal_stack <- stack_field(source_calib_list, "Z_site")
      W_cal_stack <- stack_field(source_calib_list, "W_outcome")
      A_cal_stack <- stack_field(source_calib_list, "A")
      Y_cal_stack <- stack_field(source_calib_list, "Y")

      mean_grad_list <- lapply(k2_values, function(k2) {
        key <- paste0("k2_", k2)
        target_k2 <- materialize_fold(folds$target_folds, k2)
        .mean_glm_gradient_site_basis(
          target_k2$W_outcome, target_k2$Z_site,
          src$per_k2_alpha[[key]], FAMILY_BINOMIAL, LINK_LOGIT
        )
      })
      mean_grad_avg <- mean_list(mean_grad_list)

      W_plugin_block <- .make_plugin_block_design(
        lapply(source_calib_list, `[[`, "W_outcome"),
        caller = "c2_nuisance_refit_probe(gamma)"
      )
      alpha_plugin_block <- c(0, unlist(src$per_k2_alpha, use.names = FALSE))

      Z_plugin_block <- .make_plugin_block_design(
        lapply(source_calib_list, `[[`, "Z_site"),
        caller = "c2_nuisance_refit_probe(alpha)"
      )
      gamma_plugin_block <- c(0, unlist(src$per_k2_gamma, use.names = FALSE))

      gamma_specs <- make_lambda_specs(
        attr(src$gamma_s, "lambda_used"),
        attr(src$gamma_s, "lambda_min"),
        attr(src$gamma_s, "lambda_1se")
      )
      alpha_specs <- make_lambda_specs(
        attr(src$alpha_ts, "lambda_used"),
        attr(src$alpha_ts, "lambda_min"),
        attr(src$alpha_ts, "lambda_1se")
      )

      gamma_fits <- list(current = src$gamma_s)
      for (i in seq_len(nrow(gamma_specs))) {
        lab <- gamma_specs$label[i]
        if (identical(lab, "current")) next
        fit <- fit_or_error(fit_unified_density_ratio(
          Z_site = Z_cal_stack, A = A_cal_stack,
          mean_grad_psi = mean_grad_avg, alpha_init = alpha_plugin_block,
          lambda = gamma_specs$lambda[i],
          calibrated = TRUE, M_tau = M_TAU_DEFAULT,
          W_outcome = W_plugin_block, A_val = 1L,
          family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
        ))
        refit_rows[[length(refit_rows) + 1L]] <- data.frame(
          tag = tag, mode = mode, fold = k1, source = s,
          nuisance = "gamma", label = lab, lambda = gamma_specs$lambda[i],
          ok = fit$ok, error = fit$error,
          support = if (fit$ok) safe_coef_support(fit$value) else NA_integer_,
          l2 = if (fit$ok) sqrt(mean(as.numeric(fit$value)^2)) else NA_real_,
          stringsAsFactors = FALSE
        )
        if (fit$ok) gamma_fits[[lab]] <- fit$value
      }
      refit_rows[[length(refit_rows) + 1L]] <- data.frame(
        tag = tag, mode = mode, fold = k1, source = s,
        nuisance = "gamma", label = "current",
        lambda = as.numeric(attr(src$gamma_s, "lambda_used") %||% NA_real_),
        ok = TRUE, error = NA_character_,
        support = safe_coef_support(src$gamma_s),
        l2 = sqrt(mean(as.numeric(src$gamma_s)^2)),
        stringsAsFactors = FALSE
      )

      alpha_fits <- list(current = src$alpha_ts)
      for (i in seq_len(nrow(alpha_specs))) {
        lab <- alpha_specs$label[i]
        if (identical(lab, "current")) next
        fit <- fit_or_error(fit_unified_outcome(
          W_outcome = W_cal_stack, Y = Y_cal_stack, A = A_cal_stack,
          A_val = 1L, gamma_s = gamma_plugin_block,
          lambda = alpha_specs$lambda[i],
          family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT,
          calibrated = TRUE, M_tau = M_TAU_DEFAULT, Z_site = Z_plugin_block
        ))
        refit_rows[[length(refit_rows) + 1L]] <- data.frame(
          tag = tag, mode = mode, fold = k1, source = s,
          nuisance = "alpha", label = lab, lambda = alpha_specs$lambda[i],
          ok = fit$ok, error = fit$error,
          support = if (fit$ok) safe_coef_support(fit$value) else NA_integer_,
          l2 = if (fit$ok) sqrt(mean(as.numeric(fit$value)^2)) else NA_real_,
          stringsAsFactors = FALSE
        )
        if (fit$ok) alpha_fits[[lab]] <- fit$value
      }
      refit_rows[[length(refit_rows) + 1L]] <- data.frame(
        tag = tag, mode = mode, fold = k1, source = s,
        nuisance = "alpha", label = "current",
        lambda = as.numeric(attr(src$alpha_ts, "lambda_used") %||% NA_real_),
        ok = TRUE, error = NA_character_,
        support = safe_coef_support(src$alpha_ts),
        l2 = sqrt(mean(as.numeric(src$alpha_ts)^2)),
        stringsAsFactors = FALSE
      )

      gamma_labels <- intersect(c("current", "zero", "x0.01", "x0.1", "x10", "lambda_1se"),
                                names(gamma_fits))
      alpha_labels <- intersect(c("current", "zero", "x0.01", "x0.1", "x10", "lambda_1se"),
                                names(alpha_fits))
      for (g_lab in gamma_labels) {
        for (a_lab in alpha_labels) {
          ev <- source_eval(
            target_main$W_outcome, source_main$W_outcome,
            source_main$Z_site, source_main$A, source_main$Y,
            gamma_fits[[g_lab]], alpha_fits[[a_lab]], fold_truth,
            M_TAU_INFERENCE_DEFAULT, gamma_basis = "fitted"
          )
          ev$tag <- tag
          ev$mode <- mode
          ev$fold <- k1
          ev$source <- s
          ev$gamma_label <- g_lab
          ev$alpha_label <- a_lab
          ev$fold_truth <- fold_truth
          eval_rows[[length(eval_rows) + 1L]] <- ev
        }
      }

      true_gamma <- c2_true_source_calibration_gamma(data, s, A_val = 1L)
      if (!is.null(true_gamma)) {
        true_alpha_ev <- source_eval(
          target_true_W, source_true_W, source_true_Z,
          source_main$A, source_main$Y,
          true_gamma, data$alpha1_true, fold_truth,
          M_TAU_INFERENCE_DEFAULT, gamma_basis = "true"
        )
        true_alpha_ev$variant <- "true_gamma_true_alpha"
        oracle_rows[[length(oracle_rows) + 1L]] <- true_alpha_ev

        true_gamma_est_alpha_ev <- source_eval(
          target_main$W_outcome, source_main$W_outcome, source_true_Z,
          source_main$A, source_main$Y,
          true_gamma, src$alpha_ts, fold_truth,
          M_TAU_INFERENCE_DEFAULT, gamma_basis = "true"
        )
        true_gamma_est_alpha_ev$variant <- "true_gamma_estimated_alpha"
        oracle_rows[[length(oracle_rows) + 1L]] <- true_gamma_est_alpha_ev

        est_gamma_true_alpha_ev <- source_eval(
          target_true_W, source_true_W, source_main$Z_site,
          source_main$A, source_main$Y,
          src$gamma_s, data$alpha1_true, fold_truth,
          M_TAU_INFERENCE_DEFAULT, gamma_basis = "fitted"
        )
        est_gamma_true_alpha_ev$variant <- "estimated_gamma_true_alpha"
        oracle_rows[[length(oracle_rows) + 1L]] <- est_gamma_true_alpha_ev

        last <- length(oracle_rows)
        for (idx in (last - 2L):last) {
          oracle_rows[[idx]]$tag <- tag
          oracle_rows[[idx]]$mode <- mode
          oracle_rows[[idx]]$fold <- k1
          oracle_rows[[idx]]$source <- s
          oracle_rows[[idx]]$fold_truth <- fold_truth
        }
      }
    }
  }
}

write_csv(do.call(rbind, run_rows), "run_summary")
write_csv(do.call(rbind, refit_rows), "nuisance_refits")
write_csv(do.call(rbind, eval_rows), "source_eval_grid")
write_csv(do.call(rbind, oracle_rows), "oracle_checks")

cat("[nuisance-refit] DONE\n")
