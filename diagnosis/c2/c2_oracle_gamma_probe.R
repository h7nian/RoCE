#!/usr/bin/env Rscript
# Diagnosis-only probe for oracle density-ratio gamma normalization.
#
# The RoCE DGP stores gamma_params as joint multinomial logits
# log P(R=s_j,A=a|X) / P(R=t|X).  The source-assisted oracle estimator uses
# source-conditional averages of I(A=a) exp(-Z gamma), so this probe compares
# the current raw convention with the source-conditional intercept-normalized
# convention without modifying estimator source.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(sprintf(
    "expected 5 args: tag n_total K p seed; got %d", length(args)
  ), call. = FALSE)
}

tag <- args[[1L]]
n_total <- as.integer(args[[2L]])
K_sites <- as.integer(args[[3L]])
p <- as.integer(args[[4L]])
seed <- as.integer(args[[5L]])
stopifnot(n_total > 0L, K_sites > 0L, p > 0L, seed > 0L)

suppressPackageStartupMessages({
  if (requireNamespace("RoCE", quietly = FALSE)) {
    library(RoCE)
  } else if (requireNamespace("devtools", quietly = FALSE)) {
    devtools::load_all(".")
  } else {
    stop("Neither installed RoCE nor devtools available.", call. = FALSE)
  }
})

roce_function <- function(name) {
  ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("RoCE namespace does not contain required function '%s'.", name),
         call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}

needed_functions <- c(
  "generate_simulation_data", "split_data_by_site", "estimate_oracle_dr",
  "estimate_target_only", "calculate_site_probabilities"
)
for (fn_name in needed_functions) {
  if (!exists(fn_name, mode = "function")) {
    assign(fn_name, roce_function(fn_name))
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path("diagnosis", "c2", "c2_true_gamma_utils.R"))

out_root <- Sys.getenv(
  "C2_ORACLE_GAMMA_OUTPUT_ROOT",
  file.path("diagnosis", "c2", "oracle_gamma_probe")
)
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf(
  "%s_n%d_K%d_p%d_seed%d", tag, n_total, K_sites, p, seed
))

write_csv <- function(df, suffix) {
  path <- paste0(prefix, "_", suffix, ".csv")
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("[oracle-gamma] wrote %s (%d rows)\n", path, nrow(df)))
}

make_calibrated_gamma_params <- function(data, A_val = 1L) {
  calibrated <- data$gamma_params
  for (source_idx in seq_len(data$K)) {
    source_name <- paste0("s", source_idx)
    key <- paste0(source_name, "_", A_val)
    calibrated[[key]] <- c2_true_source_calibration_gamma(
      data, source_name, A_val = A_val
    )
  }
  calibrated
}

oracle_row <- function(label, fit, truth) {
  weights <- fit$weights
  has_weights <- !is.null(weights) && length(weights) > 0L &&
    any(is.finite(as.numeric(weights)))
  data.frame(
    tag = tag,
    seed = seed,
    n_total = n_total,
    K = K_sites,
    p = p,
    method = label,
    estimate = as.numeric(fit$estimate),
    truth = truth,
    bias = as.numeric(fit$estimate) - truth,
    se = as.numeric(fit$se %||% NA_real_),
    ci_lower = as.numeric((fit$ci_lower %||% (fit$estimate - 1.96 * fit$se))),
    ci_upper = as.numeric((fit$ci_upper %||% (fit$estimate + 1.96 * fit$se))),
    covered = truth >= as.numeric((fit$ci_lower %||% (fit$estimate - 1.96 * fit$se))) &&
      truth <= as.numeric((fit$ci_upper %||% (fit$estimate + 1.96 * fit$se))),
    weight_sum = if (has_weights) sum(as.numeric(weights)) else NA_real_,
    weight_min = if (has_weights) min(as.numeric(weights)) else NA_real_,
    weight_max = if (has_weights) max(as.numeric(weights)) else NA_real_,
    stringsAsFactors = FALSE
  )
}

gamma_moment_rows <- function(data, split) {
  Z_true <- as.matrix(data$Z_site_true %||% data$Z_site)
  probs <- calculate_site_probabilities(Z_true, data$gamma_params, data$K)

  rows <- list()
  source_sites <- setdiff(names(split), "t")
  for (source_name in source_sites) {
    source_idx <- c2_source_index(source_name)
    source_rows <- which(data$R == source_name)
    source_prob <- probs[[paste0("p_s", source_idx, "_0")]] +
      probs[[paste0("p_s", source_idx, "_1")]]
    source_pi <- c2_true_source_treatment_propensity(data, source_name, A_val = 1L)

    gamma_raw <- c2_true_source_joint_gamma(data, source_name, A_val = 1L)
    gamma_cal <- c2_true_source_calibration_gamma(data, source_name, A_val = 1L)
    for (label in c("raw_joint", "source_conditional")) {
      gamma <- if (identical(label, "raw_joint")) gamma_raw else gamma_cal
      weights_all <- exp(-as.numeric(cbind(1, Z_true) %*% gamma))
      weights_source <- weights_all[source_rows]
      treated_source <- as.numeric(data$A[source_rows] == 1L)
      rows[[length(rows) + 1L]] <- data.frame(
        tag = tag,
        seed = seed,
        source = source_name,
        gamma_label = label,
        population_intercept = mean(source_prob * source_pi * weights_all) /
          mean(source_prob),
        sample_intercept = mean(treated_source * weights_source),
        sample_treated_weight_mean = if (any(treated_source == 1L)) {
          mean(weights_source[treated_source == 1L])
        } else {
          NA_real_
        },
        sample_weight_max = max(weights_source),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
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
truth <- as.numeric(data$mu1_true)
target_propensity_true <- data$p_treat_true[which(data$R == "t")]

target_only <- estimate_target_only(split, family = "binomial", A_val = 1L)
oracle_raw <- estimate_oracle_dr(
  split, data$alpha1_true, data$gamma_params,
  outcome_type = data$outcome_type,
  target_propensity_true = target_propensity_true
)
oracle_calibrated <- estimate_oracle_dr(
  split, data$alpha1_true, make_calibrated_gamma_params(data, A_val = 1L),
  outcome_type = data$outcome_type,
  target_propensity_true = target_propensity_true
)

write_csv(
  do.call(rbind, list(
    oracle_row("target_only", target_only, truth),
    oracle_row("oracle_raw_joint_gamma", oracle_raw, truth),
    oracle_row("oracle_source_conditional_gamma", oracle_calibrated, truth)
  )),
  "oracle_runs"
)
write_csv(gamma_moment_rows(data, split), "source_gamma_moments")

cat("[oracle-gamma] DONE\n")
