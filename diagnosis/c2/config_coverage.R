#!/usr/bin/env Rscript
# General config coverage probe (regression check for the density-ratio CV-scale fix).
#
# Diagnosis-only. Runs the production cross-fitting estimator for a given config
# (C1/C2/C3/C4) and records the final estimate, se, bias, and CI coverage. Used to
# confirm the cv_utils.h density-ratio CV-scale fix does not regress configs other
# than C2. Mirrors the data-generation + run_crossfit path of c2_score_moment_audit.R;
# does NOT compute C2-specific moments (so it is config-agnostic).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(sprintf(
    "expected 9 args: tag config n_total K p seed n_folds n_cores modes; got %d",
    length(args)
  ), call. = FALSE)
}

tag <- args[[1L]]
config <- args[[2L]]
n_total <- as.integer(args[[3L]])
K_sites <- as.integer(args[[4L]])
p <- as.integer(args[[5L]])
seed <- as.integer(args[[6L]])
n_folds <- as.integer(args[[7L]])
n_cores <- as.integer(args[[8L]])
modes <- strsplit(args[[9L]], ",", fixed = TRUE)[[1L]]
modes <- modes[nzchar(modes)]
nlambda_init <- as.integer(Sys.getenv("C2_CFG_NLAMBDA_INIT", "100"))
use_lambda_cache <- tolower(Sys.getenv("C2_CFG_USE_LAMBDA_CACHE", "true")) %in%
  c("1", "true", "yes")

stopifnot(n_total > 0L, K_sites > 0L, p > 0L, seed > 0L, n_folds >= 3L, n_cores >= 1L)
if (!(config %in% c("C1", "C2", "C3", "C4"))) {
  stop("config must be one of C1/C2/C3/C4.", call. = FALSE)
}
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

roce_function <- function(name) {
  ns <- asNamespace("RoCE")
  if (!exists(name, envir = ns, inherits = FALSE)) {
    stop(sprintf("RoCE namespace does not contain '%s'.", name), call. = FALSE)
  }
  get(name, envir = ns, inherits = FALSE)
}
for (fn in c("generate_simulation_data", "split_data_by_site",
             "build_crossfit_folds", "run_crossfit")) {
  if (!exists(fn, mode = "function")) assign(fn, roce_function(fn))
}

out_root <- Sys.getenv("C2_CFG_OUTPUT_ROOT",
                       file.path("diagnosis", "c2", "config_coverage"))
out_dir <- file.path(out_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf("%s_%s_n%d_K%d_p%d_seed%d_kf%d",
                                     tag, config, n_total, K_sites, p, seed, n_folds))

set.seed(seed)
data <- generate_simulation_data(
  n_total = n_total, K = K_sites, p = p, config = config,
  estimand_type = "superpopulation", site_allocation = "model",
  transform_type = "mild", outcome_type = "binary",
  heterogeneity_type = "none", shift_strength = 0.5,
  dgp_type = "roce", warn_ignored = FALSE
)
split <- split_data_by_site(data)
folds <- build_crossfit_folds(split, n_folds)
truth <- as.numeric(data$mu1_true)

rows <- list()
for (mode in modes) {
  set.seed(seed + if (identical(mode, "one_round")) 200000L else 300000L)
  res <- run_crossfit(
    split, n_folds = n_folds, communication_mode = mode,
    lambda_selection = "cv", lambda_rule = "min", verbose = FALSE,
    n_cores = n_cores, nlambda_init = nlambda_init, family = "binomial",
    A_val = 1L, use_lambda_cache = use_lambda_cache,
    precomputed_folds = folds, nuisance_lambda_rule = "min"
  )
  est <- as.numeric(res$estimate); se <- as.numeric(res$se)
  rows[[length(rows) + 1L]] <- data.frame(
    tag = tag, config = config, mode = mode,
    n_total = n_total, K = K_sites, p = p, seed = seed,
    estimate = est, truth = truth, bias = est - truth, se = se,
    covered = abs(est - truth) <= 1.96 * se,
    target_estimate = as.numeric(res$target_only$estimate),
    target_bias = as.numeric(res$target_only$estimate) - truth,
    stringsAsFactors = FALSE
  )
}
df <- do.call(rbind, rows)
write.csv(df, paste0(prefix, "_coverage.csv"), row.names = FALSE)
cat(sprintf("[config-cov] %s config=%s n=%d K=%d wrote %d rows\n",
            tag, config, n_total, K_sites, nrow(df)))
print(df[, c("config", "mode", "bias", "se", "covered")], row.names = FALSE)
cat("[config-cov] DONE\n")
