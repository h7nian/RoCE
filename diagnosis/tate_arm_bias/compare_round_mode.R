#!/usr/bin/env Rscript
# One-round vs two-round on the same seed.
#
# In one-round mode get_fold_inputs() ignores its `site` argument and hands every
# source the SAME target initial outcome model and target summaries
# (R/cross_fitting_algorithms.R:1380-1387). Any error in those target-derived
# quantities is copied into every source-assisted term and cannot average away
# as sources are added, which matches the observed bias that does not move with
# K. Two-round mode gives each source its own initial model (lines 1313-1320).
#
# Usage: compare_round_mode.R CONFIG K SEED MODE
args <- commandArgs(trailingOnly = TRUE)
configuration <- args[[1L]]
source_count <- as.integer(args[[2L]])
seed <- as.integer(args[[3L]])
mode <- args[[4L]]

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
loaded_from <- getNamespaceInfo("RoCE", "path")
if (!nzchar(project_library) ||
    !identical(normalizePath(dirname(loaded_from), mustWork = FALSE),
               normalizePath(project_library, mustWork = FALSE))) {
  stop(sprintf("RoCE loaded from %s but ROCE_PROJECT_LIB is '%s'",
               loaded_from, project_library), call. = FALSE)
}

set.seed(seed)
simulated <- RoCE:::generate_simulation_data(
  n_total = 1000L * (source_count + 1L), K = source_count, p = 100L,
  config = configuration, estimand_type = "superpopulation",
  outcome_type = "binary", dgp_type = "face", ate_deviation = 0,
  n_deviated_sites = 0L, deviation_mechanism = "treated_arm", warn_ignored = FALSE
)
split <- RoCE:::split_data_by_site(simulated)
folds <- RoCE:::build_crossfit_folds(split, 10L)

started <- Sys.time()
fit <- RoCE::run_tate_crossfit(
  split, n_folds = 10L, communication_mode = mode, verbose = FALSE,
  nlambda_init = 100L, family = "binomial", precomputed_folds = folds,
  nuisance_lambda_rule = "min", screening_rule = "soft_penalty"
)
arm <- function(name) {
  v <- tryCatch(as.numeric(fit$arm_results[[name]]$estimate), error = function(e) NA_real_)
  if (length(v) != 1L) NA_real_ else v
}
mu1_true <- as.numeric(simulated$mu1_true)
mu0_true <- as.numeric(simulated$mu0_true)
cat(sprintf(
  "RESULT\t%s\t%s\tK%d\tseed%d\ttate_bias=%.6f\tmu1_bias=%.6f\tmu0_bias=%.6f\tse=%.6f\tminutes=%.1f\n",
  mode, configuration, source_count, seed,
  as.numeric(fit$estimate) - (mu1_true - mu0_true),
  arm("mu1") - mu1_true, arm("mu0") - mu0_true, as.numeric(fit$se),
  as.numeric(difftime(Sys.time(), started, units = "mins"))
))
