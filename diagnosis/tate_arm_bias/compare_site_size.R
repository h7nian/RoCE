#!/usr/bin/env Rscript
# Per-site sample size test for the TATE bias (HISTORY #0019).
#
# The paper's first-order theory needs s^2 log^2(p) / N -> 0, where N is the
# sample used to FIT the nuisances, i.e. the per-site n, not the total across
# sites. Every production and legacy run so far used n = 1000 per site, so the
# only way to tell a sample-size regime problem from a defect in the method is
# to raise N_SITE and see whether the bias falls roughly like 1/n.
#
# CACHE_FLAG toggles use_lambda_cache (production default TRUE); it is kept as a
# parameter so the same script also runs the paired cache on/off comparison.
#
# Usage: compare_site_size.R CONFIG K SEED CACHE_FLAG N_SITE
#
# N_SITE is required: an optional argument that silently fell back to 1000 once
# turned an n = 4000 batch into a second copy of the n = 1000 batch.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("usage: compare_site_size.R CONFIG K SEED CACHE_FLAG N_SITE", call. = FALSE)
}
configuration <- args[[1L]]
source_count <- as.integer(args[[2L]])
seed <- as.integer(args[[3L]])
use_cache <- identical(args[[4L]], "TRUE")
n_site <- as.integer(args[[5L]])

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
loaded_from <- getNamespaceInfo("RoCE", "path")
if (nzchar(project_library) &&
    !identical(normalizePath(dirname(loaded_from), mustWork = FALSE),
               normalizePath(project_library, mustWork = FALSE))) {
  stop(sprintf("RoCE loaded from %s but ROCE_PROJECT_LIB is %s",
               loaded_from, project_library), call. = FALSE)
}

dimension <- 100L
n_folds <- 10L

set.seed(seed)
data <- RoCE:::generate_simulation_data(
  n_total = n_site * (source_count + 1L), K = source_count, p = dimension,
  config = configuration, estimand_type = "superpopulation",
  outcome_type = "binary", dgp_type = "face", ate_deviation = 0,
  n_deviated_sites = 0L, deviation_mechanism = "treated_arm", warn_ignored = FALSE
)
truth_tate <- as.numeric(data$mu1_true - data$mu0_true)
truth_mu1 <- as.numeric(data$mu1_true)
truth_mu0 <- as.numeric(data$mu0_true)

split <- RoCE:::split_data_by_site(data)
folds <- RoCE:::build_crossfit_folds(split, n_folds)

started <- Sys.time()
fit <- RoCE::run_tate_crossfit(
  split, n_folds = n_folds, communication_mode = "one_round", verbose = FALSE,
  nlambda_init = 100L, family = "binomial", precomputed_folds = folds,
  nuisance_lambda_rule = "min", screening_rule = "soft_penalty",
  use_lambda_cache = use_cache
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))

arm_estimate <- function(name) {
  value <- tryCatch(as.numeric(fit$arm_results[[name]]$estimate), error = function(e) NA_real_)
  if (length(value) != 1L) NA_real_ else value
}
mu1 <- arm_estimate("mu1")
mu0 <- arm_estimate("mu0")

cat(sprintf(
  "RESULT\tcache=%s\t%s\tK%d\tn%d\tseed%d\ttate=%.6f\ttate_bias=%.6f\tmu1_bias=%.6f\tmu0_bias=%.6f\tse=%.6f\tminutes=%.1f\n",
  if (use_cache) "on" else "off", configuration, source_count, n_site, seed,
  as.numeric(fit$estimate), as.numeric(fit$estimate) - truth_tate,
  mu1 - truth_mu1, mu0 - truth_mu0, as.numeric(fit$se), elapsed
))
