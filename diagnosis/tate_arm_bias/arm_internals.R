#!/usr/bin/env Rscript
# What does borrowing do differently to the control arm?
#
# Production records the aggregation diagnostics only for the treated arm, so
# the control arm's anchor weight, source weights and per-source estimates are
# not on disk. This runs one replicate and dumps both arms' internals side by
# side: the target-only anchor, each source-assisted component, and the weights
# that combine them. Comparing the two arms isolates whether the control arm's
# extra bias enters through the anchor, through the source components, or
# through the weights.
#
# Usage: arm_internals.R CONFIG K SEED
args <- commandArgs(trailingOnly = TRUE)
configuration <- args[[1L]]
source_count <- as.integer(args[[2L]])
seed <- as.integer(args[[3L]])

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

set.seed(seed)
simulated <- RoCE:::generate_simulation_data(
  n_total = 1000L * (source_count + 1L), K = source_count, p = 100L,
  config = configuration, estimand_type = "superpopulation",
  outcome_type = "binary", dgp_type = "face", ate_deviation = 0,
  n_deviated_sites = 0L, deviation_mechanism = "treated_arm", warn_ignored = FALSE
)
split <- RoCE:::split_data_by_site(simulated)
folds <- RoCE:::build_crossfit_folds(split, 10L)

fit <- RoCE::run_tate_crossfit(
  split, n_folds = 10L, communication_mode = "one_round", verbose = FALSE,
  nlambda_init = 100L, family = "binomial", precomputed_folds = folds,
  nuisance_lambda_rule = "min", screening_rule = "soft_penalty"
)

truth <- c(mu1 = as.numeric(simulated$mu1_true), mu0 = as.numeric(simulated$mu0_true))
for (arm in c("mu1", "mu0")) {
  a <- fit$arm_results[[arm]]
  if (is.null(a)) { cat(sprintf("%s: absent\n", arm)); next }
  cat(sprintf("\n===== %s  (truth %.6f) =====\n", arm, truth[[arm]]))
  cat(sprintf("  estimate            %.6f   bias %+.6f\n",
              as.numeric(a$estimate), as.numeric(a$estimate) - truth[[arm]]))
  keys <- names(a)
  cat(sprintf("  fields: %s\n", paste(keys, collapse = ", ")))
  for (k in c("weights", "target_anchor_weight", "source_weights",
              "fold_aggregated_estimates", "target_only_estimate",
              "source_assisted_estimates", "se")) {
    if (!k %in% keys) next
    v <- a[[k]]
    if (is.numeric(v) && length(v) <= 12L) {
      cat(sprintf("  %-26s %s\n", k, paste(sprintf("%.6f", v), collapse = " ")))
    } else if (is.numeric(v)) {
      cat(sprintf("  %-26s n=%d mean %.6f sd %.6f\n", k, length(v), mean(v), stats::sd(v)))
    }
  }
}
