#!/usr/bin/env Rscript
# Does double robustness actually buy anything here?
#
# The augmented estimator's bias is the PRODUCT of the two nuisance errors,
# whereas the plug-in carries the outcome error alone and inverse-probability
# weighting carries the propensity error alone. Feeding all three exactly the
# same cross-fitted nuisance estimates isolates the benefit: if the augmented
# bias is an order of magnitude below the other two, double robustness is doing
# its job and the residual is genuinely second order.
#
# Usage: dr_vs_naive.R CONFIG SEED [N_SITE]
suppressWarnings(suppressPackageStartupMessages(library(glmnet)))
args <- commandArgs(trailingOnly = TRUE)
configuration <- args[[1L]]
seed <- as.integer(args[[2L]])
n_site <- if (length(args) >= 3L) as.integer(args[[3L]]) else 1000L

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

n_folds <- 10L
crossfit_predict <- function(basis, response, folds, subset_mask) {
  prediction <- rep(NA_real_, length(response))
  for (k in seq_len(n_folds)) {
    train <- which(folds != k & subset_mask)
    test <- which(folds == k)
    if (length(train) < 30L || length(unique(response[train])) < 2L) next
    fit <- glmnet::cv.glmnet(basis[train, , drop = FALSE], response[train],
                             family = "binomial", standardize = TRUE, nfolds = 5L)
    prediction[test] <- as.numeric(stats::predict(
      fit, newx = basis[test, , drop = FALSE], s = "lambda.min", type = "response"))
  }
  prediction
}

set.seed(seed)
simulated <- RoCE:::generate_simulation_data(
  n_total = n_site * 3L, K = 2L, p = 100L, config = configuration,
  estimand_type = "superpopulation", outcome_type = "binary", dgp_type = "face",
  ate_deviation = 0, n_deviated_sites = 0L, deviation_mechanism = "treated_arm",
  warn_ignored = FALSE
)
target <- simulated$R == "t"
basis <- as.matrix(simulated$W_outcome)[target, , drop = FALSE]
y <- simulated$Y[target]
a <- simulated$A[target]
set.seed(seed + 1000000L)
folds <- sample(rep_len(seq_len(n_folds), length(y)))

m1 <- crossfit_predict(basis, y, folds, a == 1L)
m0 <- crossfit_predict(basis, y, folds, a == 0L)
e <- pmin(pmax(crossfit_predict(basis, a, folds, rep(TRUE, length(a))), 0.01), 0.99)

truth <- as.numeric(simulated$mu1_true - simulated$mu0_true)
plugin <- mean(m1) - mean(m0)
ipw <- mean(a * y / e) - mean((1 - a) * y / (1 - e))
augmented <- mean(m1 + a * (y - m1) / e) - mean(m0 + (1 - a) * (y - m0) / (1 - e))

cat(sprintf("RESULT\t%s\tn%d\tseed%d\tplugin=%.6f\tipw=%.6f\taipw=%.6f\ttruth=%.6f\n",
            configuration, n_site, seed, plugin - truth, ipw - truth,
            augmented - truth, truth))
