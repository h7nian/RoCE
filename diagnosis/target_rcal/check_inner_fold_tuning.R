#!/usr/bin/env Rscript
# A held-out inner fold may not affect its nuisance fit through a reused penalty.
# This checks that property separately from outer-fold independence.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("usage: check_inner_fold_tuning.R OUTPUT_CSV")
.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
set.seed(314)
data <- RoCE:::generate_simulation_data(
  n_total = 2000L, K = 1L, p = 4L, config = "C1",
  estimand_type = "superpopulation", outcome_type = "binary", dgp_type = "face",
  ate_deviation = 0, n_deviated_sites = 0L, warn_ignored = FALSE
)
original <- RoCE:::split_data_by_site(data)
original_folds <- RoCE:::build_crossfit_folds(original, 3L)
rows <- lapply(c("one_round", "two_round"), function(communication_mode) {
  nuisance_site <- if (communication_mode == "one_round") "t" else "s1"
  nuisance_folds <- if (nuisance_site == "t") original_folds$target_folds else
    original_folds$source_folds[[nuisance_site]]
  changed <- original
  held_out <- nuisance_folds[[3L]]$original_idx
  changed[[nuisance_site]]$Y[held_out] <-
    RoCE:::with_seed(719, sample(changed[[nuisance_site]]$Y[held_out]))
  changed_folds <- RoCE:::build_crossfit_folds(changed, 3L)
  changed_nuisance_folds <- if (nuisance_site == "t") changed_folds$target_folds else
    changed_folds$source_folds[[nuisance_site]]
  stopifnot(identical(lapply(nuisance_folds, `[[`, "original_idx"),
                      lapply(changed_nuisance_folds, `[[`, "original_idx")))
  comparisons <- lapply(c(TRUE, FALSE), function(use_cache) {
    fit <- function(split, folds) RoCE:::with_seed(501, RoCE::run_crossfit(
      split, n_folds = 3L, communication_mode = communication_mode, A_val = 1L,
      nlambda_init = 20L, n_cores = 1L, precomputed_folds = folds,
      use_lambda_cache = use_cache, family = "binomial", verbose = FALSE
    ))
    reference <- fit(original, original_folds)
    perturbed <- fit(changed, changed_folds)
    before <- reference$fold_results[[1L]]$source_results$s1$per_k2_alpha$k2_3
    after <- perturbed$fold_results[[1L]]$source_results$s1$per_k2_alpha$k2_3
    # Both coefficient fits use fold 2 only. Only evaluation fold 3 changed.
    training_rows <- nuisance_folds[[2L]]$original_idx
    stopifnot(identical(original[[nuisance_site]]$Y[training_rows],
                        changed[[nuisance_site]]$Y[training_rows]))
    difference <- max(abs(as.numeric(before) - as.numeric(after)))
    if (!use_cache) stopifnot(difference < 1e-12)
    data.frame(communication_mode, nuisance_site, use_lambda_cache = use_cache,
               lambda_before = as.numeric(attr(before, "lambda_used")),
               lambda_after = as.numeric(attr(after, "lambda_used")),
               max_coefficient_change = difference,
               outer_fold = 1L, inner_evaluation_fold = 3L,
               coefficient_training_fold = 2L)
  })
  do.call(rbind, comparisons)
})
rows <- do.call(rbind, rows)
write.csv(rows, args[[1L]], row.names = FALSE)
print(rows)
