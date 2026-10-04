.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
suppressPackageStartupMessages(library(testthat))
source("diagnosis/target_rcal/target_rcal_helpers.R")

# Low-dimensional fixtures test the production FACE branch and both protocols
# at 1000/site. They are not included in the p=100 performance experiment.
for (configuration in c("C1", "C2", "C3")) {
  test_that(paste("FACE", configuration, "preserves both protocol contracts"), {
    set.seed(510)
    data <- RoCE:::generate_simulation_data(
      n_total = 2000L, K = 1L, p = 4L, config = configuration,
      estimand_type = "superpopulation", outcome_type = "binary", dgp_type = "face",
      ate_deviation = 0, n_deviated_sites = 0L, warn_ignored = FALSE
    )
    split <- RoCE:::split_data_by_site(data)
    original_split <- serialize(split, NULL)
    folds <- RoCE:::build_crossfit_folds(split, 3L)
    for (mode in c("one_round", "two_round")) {
      fit <- RoCE::run_tate_crossfit(
        split, n_folds = 3L, communication_mode = mode,
        precomputed_folds = folds, nlambda_init = 5L, n_cores = 1L,
        family = "binomial", verbose = FALSE
      )
      expect_true(is.finite(fit$estimate))
      expect_gt(fit$se, 0)
      expect_equal(fit$estimate, mean(fit$all_phi_agg), tolerance = 1e-12)
      expect_equal(fit$se^2, fit$variance, tolerance = 1e-12)
      variants <- lapply(fit$arm_results, function(arm) lapply(arm$fold_results, function(fold) {
        list(outer = list(lasso = fold$target_only),
             inner = lapply(fold$target_only_inner, function(target) list(lasso = target)))
      }))
      reconstructed <- reaggregate_target_variant(split, folds, fit, variants, "lasso")
      expect_equal(reconstructed$estimate, fit$estimate, tolerance = 1e-12)
      expect_equal(reconstructed$se, fit$se, tolerance = 1e-10)
      expect_equal(reconstructed$fold_weights, fit$fold_weights, tolerance = 1e-10)
      expect_identical(serialize(split, NULL), original_split)
    }
  })
}
