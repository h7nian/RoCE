# ============================================================================
# Tests for the multiplier (wild) bootstrap variance of the comparison baselines.
#
# Covers:
#   1. .multiplier_bootstrap_se   - reproduces the fixed-effects IF variance,
#                                   is RNG-neutral, and handles empty blocks.
#   2. .resolve_comparison_variance - bootstrap-vs-analytic selection + fallback.
#   3. run_all_comparisons         - all five comparison baselines report a
#                                   bootstrap SE by default (target_only stays
#                                   analytic), analytic mode is recoverable, and
#                                   the bootstrap SE matches the fixed-effects
#                                   variance (poolers) / analytic IF variance (DR).
# ============================================================================

# ----------------------------------------------------------------------------
# 1. .multiplier_bootstrap_se
# ----------------------------------------------------------------------------

test_that(".multiplier_bootstrap_se reproduces the fixed-effects IF variance", {
  set.seed(1L)
  infl_a <- rnorm(600L)
  infl_b <- rnorm(400L, sd = 2)
  blocks <- list(
    list(influence = infl_a, weight = 0.6),
    list(influence = infl_b, weight = 0.4)
  )

  center_var <- function(v) mean((v - mean(v))^2)
  analytic_se <- sqrt(
    0.6^2 * center_var(infl_a) / length(infl_a) +
    0.4^2 * center_var(infl_b) / length(infl_b)
  )

  set.seed(2L)
  boot_se <- .multiplier_bootstrap_se(blocks, n_bootstrap = 5000L)

  expect_true(is.finite(boot_se))
  # Wild bootstrap is unbiased for the fixed-effects variance; allow Monte-Carlo slack.
  expect_equal(boot_se, analytic_se, tolerance = 0.08)
})

test_that(".multiplier_bootstrap_se is RNG-neutral (restores the global seed)", {
  blocks <- list(list(influence = rnorm(100L), weight = 1))

  set.seed(99L)
  expected_next <- runif(1L)

  set.seed(99L)
  invisible(.multiplier_bootstrap_se(blocks, n_bootstrap = 200L))
  actual_next <- runif(1L)

  expect_equal(actual_next, expected_next)
})

test_that(".multiplier_bootstrap_se returns NA when there is nothing to resample", {
  expect_true(is.na(.multiplier_bootstrap_se(list())))
  expect_true(is.na(.multiplier_bootstrap_se(
    list(list(influence = numeric(0), weight = 1))
  )))
})

# ----------------------------------------------------------------------------
# 2. .resolve_comparison_variance
# ----------------------------------------------------------------------------

test_that(".resolve_comparison_variance selects bootstrap by default, analytic on request", {
  set.seed(7L)
  blocks <- list(list(influence = rnorm(300L), weight = 1))

  res_boot <- .resolve_comparison_variance(
    analytic_variance = 5, blocks = blocks, variance_method = "bootstrap"
  )
  expect_identical(res_boot$variance_method, "bootstrap")
  expect_true(is.finite(res_boot$se_bootstrap))
  expect_equal(res_boot$variance, res_boot$se_bootstrap^2, tolerance = 1e-10)
  expect_equal(res_boot$se, sqrt(res_boot$variance), tolerance = 1e-10)
  expect_equal(res_boot$variance_analytic, 5)

  res_analytic <- .resolve_comparison_variance(
    analytic_variance = 5, blocks = blocks, variance_method = "analytic"
  )
  expect_identical(res_analytic$variance_method, "analytic")
  expect_equal(res_analytic$variance, 5)
  expect_true(is.na(res_analytic$se_bootstrap))
})

test_that(".resolve_comparison_variance falls back to analytic when the bootstrap is undefined", {
  res <- .resolve_comparison_variance(
    analytic_variance = 3, blocks = list(), variance_method = "bootstrap"
  )
  expect_identical(res$variance_method, "analytic")
  expect_equal(res$variance, 3)
})

# ----------------------------------------------------------------------------
# 3. run_all_comparisons integration
# ----------------------------------------------------------------------------

.bootstrap_test_cache <- new.env(parent = emptyenv())

.get_bootstrap_test_data <- function() {
  if (!exists("data_split", envir = .bootstrap_test_cache, inherits = FALSE)) {
    .bootstrap_test_cache$data_split <- make_small_data_split(
      n_per_site = 80L, K = 2L, p = 4L, seed = 2026L,
      n_folds = 3L, outcome_type = "binary"
    )
  }
  .bootstrap_test_cache$data_split
}

.comparison_baselines <- c(
  "sample_size", "inverse_variance", "federated_dr", "pooled_dr", "tilted_aipw"
)

test_that("run_all_comparisons: all baselines report a bootstrap SE by default", {
  data_split <- .get_bootstrap_test_data()

  set.seed(11L)
  res <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L, A_val = 1L
  )

  for (method in .comparison_baselines) {
    r <- res[[method]]
    expect_identical(r$components$variance_method, "bootstrap",
                     info = sprintf("[%s] variance_method should default to bootstrap", method))
    expect_true(is.finite(r$se) && r$se > 0,
                info = sprintf("[%s] bootstrap se must be finite and positive", method))
    expect_true(is.finite(r$components$se_bootstrap),
                info = sprintf("[%s] components$se_bootstrap must be finite", method))
    expect_true(is.finite(r$components$variance_analytic),
                info = sprintf("[%s] analytic variance must be retained", method))
    expect_equal(r$se, sqrt(r$variance), tolerance = 1e-10,
                 info = sprintf("[%s] se must equal sqrt(variance)", method))
  }

  # target_only is the benchmark and keeps its analytic influence-function variance.
  expect_true(is.finite(res$target_only$se) && res$target_only$se > 0)
  expect_null(res$target_only$components$variance_method)
})

test_that("run_all_comparisons: bootstrap SE matches the fixed-effects / IF variance", {
  data_split <- .get_bootstrap_test_data()

  set.seed(11L)
  res <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L, A_val = 1L
  )

  # Poolers: the wild bootstrap reproduces the fixed-effects variance (drops tau^2).
  for (method in c("sample_size", "inverse_variance")) {
    r <- res[[method]]
    expect_equal(unname(r$components$se_bootstrap),
                 unname(sqrt(r$components$var_fixed_effects)),
                 tolerance = 0.25,
                 info = sprintf("[%s] bootstrap SE should match fixed-effects SE", method))
  }

  # DR methods: the bootstrap reproduces their analytic influence-function variance.
  for (method in c("federated_dr", "pooled_dr", "tilted_aipw")) {
    r <- res[[method]]
    expect_equal(r$components$se_bootstrap, r$components$se_analytic,
                 tolerance = 0.25,
                 info = sprintf("[%s] bootstrap SE should match analytic IF SE", method))
  }
})

test_that("run_all_comparisons: variance_method = 'analytic' recovers analytic SEs", {
  data_split <- .get_bootstrap_test_data()

  set.seed(11L)
  res <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L, A_val = 1L,
    variance_method = "analytic"
  )

  for (method in .comparison_baselines) {
    r <- res[[method]]
    expect_identical(r$components$variance_method, "analytic",
                     info = sprintf("[%s] variance_method should be analytic", method))
    expect_equal(r$se, r$components$se_analytic, tolerance = 1e-8,
                 info = sprintf("[%s] analytic se must equal se_analytic", method))
  }
})

test_that("run_all_comparisons: bootstrap SEs are reproducible given a seed", {
  data_split <- .get_bootstrap_test_data()

  set.seed(123L)
  res_a <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L, A_val = 1L
  )
  set.seed(123L)
  res_b <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L, A_val = 1L
  )

  for (method in .comparison_baselines) {
    expect_equal(res_a[[method]]$components$se_bootstrap,
                 res_b[[method]]$components$se_bootstrap,
                 tolerance = 1e-8,
                 info = sprintf("[%s] bootstrap SE must be reproducible", method))
  }
})
