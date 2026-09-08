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

test_that("bootstrap replicate counts are validated consistently", {
  for (invalid in list(1L, 2.5, NA_real_, Inf, .Machine$integer.max + 1)) {
    expect_error(
      .validate_bootstrap_replicates(invalid, "test"),
      "integer >= 2",
      fixed = TRUE
    )
  }
  expect_identical(.validate_bootstrap_replicates(37, "test"), 37L)
})

# ----------------------------------------------------------------------------
# 2. .resolve_comparison_variance
# ----------------------------------------------------------------------------

test_that(".resolve_comparison_variance selects bootstrap by default, analytic on request", {
  set.seed(7L)
  blocks <- list(list(influence = rnorm(300L), weight = 1))

  res_boot <- .resolve_comparison_variance(
    analytic_variance = 5, blocks = blocks, variance_method = "bootstrap",
    n_bootstrap = 37L
  )
  expect_identical(res_boot$variance_method, "bootstrap")
  expect_true(is.finite(res_boot$se_bootstrap))
  expect_equal(res_boot$variance, res_boot$se_bootstrap^2, tolerance = 1e-10)
  expect_equal(res_boot$se, sqrt(res_boot$variance), tolerance = 1e-10)
  expect_equal(res_boot$variance_analytic, 5)
  expect_identical(res_boot$n_bootstrap, 37L)

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

test_that("run_all_comparisons propagates an explicit bootstrap count", {
  data_split <- .get_bootstrap_test_data()

  set.seed(12L)
  res <- run_all_comparisons(
    data_split,
    family = "binomial",
    use_crossfit = TRUE,
    n_folds = 3L,
    A_val = 1L,
    methods = c("sample_size", "inverse_variance"),
    n_bootstrap = 37L
  )

  expect_true(all(vapply(
    res,
    function(method) identical(method$components$n_bootstrap, 37L),
    logical(1L)
  )))
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

test_that("site-level SS/IVW fits are invariant to the requested core count", {
  skip_on_os("windows")
  data_split <- .get_bootstrap_test_data()
  methods <- c("sample_size", "inverse_variance")

  set.seed(812L)
  sequential <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L,
    A_val = 1L, variance_method = "analytic", methods = methods,
    n_cores = 1L
  )
  set.seed(812L)
  parallel <- run_all_comparisons(
    data_split, family = "binomial", use_crossfit = TRUE, n_folds = 3L,
    A_val = 1L, variance_method = "analytic", methods = methods,
    n_cores = 2L
  )

  for (method in methods) {
    expect_equal(parallel[[method]]$estimate, sequential[[method]]$estimate)
    expect_equal(parallel[[method]]$se, sequential[[method]]$se)
    expect_equal(
      parallel[[method]]$influence_blocks,
      sequential[[method]]$influence_blocks
    )
  }
})

test_that("DR density-ratio weights are reused across methods and treatment arms", {
  data_split <- .get_bootstrap_test_data()
  methods <- c("federated_dr", "pooled_dr")

  shared_weights <- .resolve_dr_weights_by_site(
    data_split,
    n_cores = 1L,
    caller = "test"
  )
  expect_setequal(names(shared_weights), setdiff(names(data_split), "t"))
  if (.Platform$OS.type != "windows") {
    parallel_weights <- .resolve_dr_weights_by_site(
      data_split,
      n_cores = 2L,
      caller = "test"
    )
    expect_equal(parallel_weights, shared_weights)
  }

  set.seed(778L)
  expected_next_draw <- stats::runif(1L)
  set.seed(778L)
  invisible(.resolve_dr_weights_by_site(
    data_split,
    n_cores = 1L,
    caller = "test"
  ))
  expect_equal(stats::runif(1L), expected_next_draw)

  set.seed(777L)
  uncached <- run_all_comparisons_tate(
    data_split,
    family = "binomial",
    variance_method = "analytic",
    methods = methods,
    n_cores = 1L
  )
  set.seed(777L)
  cached <- run_all_comparisons_tate(
    data_split,
    family = "binomial",
    variance_method = "analytic",
    methods = methods,
    n_cores = 1L,
    dr_weights_by_site = shared_weights
  )

  for (method in methods) {
    expect_equal(cached[[method]]$estimate, uncached[[method]]$estimate)
    expect_equal(cached[[method]]$se, uncached[[method]]$se)
    expect_equal(
      cached[[method]]$influence_blocks,
      uncached[[method]]$influence_blocks
    )
    component <- cached[[method]]$components
    expect_equal(
      component$se_analytic^2,
      component$mu1_influence_se^2 + component$mu0_influence_se^2 -
        2 * component$cross_arm_covariance,
      tolerance = 1e-12
    )
    expect_lte(abs(component$cross_arm_correlation), 1 + 1e-10)
  }
})

test_that("DR IF correction matches the full-source density-ratio score", {
  data_split <- .get_bootstrap_test_data()
  source_site <- setdiff(names(data_split), "t")[[1L]]
  shared_weights <- .resolve_dr_weights_by_site(
    data_split, n_cores = 1L, caller = "test"
  )
  components <- .resolve_dr_site_components(
    data_split = data_split,
    dr_weights_by_site = shared_weights,
    family = "binomial",
    A_val = 1L,
    required_sites = source_site,
    n_cores = 1L,
    caller = "test"
  )

  component <- components[[source_site]]
  source_data <- data_split[[source_site]]
  density_weights <- as.numeric(shared_weights[[source_site]])
  target_mean <- c(1, colMeans(as.matrix(data_split$t$Z_site)))
  centered_design <- sweep(
    cbind(1, as.matrix(source_data$Z_site)), 2L, target_mean, "-"
  )
  centered_pseudo_outcome <-
    component$base_result$phi - component$base_result$estimate
  derivative <- -colMeans(
    density_weights * centered_design * centered_pseudo_outcome
  ) / max(mean(density_weights), DIVISION_FLOOR)
  jacobian <- t(centered_design) %*%
    (centered_design * density_weights) / source_data$n
  adjustment <- solve_with_ridge(jacobian) %*% derivative
  expected_correction <- as.numeric(
    density_weights * (centered_design %*% adjustment)
  )
  expected_influence <-
    as.numeric(component$base_result$varphi_ot) + expected_correction
  expected_influence <- expected_influence - mean(expected_influence)

  expect_equal(component$source_if_correction, expected_correction)
  expect_equal(component$varphi_ot, expected_influence)
})

test_that("DR density-ratio cache validates source names and vector lengths", {
  data_split <- .get_bootstrap_test_data()
  source_sites <- setdiff(names(data_split), "t")
  malformed <- stats::setNames(
    lapply(source_sites, function(site) rep(1, data_split[[site]]$n)),
    source_sites
  )
  malformed[[source_sites[[1L]]]] <- 1

  expect_error(
    .resolve_dr_weights_by_site(
      data_split,
      dr_weights_by_site = malformed,
      caller = "test"
    ),
    "must contain"
  )
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
