# test-estimators-smoke.R - End-to-end happy-path smoke tests for exported
# estimators that are otherwise only covered by entry-validation tests.
#
# Covers (one happy-path test_that block each):
#   - estimate_tilted_aipw
#   - estimate_federated_dr
#   - estimate_pooled_dr
#   - estimate_sample_size_weighted
#   - estimate_inverse_variance_weighted
#   - estimate_oracle_dr
#   - run_all_comparisons (wrapper)
#
# Each test verifies that the estimator runs on a small synthetic data_split
# without error/warning and produces the documented output fields with finite,
# well-defined values (estimate, variance, se, ci_lower/ci_upper when present).
#
# Failures here mean a regression in the comparison-method code path that the
# entry-validation tests in test-utils.R cannot catch (those only check that
# *bad* inputs are rejected; they do not check that good inputs succeed).

library(testthat)

# Uses shared make_small_data_split() from helper-data.R. The K=2 split gives
# SSW / IVW more than one source site so the random-effects variance has a
# well-defined Cochran's-Q step. n_per_site is kept comfortably large so that the
# per-site treated-arm outcome CV folds are never degenerate (a single treated
# fold with one outcome class makes estimate_target_only_crossfit fail-fast); the
# comparison methods are happy-path estimators, not meant for pathologically tiny sites.

.estimator_smoke_cache <- new.env(parent = emptyenv())

.get_estimator_smoke <- function() {
  if (!exists("payload", envir = .estimator_smoke_cache, inherits = FALSE)) {
    data_split <- make_small_data_split(
      n_per_site = 200L, K = 2L, p = 4L,
      seed = 2026L, n_folds = 3L,
      outcome_type = "binary"
    )
    .estimator_smoke_cache$payload <- list(
      data_split = data_split,
      data       = attr(data_split, "raw_data")
    )
  }
  .estimator_smoke_cache$payload
}

.check_point_estimator <- function(res, method_label,
                                   require_ci = TRUE) {
  expect_true(is.list(res),
              info = sprintf("[%s] result must be a list", method_label))
  expect_true("estimate" %in% names(res),
              info = sprintf("[%s] result must include 'estimate'", method_label))
  expect_true("variance" %in% names(res),
              info = sprintf("[%s] result must include 'variance'", method_label))
  expect_true(is.finite(res$estimate),
              info = sprintf("[%s] estimate must be finite", method_label))
  expect_true(is.finite(res$variance) && res$variance >= 0,
              info = sprintf("[%s] variance must be finite and non-negative", method_label))

  if ("se" %in% names(res)) {
    expect_equal(res$se, sqrt(res$variance), tolerance = 1e-10,
                 info = sprintf("[%s] se must equal sqrt(variance)", method_label))
  }
  if (require_ci && all(c("ci_lower", "ci_upper") %in% names(res))) {
    expect_true(is.finite(res$ci_lower) && is.finite(res$ci_upper),
                info = sprintf("[%s] CI endpoints must be finite", method_label))
    expect_lte(res$ci_lower, res$estimate)
    expect_gte(res$ci_upper, res$estimate)
  }
}

# ============================================================================
# Comparison method smokes
# ============================================================================

test_that("estimate_tilted_aipw: runs on valid data and returns finite output", {
  payload <- .get_estimator_smoke()
  res <- estimate_tilted_aipw(payload$data_split,
                              family = "binomial", A_val = 1L)
  .check_point_estimator(res, "tilted_aipw", require_ci = FALSE)
})

test_that("estimate_federated_dr: runs on valid data and returns finite output", {
  payload <- .get_estimator_smoke()
  res <- estimate_federated_dr(payload$data_split,
                               dr_lambda = 0.05, A_val = 1L, family = "binomial")
  .check_point_estimator(res, "federated_dr")
})

test_that("estimate_pooled_dr: runs on valid data and returns finite output", {
  payload <- .get_estimator_smoke()
  res <- estimate_pooled_dr(payload$data_split,
                            dr_lambda = 0.05, A_val = 1L, family = "binomial")
  .check_point_estimator(res, "pooled_dr")
})

test_that("estimate_sample_size_weighted: runs on valid data and returns finite output", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  payload <- .get_estimator_smoke()
  res <- estimate_sample_size_weighted(payload$data_split,
                                       family = "binomial",
                                       use_crossfit = TRUE, n_folds = 3L,
                                       A_val = 1L)
  .check_point_estimator(res, "sample_size_weighted", require_ci = FALSE)
})

test_that("estimate_inverse_variance_weighted: runs on valid data and returns finite output", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  payload <- .get_estimator_smoke()
  res <- estimate_inverse_variance_weighted(payload$data_split,
                                            family = "binomial",
                                            use_crossfit = TRUE, n_folds = 3L,
                                            A_val = 1L)
  .check_point_estimator(res, "inverse_variance_weighted", require_ci = FALSE)
})

# ============================================================================
# Oracle estimator smoke
# ============================================================================

test_that("estimate_oracle_dr: runs with true nuisance parameters and is finite", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  payload <- .get_estimator_smoke()
  data <- payload$data
  skip_if(is.null(data$gamma_params) || is.null(data$alpha1_true),
          message = "RoCE truth (gamma_params + alpha1_true) not exposed.")

  res <- estimate_oracle_dr(
    payload$data_split,
    alpha1_true  = data$alpha1_true,
    gamma_params = data$gamma_params,
    outcome_type = "binary",
    A_val        = 1L,
    lambda_selection = 0.05
  )
  .check_point_estimator(res, "oracle_dr", require_ci = FALSE)
})

test_that("estimate_oracle_dr: default CV uses legacy min aggregation rule", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  payload <- .get_estimator_smoke()
  data <- payload$data
  skip_if(is.null(data$gamma_params) || is.null(data$alpha1_true),
          message = "RoCE truth (gamma_params + alpha1_true) not exposed.")

  res <- estimate_oracle_dr(
    payload$data_split,
    alpha1_true  = data$alpha1_true,
    gamma_params = data$gamma_params,
    outcome_type = "binary",
    A_val        = 1L,
    lambda_grid  = c(0.01, 0.05, 0.1)
  )
  .check_point_estimator(res, "oracle_dr", require_ci = FALSE)
  expect_equal(res$aggregation_lambda_rule, "min")
})

# ============================================================================
# run_all_comparisons wrapper smoke
# ============================================================================

test_that("run_all_comparisons returns a list with estimates for every selected method", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")
  payload <- .get_estimator_smoke()

  res <- run_all_comparisons(payload$data_split,
                             use_rcal = FALSE,
                             use_crossfit = TRUE,
                             n_folds = 3L,
                             family = "binomial",
                             A_val = 1L)

  required_methods <- c("sample_size", "inverse_variance",
                        "federated_dr", "pooled_dr", "tilted_aipw")
  missing_methods <- setdiff(required_methods, names(res))
  expect_equal(length(missing_methods), 0L,
               info = sprintf("run_all_comparisons missing methods: %s",
                              paste(missing_methods, collapse = ", ")))

  for (m in intersect(required_methods, names(res))) {
    .check_point_estimator(res[[m]], paste0("run_all_comparisons$", m),
                           require_ci = FALSE)
  }
})
