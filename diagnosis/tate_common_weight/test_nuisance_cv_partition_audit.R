#!/usr/bin/env Rscript
# Run from repository root with R_LIBS pointing at the tested v19 library.

library(testthat)
library(RoCE)
source("diagnosis/tate_common_weight/nuisance_cv_partition_audit.R")
source("diagnosis/tate_common_weight/full_refit_resampling.R")

partition_record <- function(groups, folds, n_folds = 3L, valid = NULL) {
  record <- list(
    caller = "unit-test", n_folds = n_folds,
    cv_group_id = groups, cv_fold_id = folds
  )
  if (!is.null(valid)) record$valid <- valid
  record
}

test_that("partition record validation fails closed", {
  expect_true(roce_cv_partition_record_valid(
    partition_record(c(1L, 1L, 2L, 2L, 3L, 3L), c(1L, 1L, 2L, 2L, 3L, 3L))
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(c(1L, 1L, 2L, 2L), c(1L, 2L, 2L, 3L))
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(1:6, c(1L, 2L, NA, 1L, 2L, 3L))
  ))
  expect_true(roce_cv_partition_record_valid(partition_record(1:6, NULL)))
  expect_false(roce_cv_partition_record_valid(
    partition_record(c(1L, 1L, 2:5), NULL)
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(1:6, c(1L, 1L, 2L, 2L, 2L, 2L))
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(c(1, 1, 2, 2.5, 3, 3), c(1, 1, 2, 2, 3, 3))
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(matrix(1:6, ncol = 1L), rep(1:3, each = 2L))
  ))
  expect_false(roce_cv_partition_record_valid(
    partition_record(c(1L, 1L, 2L, 2L), c(1L, 2L, 2L, 3L), valid = TRUE)
  ))
})

test_that("runtime recorder preserves semantics, RNG, and namespace state", {
  namespace <- asNamespace("RoCE")
  symbol <- ".make_nuisance_cv_fold_id"
  original <- get(symbol, envir = namespace, inherits = FALSE)
  fit <- function() RoCE:::.make_nuisance_cv_fold_id(
    rep(1:6, each = 2L), 3L, "runtime-test"
  )

  set.seed(8101)
  direct <- fit()
  rng_direct <- .Random.seed
  set.seed(8101)
  audited <- roce_with_nuisance_cv_audit(fit)
  expect_identical(audited$fitted, direct)
  expect_identical(.Random.seed, rng_direct)
  expect_length(audited$partitions, 1L)
  expect_true(roce_cv_partition_record_valid(audited$partitions[[1L]]))
  expect_false(inherits(get(symbol, envir = namespace), "functionWithTrace"))
  expect_identical(body(get(symbol, envir = namespace)), body(original))

  failed <- roce_with_nuisance_cv_audit(function() {
    fit()
    stop("intentional fit failure")
  })
  expect_s3_class(failed$fitted, "error")
  expect_match(conditionMessage(failed$fitted), "intentional fit failure")
  expect_length(failed$partitions, 1L)
  expect_false(inherits(get(symbol, envir = namespace), "functionWithTrace"))

  trace(symbol, quote(invisible(NULL)), print = FALSE, where = namespace)
  on.exit(untrace(symbol, where = namespace), add = TRUE)
  expect_error(roce_with_nuisance_cv_audit(function() NULL), "existing nuisance-CV trace")
})

make_small_refit_fixture <- function() {
  set.seed(8102)
  sites <- lapply(seq_len(3L), function(site) {
    n <- 300L
    x <- matrix(rnorm(n * 2L), ncol = 2L)
    a <- rep(0:1, length.out = n)
    y <- 0.5 * a + 0.35 * x[, 1L] - 0.2 * x[, 2L] +
      0.1 * site + rnorm(n)
    list(
      n = n, X = x, X_dagger = x, Z_site_true = x,
      W_outcome_true = x, Z_site = x, W_outcome = x, A = a, Y = y
    )
  })
  names(sites) <- c("t", "s1", "s2")
  views <- lapply(sites, function(site) {
    fold_id <- (seq_len(site$n) - 1L) %% 3L + 1L
    folds <- lapply(seq_len(3L), function(k) {
      idx <- which(fold_id == k)
      list(original_idx = idx, n = length(idx))
    })
    attr(folds, ".data_ref") <- site
    folds
  })
  list(
    data = sites,
    folds = list(target_folds = views$t, source_folds = views[c("s1", "s2")])
  )
}

test_that("small origin-mode full refits record every nuisance stage", {
  fixture <- make_small_refit_fixture()
  set.seed(8103)
  reference <- RoCE::run_tate_crossfit(
    data_split = fixture$data, n_folds = 3L,
    communication_mode = "one_round", family = "gaussian",
    nlambda_init = 8L, n_cores = 1L, parallel_arms = FALSE,
    verbose = FALSE, precomputed_folds = fixture$folds
  )

  run_refit <- function(seed, identity) {
    set.seed(8103)
    roce_refit_once(
      fixture$data, reference, seed = seed, identity = identity,
      n_cores = 1L, parallel_arms = FALSE,
      nuisance_cv_grouping = "origin"
    )
  }
  identity <- run_refit(8104L, TRUE)
  resampled <- run_refit(8105L, FALSE)
  for (result in list(identity, resampled)) {
    expect_null(result$failure_message)
    expect_true(length(result$nuisance_cv_partitions) > 0L)
    expect_true(all(vapply(
      result$nuisance_cv_partitions,
      roce_cv_partition_record_valid, logical(1L)
    )))
    expect_setequal(
      unique(vapply(result$nuisance_cv_partitions, `[[`, character(1L), "caller")),
      c(
        "estimate_complement_fold_aipw PS",
        "estimate_complement_fold_aipw OR",
        "fit_initial_outcome", "fit_initial_density_ratio",
        "fit_unified_density_ratio", "fit_unified_outcome"
      )
    )
  }
  expect_true(identity$identity)
  expect_false(resampled$identity)
  expect_equal(identity$fitted$estimate, reference$estimate, tolerance = 1e-10)
  expect_equal(identity$fitted$se, reference$se, tolerance = 1e-10)
})
