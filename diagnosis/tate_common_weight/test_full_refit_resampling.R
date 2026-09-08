#!/usr/bin/env Rscript
# Run from repository root with R_LIBS pointing to the frozen RoCE installation.
library(testthat)
library(RoCE)
source("diagnosis/tate_common_weight/full_refit_resampling.R")

make_refit_fixture <- function() {
  data <- lapply(seq_len(3L), function(s) {
    x <- cbind(1, seq_len(18L) + 100 * s)
    list(n = 18L, X = x, X_dagger = x, A = rep(0:1, 9L),
         Y = as.numeric(seq_len(18L)), Z_site_true = x,
         W_outcome_true = x, Z_site = x, W_outcome = x)
  })
  names(data) <- c("t", "s1", "s2")
  folds <- setNames(rep(list(split(seq_len(18L), rep(1:3, each = 6L))), 3L),
                    names(data))
  list(data = data, folds = folds)
}

test_that("identity preserves rows, fold positions, and caller RNG", {
  x <- make_refit_fixture()
  set.seed(71)
  old_rng <- .Random.seed
  b <- roce_refit_resample(x$data, x$folds, 811L, identity = TRUE)
  expect_identical(.Random.seed, old_rng)
  expect_identical(b$data_split, x$data)
  for (site in names(x$data)) {
    expect_identical(b$original_row_ids[[site]], seq_len(18L))
    expect_identical(b$multiplicities[[site]], rep(1L, 18L))
  }
  expect_identical(attr(b$precomputed_folds$target_folds, ".data_ref"), x$data$t)
})

test_that("row maps are shared across fields and remain inside original folds", {
  x <- make_refit_fixture()
  set.seed(19)
  old_rng <- .Random.seed
  b <- roce_refit_resample(x$data, x$folds, 901L)
  expect_identical(.Random.seed, old_rng)
  expect_identical(b, roce_refit_resample(x$data, x$folds, 901L))
  expect_false(identical(b$original_row_ids$t, seq_len(18L)))
  for (site in names(x$data)) {
    rows <- b$original_row_ids[[site]]
    expect_identical(sum(b$multiplicities[[site]]), 18L)
    for (idx in x$folds[[site]]) expect_true(all(rows[idx] %in% idx))
    for (field in setdiff(names(x$data[[site]]), "n")) {
      original <- x$data[[site]][[field]]
      expected <- if (is.matrix(original)) original[rows, , drop = FALSE] else original[rows]
      expect_identical(b$data_split[[site]][[field]], expected)
    }
    folds <- if (site == "t") b$precomputed_folds$target_folds else
      b$precomputed_folds$source_folds[[site]]
    expect_identical(attr(folds, ".data_ref"), b$data_split[[site]])
    for (k in seq_along(folds)) {
      view <- RoCE:::materialize_fold(folds, k)
      expect_identical(view$Y, x$data[[site]]$Y[rows[x$folds[[site]][[k]]]])
      expect_identical(view$original_idx, x$folds[[site]][[k]])
    }
  }
  order <- c("t", "s2", "s1")
  permuted <- roce_refit_resample(x$data[order], x$folds[order], 901L)
  expect_identical(permuted$original_row_ids[names(x$data)], b$original_row_ids)
})

test_that("invalid schemas, partitions and seeds fail closed", {
  x <- make_refit_fixture()
  for (seed in list(0, -1, 1.5, NA_real_, Inf, "1", c(1, 2))) {
    expect_error(roce_refit_resample(x$data, x$folds, seed), "seed")
  }
  expect_error(roce_refit_resample(x$data, x$folds, 1, identity = 1), "identity")
  bad <- x$data
  bad$t$unknown_metadata <- 1
  expect_error(roce_refit_resample(bad, x$folds, 1), "schema")
  bad <- x$data
  bad$t$X[1, 1] <- NA_real_
  expect_error(roce_refit_resample(bad, x$folds, 1), "matrix")
  bad <- x$folds
  bad$t[[1]][1] <- bad$t[[2]][1]
  expect_error(roce_refit_resample(x$data, bad, 1), "partition")
  bad <- x$folds
  bad$t[[1]][1] <- 1.5
  expect_error(roce_refit_resample(x$data, bad, 1), "partition")
  expect_error(roce_refit_resample(x$data, x$folds[c(2, 1, 3)], 1), "names")
})

# Optional real saved-fit contract check; not a nuisance-refit/coverage test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: test_full_refit_resampling.R [artifacts.rds]")
if (length(args) == 1L) {
  artifact <- readRDS(args[[1L]])$group_result$artifacts
  test_that("saved production estimates use the same fixed-weight functional", {
    for (entry in artifact) {
      fit <- entry$direct_tate_results$one_round_crossfit
      expect_equal(roce_refit_fixed_weight_estimate(fit, fit$fold_weights),
                   fit$estimate, tolerance = 1e-12)
      zero <- fit$fold_weights * 0
      expect_equal(roce_refit_fixed_weight_estimate(fit, zero),
                   fit$target_only$estimate, tolerance = 1e-12)
      sizes <- fit$intermediates$sample_sizes
      exact <- sum(vapply(seq_along(fit$intermediates$fold_info), function(k) {
        info <- fit$intermediates$fold_info[[k]]
        eta <- fit$fold_weights[k, ]
        target_fraction <- length(info$target_idx) / sizes$n_t
        source_fractions <- lengths(info$source_idx) / sizes$n_source
        target_fraction * ((1 - sum(eta)) * info$fold_target_estimate +
          sum(eta * info$mu_pred_ts)) +
          sum(source_fractions * eta * info$delta_ts)
      }, numeric(1L)))
      expect_equal(exact, fit$estimate, tolerance = 1e-12)
      # Mean fold weights are diagnostics, not weights for the pooled estimates.
      # The pooled identity is valid if that same vector is applied in ALL folds.
      constant_weights <- matrix(rep(fit$weights, each = fit$n_folds),
                                 nrow = fit$n_folds,
                                 dimnames = list(NULL, names(fit$weights)))
      pooled <- fit$target_only$estimate + sum(fit$weights *
        (fit$source_estimates - fit$target_only$estimate))
      expect_equal(roce_refit_fixed_weight_estimate(fit, constant_weights),
                   pooled, tolerance = 1e-12)
      counts <- lapply(entry$data_split, function(site) rep(1L, site$n))
      matched <- roce_refit_matched_fixed_nuisance(fit, counts)
      expect_equal(matched$estimate_original_weights, fit$estimate, tolerance = 1e-12)
      expect_equal(matched$estimate_relearned_weights, fit$estimate, tolerance = 1e-12)
      expect_equal(unname(matched$fold_weights), unname(fit$fold_weights), tolerance = 1e-12)
      # Nonconstant counts with exactly preserved fold-cell totals.
      for (info in fit$intermediates$fold_info) {
        idx <- info$target_idx[1:2]
        counts$t[idx] <- c(0L, 2L)
      }
      changed <- roce_refit_matched_fixed_nuisance(fit, counts)
      expect_true(is.finite(changed$estimate_relearned_weights))
      expected <- mean(fit$all_phi_tau * unlist(counts, use.names = FALSE))
      expect_equal(changed$estimate_original_weights, expected, tolerance = 1e-12)
      counts$t[1] <- counts$t[1] + 1L
      expect_error(roce_refit_matched_fixed_nuisance(fit, counts), "invalid counts")
    }
  })
}
