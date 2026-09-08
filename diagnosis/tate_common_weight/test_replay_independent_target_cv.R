#!/usr/bin/env Rscript

library(testthat)
source("diagnosis/tate_common_weight/replay_independent_target_cv.R")

test_that("constant outcome fallback is not mislabeled as a CV replay", {
  expect_error(.ritcv_require_standard_cv(c(rep(0, 7), rep(1, 30)),
                                         "OR", "binomial"), "constant fallback")
  expect_true(.ritcv_require_standard_cv(rep(0:1, each = 8), "OR", "binomial"))
  expect_true(.ritcv_require_standard_cv(c(rep(0, 7), rep(1, 30)),
                                        "PS", "binomial"))
})

test_that("CV argument builder distinguishes internal and explicit folds", {
  x <- matrix(1:12, 6, 2); y <- 1:6
  internal <- .ritcv_build_glmnet_args(x, y, "gaussian", 5L, 100L, 99L)
  expect_identical(internal$nfolds, 5L)
  expect_null(internal$foldid)
  expect_true(internal$keep)
  explicit <- .ritcv_build_glmnet_args(
    x, y, "gaussian", 5L, 100L, 99L, rep(1:3, each = 2L)
  )
  expect_null(explicit$nfolds)
  expect_identical(explicit$foldid, rep(1:3, each = 2L))

  fake_maker <- function(groups, n_folds, caller) {
    if (is.null(groups) || !anyDuplicated(groups)) NULL else
      rep(seq_len(n_folds), length.out = length(groups))
  }
  unique_groups <- .ritcv_prepare_glmnet_args(
    x, y, "gaussian", 3L, 100L, 99L, seq_len(6L), "test",
    fake_maker
  )
  expect_identical(unique_groups$args$nfolds, 3L)
  expect_null(unique_groups$args$foldid)
  grouped <- .ritcv_prepare_glmnet_args(
    x, y, "gaussian", 3L, 100L, 99L, rep(1:3, each = 2L), "test",
    fake_maker
  )
  expect_null(grouped$args$nfolds)
  expect_identical(grouped$args$foldid, rep(1:3, length.out = 6L))
})

test_that("prediction identity validator fails closed", {
  expect_equal(.ritcv_validate_prediction(c(1, 2), c(1, 2)), 0)
  expect_lt(abs(
    .ritcv_validate_prediction(c(1, 2 + 1e-11), c(1, 2)) - 1e-11
  ), 1e-15)
  expect_error(.ritcv_validate_prediction(c(1, 2 + 1e-5), c(1, 2)), "tolerance")
  expect_error(.ritcv_validate_prediction(c(1, NA), c(1, 2)), "finite")
  expect_error(.ritcv_validate_prediction(1, c(1, 2)), "equal-length")
  expect_error(.ritcv_validate_prediction(character(), numeric()), "numeric")
  expect_error(.ritcv_validate_prediction(matrix(1:2), 1:2), "finite numeric")
  expect_error(.ritcv_validate_prediction(structure(1:2, class = "x"), 1:2),
               "finite numeric")
  for (tolerance in list(0, -1, NA_real_, Inf, c(1e-10, 1e-9),
                         matrix(1e-10))) {
    expect_error(.ritcv_validate_prediction(1:2, 1:2, tolerance), "tolerance")
  }
})

test_that("integer identities reject truncation and out-of-range values", {
  expect_identical(.ritcv_integer_scalar(13, "task", 1L, 100L), 13L)
  for (value in list(13.5, 0, 101, NA_real_, Inf, "13", matrix(13))) {
    expect_error(.ritcv_integer_scalar(value, "task", 1L, 100L), "integer")
  }
})

test_that("checksum gate rejects extra payloads and directories", {
  fixture <- tempfile("ritcv-sha-"); dir.create(fixture)
  on.exit(unlink(fixture, recursive = TRUE), add = TRUE)
  fake_sha <- function(path) strrep("a", 64L)
  writeLines("payload", file.path(fixture, "a.csv"))
  writeLines(paste(fake_sha("a"), "a.csv", sep = "  "),
             file.path(fixture, "sha256.txt"))
  expect_identical(.ritcv_check_sha(fixture, "a.csv", fake_sha), strrep("a", 64L))
  writeLines("extra", file.path(fixture, "extra.txt"))
  expect_error(.ritcv_check_sha(fixture, "a.csv", fake_sha), "validation")
  unlink(file.path(fixture, "extra.txt"))
  dir.create(file.path(fixture, "extra_dir"))
  expect_error(.ritcv_check_sha(fixture, "a.csv", fake_sha), "validation")
})

test_that("failure publication retains completed context without passed marker", {
  output <- tempfile("ritcv-failure-")
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  state <- new.env(parent = emptyenv())
  state$bundle <- "/bundle"; state$bundle_hash <- "bundle-hash"
  state$audit_hash <- "audit-hash"
  state$task_id <- 13L; state$sim_id <- 10013L
  state$current <- list(model = "PS", arm = 1L, k1 = 5L, k2 = NA_integer_)
  state$diagnostics <- list(data.frame(model = "OR", prediction_max_error = 0))
  state$warning_events <- list(data.frame(
    event_index = 1L, model = "OR", arm = 1L, k1 = 1L, k2 = NA_integer_,
    message = "warning", condition_class = "simpleWarning;warning;condition",
    condition_call = "fit()"
  ))
  fake_atomic <- function(path, writer, caller) { dir.create(path); writer(path) }
  fake_sha <- function(path) strrep("b", 64L)
  .ritcv_publish_failure(output, simpleError("boom"), state, fake_sha, fake_atomic)
  attempt <- read.csv(file.path(output, "attempt_failure.csv"))
  expect_identical(attempt$completed_call_count, 1L)
  expect_identical(attempt$current_model, "PS")
  expect_equal(nrow(read.csv(file.path(output, "completed_calls.csv"))), 1L)
  expect_equal(nrow(read.csv(file.path(output, "warning_events.csv"))), 1L)
  expect_false(any(grepl("passed", readLines(file.path(output, "metadata.txt")))))
})

test_that("output validation rejects malformed and existing paths", {
  expect_error(.ritcv_output_path(NA_character_), "concrete")
  expect_error(.ritcv_output_path("."), "concrete")
  existing <- tempfile("replay-existing-"); dir.create(existing)
  on.exit(unlink(existing, recursive = TRUE), add = TRUE)
  expect_error(.ritcv_output_path(existing), "already exists")
  candidate <- tempfile("new-replay-output-")
  expect_identical(.ritcv_output_path(candidate), candidate)
})
