#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
source(arguments[[1L]])
library(testthat)

test_that("the original experiment preserves its frozen settings and seeds", {
  original <- make_gaussian_validation_settings("original")
  frozen <- read.csv(arguments[[2L]], stringsAsFactors = FALSE)
  expect_equal(original[names(frozen)], frozen)
  expect_equal(nrow(original), 145)
  expect_true(all(original$votes_required <= original$valid_minimum))
})

test_that("validity-fraction experiments keep assumptions and paired seeds explicit", {
  settings <- make_gaussian_validation_settings("validity_fraction")
  expect_equal(nrow(settings), 232)
  expect_setequal(settings$unshifted_fraction, c(.25, .5, .75))
  expect_true(all(settings$votes_required >= 1 & settings$votes_required <= settings$valid_minimum))
  expect_true(all(settings$assumption_valid))
  data_key <- paste(settings$num_sources, settings$shared_correlation, settings$bias_scale,
    ifelse(settings$bias_scale == 0, settings$num_sources, settings$valid_minimum), sep = ":")
  expect_true(all(vapply(split(settings$seed, data_key), function(values) length(unique(values)) == 1L, logical(1))))
  expect_equal(length(unique(settings$seed)), length(unique(data_key)))
  expect_error(make_gaussian_validation_settings("unrecognized"), "arg")
})

cat("GAUSSIAN_VALIDATION_SETTINGS_PASSED\n")
