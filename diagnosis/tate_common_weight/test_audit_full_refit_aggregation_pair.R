#!/usr/bin/env Rscript

library(testthat)
source("diagnosis/tate_common_weight/audit_full_refit_aggregation_intermediates.R")

make_pair <- function(rho = 1, sim = 7, positive_draw = 1,
                      source = "/tmp/source", positive_identity = FALSE) {
  rows <- list(
    data.frame(sim_id = sim, rho = rho, draw_id = 0, identity = TRUE,
               status = "completed"),
    data.frame(sim_id = sim, rho = rho, draw_id = positive_draw,
               identity = positive_identity, status = "completed")
  )
  saved <- list(
    list(source_bundle = source, summary = rows[[1L]]),
    list(source_bundle = source, summary = rows[[2L]])
  )
  list(rows = rows, saved = saved)
}

test_that("pair validator accepts only matching identity-positive draws", {
  pair <- make_pair()
  observed <- .afri_validate_draw_pair(pair$rows, pair$saved)
  expect_identical(observed, list(sim_id = 7, rho = 1,
                                  source_bundle = "/tmp/source"))

  mutations <- list(
    function(x) { x$rows[[2]]$sim_id <- 8; x },
    function(x) { x$rows[[2]]$rho <- 0; x },
    function(x) { x$saved[[2]]$source_bundle <- "/tmp/other"; x },
    function(x) { x$rows[[1]]$identity <- FALSE; x },
    function(x) { x$rows[[1]]$draw_id <- 1; x },
    function(x) { x$rows[[2]]$identity <- TRUE; x },
    function(x) { x$rows[[2]]$draw_id <- 0; x },
    function(x) { x$rows[[2]]$status <- "failed"; x },
    function(x) {
      x$rows[[1]]$sim_id <- x$saved[[1]]$summary$sim_id <- 1.5
      x
    },
    function(x) {
      x$rows[[2]]$draw_id <- x$saved[[2]]$summary$draw_id <- 1.5
      x
    }
  )
  for (mutate in mutations) {
    bad <- mutate(make_pair())
    expect_error(.afri_validate_draw_pair(bad$rows, bad$saved))
  }
  bad <- make_pair()
  bad$rows[[1]] <- rbind(bad$rows[[1]], bad$rows[[1]])
  expect_error(.afri_validate_draw_pair(bad$rows, bad$saved), "one-row")
  bad <- make_pair()
  bad$saved[[2L]]$summary$rho <- 0
  expect_error(.afri_validate_draw_pair(bad$rows, bad$saved), "CSV and RDS")
  expect_error(.afri_scalar(matrix(1, 1, 1), "scalar"), "scalar")
  expect_error(.afri_scalar(1 + 1i, "scalar"), "numeric scalar")
})

test_that("source fit selection uses the actual numeric rho key", {
  fits <- list(zero = list(marker = "rho0"), one = list(marker = "rho1"))
  source <- list(group_result = list(artifacts = list(
    "0" = list(direct_tate_results = list(one_round_crossfit = fits$zero)),
    "1" = list(direct_tate_results = list(one_round_crossfit = fits$one))
  )))
  expect_identical(.afri_select_source_fit(source, 0), fits$zero)
  expect_identical(.afri_select_source_fit(source, 1), fits$one)
  expect_error(.afri_select_source_fit(source, 0.5), "missing rho artifact")
})
