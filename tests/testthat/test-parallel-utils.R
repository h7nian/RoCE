library(testthat)

test_that("parallel_lapply promotes Unix worker failures with item context", {
  skip_on_os("windows")

  expect_error(
    suppressWarnings(parallel_lapply(
      c(good = 1L, bad = 2L),
      function(value) {
        if (value == 2L) stop("intentional worker failure")
        value
      },
      n_cores = 2L
    )),
    "item\\(s\\) bad failed:.*intentional worker failure"
  )
})

test_that("parallel result validation preserves successful named results", {
  results <- list(first = 1L, second = 2L)

  expect_identical(
    RoCE:::.validate_parallel_results(results, results),
    results
  )
})
