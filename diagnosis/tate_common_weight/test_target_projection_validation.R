library(testthat)
source(if (file.exists("target_projection_validation.R")) "target_projection_validation.R" else
         "diagnosis/tate_common_weight/target_projection_validation.R")
source(if (file.exists("sparse_moment_projection.R")) "sparse_moment_projection.R" else
         "diagnosis/tate_common_weight/sparse_moment_projection.R")

loss_fixture <- function() {
  x <- expand.grid(candidate = c("null", "c0.5", "c1", "c2"), fold = 2:5,
                   stringsAsFactors = FALSE)
  x$n_validation <- c(199, 200, 200, 201)[x$fold-1L]
  x$loss <- 0; x$failed <- FALSE
  x
}

test_that("ties prefer null and then larger regularization constants", {
  x <- loss_fixture()
  expect_identical(.select_target_projection_candidate(x)$selected, "null")
  x$loss[x$candidate != "null"] <- -.1
  expect_identical(.select_target_projection_candidate(x)$selected, "c2")
  expect_identical(.select_target_projection_candidate(x[nrow(x):1, ])$selected, "c2")
})

test_that("all folds count and failed candidates are not partially averaged", {
  x <- loss_fixture()
  x$loss[x$candidate == "c0.5"] <- -2
  x$failed[x$candidate == "c0.5" & x$fold == 2] <- TRUE
  x$loss[x$candidate == "c1"] <- -.2
  result <- .select_target_projection_candidate(x)
  expect_identical(result$selected, "c1")
  expect_false(result$scores$eligible[result$scores$candidate == "c0.5"])
  expect_true(is.na(result$scores$validation_risk[result$scores$candidate == "c0.5"]))
  x$failed[x$candidate != "null"] <- TRUE
  expect_true(.select_target_projection_candidate(x)$all_nonnull_ineligible)
})

test_that("validation losses use actual fold sizes", {
  x <- loss_fixture()
  x$loss[x$candidate == "c1"] <- c(-1, 0, 0, 0)
  result <- .select_target_projection_candidate(x)
  expect_equal(result$scores$validation_risk[result$scores$candidate == "c1"], -199/800)
})

test_that("three-fold nested selection uses the same risk and failure policy", {
  x <- loss_fixture()
  x <- x[x$fold %in% 3:5, ]
  x$loss[x$candidate == "c1"] <- c(-1, 0, 0)
  result <- .select_target_projection_candidate(x, 3:5)
  expect_identical(result$selected, "c1")
  expect_equal(result$scores$validation_risk[result$scores$candidate == "c1"], -200/601)
  x$failed[x$candidate == "c1" & x$fold == 3] <- TRUE
  expect_identical(.select_target_projection_candidate(x, 3:5)$selected, "null")
  expect_error(.select_target_projection_candidate(x[-1, ], 3:5), "complete")
  expect_error(.select_target_projection_candidate(x, 4:5), "validation folds")
})

test_that("selection supports every outer complement without changing the policy", {
  for (outer_fold in 1:5) {
    x <- loss_fixture()
    validation_folds <- setdiff(1:5, outer_fold)
    x$fold <- validation_folds[x$fold-1L]
    x$loss[x$candidate == "c1"] <- -.2
    expect_identical(.select_target_projection_candidate(x, validation_folds)$selected, "c1")
    expect_error(.select_target_projection_candidate(x, c(outer_fold, outer_fold, 1, 2)),
                 "validation folds")
  }
})

test_that("incomplete, nonfinite and inconsistent inputs fail closed", {
  x <- loss_fixture()
  expect_error(.select_target_projection_candidate(x[-1, ]), "complete")
  x$loss[2] <- NA_real_
  expect_error(.select_target_projection_candidate(x), "nonfinite")
  x <- loss_fixture(); x$n_validation[2] <- 1
  expect_error(.select_target_projection_candidate(x), "fold sizes")
  x <- loss_fixture(); x$loss[x$candidate == "null"] <- -.01
  expect_error(.select_target_projection_candidate(x), "zero loss")
})

test_that("projection fitting validates its partition and noise-scale inputs", {
  system <- list(gradient = c(.3, -.4, .2), jacobian = -diag(3),
    index = list(propensity = 1L, outcome_treated = 2L, outcome_control = 3L))
  rows <- matrix(c(-.1, 0, .1, 0, 0, -.1, .1, 0, 0, -.2, .2, 0), 4, 3)
  fit <- .fit_target_projection_candidate(system, rows, penalty_scale = .5)
  expected_penalty <- .5*apply(rows, 2L, sd)*sqrt(log(6)/4)
  expect_equal(fit$penalties, expected_penalty)
  expect_equal(fit$coefficients, -sign(system$gradient)*
                 pmax(abs(system$gradient)-expected_penalty, 0))
  expect_error(.fit_target_projection_candidate(system, rows, penalty_scale = -1), "penalty_scale")
  expect_error(.fit_target_projection_candidate(system, as.numeric(rows), 1), "derivative rows")
  bad_rows <- rows; bad_rows[1, 1] <- NA_real_
  expect_error(.fit_target_projection_candidate(system, bad_rows, 1), "derivative rows")
  bad <- system; bad$index$outcome_control <- 2L
  expect_error(.fit_target_projection_candidate(bad, rows, 1), "block partition")
  expect_error(.fit_target_projection_candidate(1, rows, 1), "target nuisance system")
})
