test_that("standard source fitting solves the unweighted outcome objective", {
  set.seed(4601)
  site <- function(shift) {
    x <- matrix(rnorm(1800), 600, 3)
    x[, 1] <- x[, 1] + shift
    a <- rep(0:1, 300)
    list(n = 600L, W_outcome = x, Z_site = x, A = a,
         Y = .2 + .8 * a + .5 * x[, 1] - .2 * x[, 2] + rnorm(600))
  }
  data <- list(t = site(.2), s1 = site(-.1))
  folds <- build_crossfit_folds(data, 4L)
  message <- .prepare_source_training_message(folds$target_folds, 3:4)
  train <- combine_folds(folds$source_folds$s1, 3:4)
  for (arm in 0:1) {
    fit <- .fit_standard_source_nuisances(folds$source_folds$s1, "s1", 3:4, message,
      arm, 0L, 0L, 12, 6L, 1000L, "min", tol = 1e-10, weight_lambda = 0, outcome_lambda = 0)
    selected <- train$A == arm
    design <- cbind(1, train$W_outcome[selected, , drop = FALSE])
    expected <- solve(crossprod(design), crossprod(design, train$Y[selected]))
    expect_equal(as.numeric(fit$outcome), as.numeric(expected), tolerance = 1e-7)
    expect_identical(fit$training_folds, 3:4)
    expect_length(fit$calibration_folds, 0L)
    expect_identical(attr(fit$weight, "cv_seed"), .nuisance_training_seed(
      .nuisance_training_key("fit_unified_density_ratio", "s1", 3:4, arm)))
    expect_true(isTRUE(attr(fit$outcome, "converged")))
  }
})

test_that("complete standard training excludes outer and inner validation observations", {
  set.seed(4602)
  site <- function() {
    x <- matrix(rnorm(1800), 600, 3)
    list(n = 600L, W_outcome = x, Z_site = x, A = rep(0:1, 300),
         Y = .5 + x[, 1] + rnorm(600))
  }
  data <- list(t = site(), s1 = site())
  folds <- build_crossfit_folds(data, 4L)
  changed <- data
  for (name in names(data)) {
    partition <- if (name == "t") folds$target_folds else folds$source_folds[[name]]
    held_out <- unlist(lapply(partition[1:2], `[[`, "original_idx"))
    changed[[name]]$W_outcome[held_out, ] <- changed[[name]]$W_outcome[held_out, ] + 4
    changed[[name]]$Z_site[held_out, ] <- changed[[name]]$Z_site[held_out, ] + 4
    changed[[name]]$Y[held_out] <- changed[[name]]$Y[held_out] + 20
  }
  # Target outcomes are not inputs to the standard source program, even on
  # the allowed training rows. Only its feature moment is sent to the source.
  changed$t$Y <- changed$t$Y + 7
  changed_folds <- build_crossfit_folds(changed, 4L)
  fit <- function(partition) .fit_complete_source_program(partition$target_folds,
    partition$source_folds$s1, "s1", 3:4, "one_round", 1L, "gaussian", 12, 6L,
    1000L, "min", tol = 1e-10, calibration_control = list(source_nuisance_method = "standard"))
  original <- fit(folds)
  perturbed <- fit(changed_folds)
  expect_identical(as.numeric(original$weight), as.numeric(perturbed$weight))
  expect_identical(as.numeric(original$outcome), as.numeric(perturbed$outcome))
  expect_error(.fit_standard_source_nuisances(folds$source_folds$s1, "s1", 2:3,
    .prepare_source_training_message(folds$target_folds, 3:4), 1, 0, 0, 12, 6, 1000, "min"),
    "different folds")
})

test_that("constant outcome derivatives give the same paired weight objective", {
  set.seed(4603)
  site <- function() {
    x <- matrix(rnorm(1800), 600, 3)
    list(n = 600L, W_outcome = x, Z_site = x, A = rep(0:1, 300), Y = x[, 1] + rnorm(600))
  }
  folds <- build_crossfit_folds(list(t = site(), s1 = site()), 4)
  train <- function(program) .fit_complete_source_program(folds$target_folds,
    folds$source_folds$s1, "s1", 3:4, "one_round", 1, "gaussian", 12, 8, 1000, "min",
    tol = 1e-10, calibration_control = list(recipe = "score_derivative", source_nuisance_method = program))
  calibrated <- train("calibrated")
  standard <- train("standard")
  expect_identical(attr(calibrated$weight, "cv_seed"), attr(standard$weight, "cv_seed"))
  expect_equal(attr(calibrated$weight, "lambda_used"), attr(standard$weight, "lambda_used"), tolerance = 1e-12)
  expect_equal(as.numeric(calibrated$weight), as.numeric(standard$weight), tolerance = 1e-7)
})

test_that("source program controls distinguish complete validation from initial-model validation", {
  calibrated <- .validate_calibration_control()
  standard <- .validate_calibration_control(list(source_nuisance_method = "standard"))
  expect_identical(calibrated$source_nuisance_method, "calibrated")
  expect_identical(.match_source_validation("calibrated", calibrated, "one_round"), "calibrated")
  expect_identical(.match_source_validation("complete", calibrated, "one_round"), "complete")
  expect_identical(.match_source_validation("calibrated", standard, "one_round"), "complete")
  expect_error(.match_source_validation("initial", standard, "one_round"), "complete")
  expect_error(.match_source_validation("complete", standard, "two_round"), "one_round")
})
