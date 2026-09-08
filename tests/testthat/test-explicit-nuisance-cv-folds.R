cv_fold_test_calls <- function() {
  set.seed(20260905)
  n <- 48L
  p <- 3L
  x <- matrix(rnorm(n * p), n, p)
  a <- rep(c(0, 1), length.out = n)
  y <- rbinom(n, 1, plogis(x[, 1] - 0.25 * x[, 2]))
  lambda <- c(0.5, 0.1, 0.02)
  zero <- numeric(p + 1L)
  arm_n <- sum(a == 1)
  # Deliberately unequal groups exercise the documented equal-fold-score
  # averaging behavior without asking C++ to rebalance central group IDs.
  fold_id <- rep(1:3, times = c(6L, 8L, arm_n - 14L))

  calls <- list(
    refined_density = function(id) select_lambda_cv_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), zero, lambda, 3L, 500L, 1e-4,
      1L, 1L, 1L, cv_fold_id = id
    ),
    initial_density = function(id) select_lambda_cv_initial_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), lambda, 3L, 500L, 1e-4, 1L, 5,
      cv_fold_id = id
    ),
    calibrated_density = function(id) select_lambda_cv_calibrated_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), zero, lambda, 3L, 500L, 1e-4,
      5, x, 1L, 1L, 1L, cv_fold_id = id
    ),
    refined_outcome = function(id) select_lambda_cv_general_refined_outcome_cpp(
      x, y, a, zero, 1L, 1L, lambda, 3L, 500L, 1e-4, 1L, x,
      cv_fold_id = id
    ),
    calibrated_outcome = function(id) select_lambda_cv_calibrated_outcome_cpp(
      x, y, a, zero, lambda, 3L, 500L, 1e-4, 1L, 5, x, 1L, 1L,
      cv_fold_id = id
    )
  )
  omitted <- list(
    refined_density = function() select_lambda_cv_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), zero, lambda, 3L, 500L, 1e-4,
      1L, 1L, 1L
    ),
    initial_density = function() select_lambda_cv_initial_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), lambda, 3L, 500L, 1e-4, 1L, 5
    ),
    calibrated_density = function() select_lambda_cv_calibrated_density_ratio_cpp(
      x, a, c(0.5, rep(0, p)), zero, lambda, 3L, 500L, 1e-4,
      5, x, 1L, 1L, 1L
    ),
    refined_outcome = function() select_lambda_cv_general_refined_outcome_cpp(
      x, y, a, zero, 1L, 1L, lambda, 3L, 500L, 1e-4, 1L, x
    ),
    calibrated_outcome = function() select_lambda_cv_calibrated_outcome_cpp(
      x, y, a, zero, lambda, 3L, 500L, 1e-4, 1L, 5, x, 1L, 1L
    )
  )
  list(calls = calls, omitted = omitted, fold_id = fold_id)
}

test_that("all nuisance selectors honor arm-filtered explicit CV folds", {
  fixture <- cv_fold_test_calls()
  for (call in fixture$calls) {
    set.seed(1)
    first <- call(fixture$fold_id)
    set.seed(999)
    second <- call(fixture$fold_id)
    expect_equal(second, first, tolerance = 0)
    expect_equal(first$cv_fold_id, fixture$fold_id)
    expect_equal(first$validation_fold_sizes, tabulate(fixture$fold_id, 3L))
  }
})

test_that("omitted and NULL explicit folds preserve the historical path", {
  fixture <- cv_fold_test_calls()
  for (name in names(fixture$calls)) {
    set.seed(42)
    explicit_null <- fixture$calls[[name]](NULL)
    set.seed(42)
    omitted <- fixture$omitted[[name]]()
    expect_equal(explicit_null, omitted, tolerance = 0)
    expect_false(any(c("cv_fold_id", "validation_fold_sizes") %in% names(omitted)))
  }
})

test_that("all nuisance selectors reject malformed explicit CV folds", {
  fixture <- cv_fold_test_calls()
  malformed <- list(
    fixture$fold_id[-1L],
    replace(fixture$fold_id, 1L, 1.5),
    replace(fixture$fold_id, 1L, NA_real_),
    replace(fixture$fold_id, 1L, Inf),
    replace(fixture$fold_id, 1L, 0),
    rep(1, length(fixture$fold_id)),
    factor(fixture$fold_id),
    matrix(fixture$fold_id, ncol = 1L)
  )
  for (call in fixture$calls) {
    for (id in malformed) {
      expect_error(call(id), "cv_fold_id")
    }
  }
})
