library(testthat)

test_that("nuisance CV group IDs are validated without coercion", {
  expect_null(.validate_nuisance_cv_group_id(NULL, 3L, "test"))
  expect_equal(
    .validate_nuisance_cv_group_id(c(1, 1, 2), 3L, "test"),
    c(1L, 1L, 2L)
  )
  for (bad in list(
    c(1, 2), c(1, 2.5, 3), c(0, 1, 2), c(-1, 1, 2),
    c(1, NA, 2), c(1, Inf, 2), c("1", "1", "2"),
    c(1, 2, .Machine$integer.max + 1), c(1 + 0i, 2 + 0i, 3 + 0i)
  )) {
    expect_error(
      .validate_nuisance_cv_group_id(bad, 3L, "test"),
      "cv_group_id.*positive finite integer vector"
    )
  }
  for (bad_n in list(0L, 2.5, NA_real_, Inf)) {
    expect_error(
      .validate_nuisance_cv_group_id(1:3, bad_n, "test"),
      "n must be one positive integer"
    )
  }
})

test_that("NULL and unique nuisance groups consume no RNG", {
  assert_rng_neutral <- function(groups) {
    set.seed(810L)
    expected <- runif(1L)
    set.seed(810L)
    result <- .make_nuisance_cv_fold_id(groups, 5L, "test")
    observed <- runif(1L)
    expect_null(result)
    expect_identical(observed, expected)
  }
  assert_rng_neutral(NULL)
  assert_rng_neutral(1:20)
})

test_that("repeated nuisance groups remain intact in near-equal group folds", {
  groups <- rep(c(101L, 205L, 307L, 409L, 511L, 613L, 719L),
                c(3L, 1L, 4L, 2L, 5L, 1L, 2L))
  set.seed(44L)
  fold_id <- .make_nuisance_cv_fold_id(groups, 3L, "test")

  expect_length(fold_id, length(groups))
  expect_setequal(unique(fold_id), 1:3)
  expect_true(all(vapply(split(fold_id, groups), function(x) {
    length(unique(x)) == 1L
  }, logical(1L))))
  groups_per_fold <- table(vapply(split(fold_id, groups), `[[`, integer(1L), 1L))
  expect_lte(max(groups_per_fold) - min(groups_per_fold), 1L)

  expect_error(
    .make_nuisance_cv_fold_id(c(1L, 1L, 2L, 2L), 3L, "test"),
    "only 2 unique groups for 3 folds"
  )
  for (bad_folds in list(1L, 2.5, NA_real_, Inf)) {
    expect_error(
      .make_nuisance_cv_fold_id(groups, bad_folds, "test"),
      "n_folds must be one integer >= 2"
    )
  }
})

test_that("fit_glmnet_cv unique groups preserve the exact legacy path", {
  skip_if_not_installed("glmnet")
  set.seed(72L)
  x <- matrix(rnorm(80L * 4L), 80L, 4L)
  y <- as.numeric(x[, 1L] + rnorm(80L) > 0)
  x_new <- x[1:7, , drop = FALSE]

  set.seed(901L)
  legacy <- fit_glmnet_cv(
    x, y, x_new, family = "binomial", nlambda = 7L,
    min_per_fold = 10L, lambda_rule = "min"
  )
  set.seed(901L)
  unique_group <- fit_glmnet_cv(
    x, y, x_new, family = "binomial", nlambda = 7L,
    min_per_fold = 10L, lambda_rule = "min", cv_group_id = seq_len(nrow(x))
  )
  expect_identical(unique_group, legacy)
})

test_that("fit_glmnet_cv fails before fitting malformed groups", {
  x <- matrix(seq_len(40), 10L, 4L)
  y <- rep(0:1, 5L)
  expect_error(
    fit_glmnet_cv(
      x, y, x[1:2, , drop = FALSE], family = "binomial",
      cv_group_id = rep(1L, 9L)
    ),
    "cv_group_id.*length 10"
  )
})

test_that("target complement wires PS and arm-subset groups and isolates caches", {
  n_train <- 20L
  train <- list(
    W_outcome = matrix(seq_len(n_train * 4L), n_train, 4L),
    Z_site = matrix(seq_len(n_train * 4L), n_train, 4L),
    A = rep(0:1, n_train / 2L),
    Y = rep(c(0, 0, 0, 1), length.out = n_train),
    n = n_train,
    original_idx = seq_len(n_train),
    cv_group_id = rep(seq_len(10L), each = 2L)
  )
  evaluation <- list(
    W_outcome = matrix(seq_len(16L), 4L, 4L),
    Z_site = matrix(seq_len(16L), 4L, 4L),
    A = c(0, 1, 0, 1), Y = c(0, 1, 1, 0), n = 4L,
    original_idx = 1:4
  )
  calls <- list()
  fake_fit <- function(x_train, y_train, x_predict, ..., model_name,
                       cv_group_id = NULL) {
    calls[[length(calls) + 1L]] <<- list(
      model = model_name, groups = cv_group_id, n = nrow(x_train)
    )
    rep(0.5, nrow(x_predict))
  }
  cache <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    combine_folds = function(...) train,
    materialize_fold = function(...) evaluation,
    fit_glmnet_cv = fake_fit,
    .package = "RoCE"
  )

  estimate_target_only_from_complement(
    target_folds = vector("list", 3L), k1 = 1L, n_folds = 3L,
    family = "binomial", A_val = 1L, propensity_cache = cache,
    nuisance_lambda_rule = "min"
  )
  expect_identical(calls[[1L]]$model, "PS")
  expect_identical(calls[[1L]]$groups, train$cv_group_id)
  expect_identical(calls[[2L]]$model, "OR")
  expect_identical(calls[[2L]]$groups, train$cv_group_id[train$A == 1L])

  # A different grouping vector must not reuse the first propensity cache.
  train$cv_group_id <- rep(seq_len(5L), each = 4L)
  estimate_target_only_from_complement(
    target_folds = vector("list", 3L), k1 = 1L, n_folds = 3L,
    family = "binomial", A_val = 1L, propensity_cache = cache,
    nuisance_lambda_rule = "min"
  )
  expect_identical(vapply(calls, `[[`, character(1L), "model"),
                   c("PS", "OR", "PS", "OR"))
})

test_that("target complement default cache behavior remains unchanged", {
  n_train <- 20L
  train <- list(
    W_outcome = matrix(seq_len(n_train * 4L), n_train, 4L),
    Z_site = matrix(seq_len(n_train * 4L), n_train, 4L),
    A = rep(0:1, n_train / 2L),
    Y = rep(c(0, 0, 0, 1), length.out = n_train),
    n = n_train, original_idx = seq_len(n_train)
  )
  evaluation <- list(
    W_outcome = matrix(seq_len(16L), 4L, 4L),
    Z_site = matrix(seq_len(16L), 4L, 4L),
    A = c(0, 1, 0, 1), Y = c(0, 1, 1, 0), n = 4L,
    original_idx = 1:4
  )
  models <- character()
  testthat::local_mocked_bindings(
    combine_folds = function(...) train,
    materialize_fold = function(...) evaluation,
    fit_glmnet_cv = function(x_train, y_train, x_predict, ..., model_name,
                             cv_group_id = NULL) {
      expect_null(cv_group_id)
      models <<- c(models, model_name)
      rep(0.5, nrow(x_predict))
    },
    .package = "RoCE"
  )
  cache <- new.env(parent = emptyenv())
  for (iteration in 1:2) {
    estimate_target_only_from_complement(
      vector("list", 3L), 1L, 3L, family = "binomial", A_val = 1L,
      propensity_cache = cache, nuisance_lambda_rule = "min"
    )
  }
  # PS is cached on the second call; OR is deliberately refit each time.
  expect_identical(models, c("PS", "OR", "OR"))
})
