test_that("plugin block design reproduces fold-specific linear predictors", {
  skip_if_not(exists(".make_plugin_block_design"), message = "internal helper not loaded")

  X1 <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  X2 <- matrix(c(5, 6), nrow = 1)
  beta1 <- c(0.5, 1.0, -0.25)
  beta2 <- c(-1.0, 0.2, 0.3)

  block <- .make_plugin_block_design(list(X1, X2), caller = "test")
  beta_block <- c(0, beta1, beta2)
  pred <- drop(cbind(1, block) %*% beta_block)
  expected <- c(drop(cbind(1, X1) %*% beta1),
                drop(cbind(1, X2) %*% beta2))

  expect_equal(pred, expected)
})

test_that("cross-fitting code uses fold-summed calibration, not parameter averaging", {
  path <- test_path("../../R/cross_fitting_algorithms.R")
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(grepl("Fold-summed SMMAL-style calibrated optimization", txt, fixed = TRUE))
  expect_false(grepl("gamma_final_k1 <- colMeans(do.call(rbind, gamma_cal_list))", txt, fixed = TRUE))
  expect_false(grepl("alpha_final_k1 <- colMeans(do.call(rbind, alpha_cal_list))", txt, fixed = TRUE))
  expect_false(grepl("mean_phi <- c(1, colMeans(target_calib_k2$Z_site))", txt, fixed = TRUE))
  expect_false(grepl("mean_phi_cache[[k2_key]] <- c(1, colMeans(target_fold_cache[[k2_key]]$Z_site))", txt, fixed = TRUE))
})

test_that("aggregation lambda cv evaluates held-out inner folds", {
  skip_if_not(exists("select_aggregation_lambda_inner_cv"), message = "inner CV selector not loaded")
  skip_if_not(exists("optimize_weights_cpp"), message = "C++ weight optimizer not compiled")

  make_component <- function(mu_t, mu_s, v_shift = 0) {
    list(
      V_ot = 1.0 + v_shift,
      V_t = c(0.8 + v_shift, 0.9 + v_shift),
      V_s = c(0.7 + v_shift, 0.6 + v_shift),
      C_ot = c(0.05, 0.04),
      C_cross = matrix(c(0, 0.02, 0.02, 0), 2, 2),
      avg_target_est = mu_t,
      avg_source_est = mu_s,
      n_t = 40,
      n_s = c(45, 50)
    )
  }

  components <- list(
    make_component(0.50, c(0.49, 0.62), 0.00),
    make_component(0.52, c(0.50, 0.65), 0.05),
    make_component(0.48, c(0.47, 0.61), 0.02)
  )

  lambda <- select_aggregation_lambda_inner_cv(
    components,
    lambda_grid = c(0.001, 0.01, 0.1),
    lambda_rule = "min"
  )

  expect_true(is.finite(lambda))
  expect_true(lambda %in% c(0.001, 0.01, 0.1))
})

test_that("aggregation lambda grid validates user input", {
  skip_if_not(exists(".aggregation_lambda_grid"), message = "aggregation grid helper not loaded")

  component <- list(
    V_t = c(0.8, 0.9),
    V_s = c(0.7, 0.6),
    avg_target_est = 0.5,
    avg_source_est = c(0.49, 0.62),
    n_t = 40,
    n_s = c(45, 50)
  )

  expect_error(
    .aggregation_lambda_grid(component, lambda_grid = c(0.01, NA_real_)),
    "positive finite"
  )
  expect_warning(
    clipped <- .aggregation_lambda_grid(component, lambda_grid = c(LAMBDA_MIN / 10, LAMBDA_MAX * 10)),
    "clipped"
  )
  expect_equal(clipped, c(LAMBDA_MIN, LAMBDA_MAX))
})

test_that("aggregation default lambda grid starts at KKT lambda_max", {
  skip_if_not(exists(".aggregation_lambda_grid"), message = "aggregation grid helper not loaded")
  skip_if_not(exists(".aggregation_lambda_max"), message = "aggregation lambda_max helper not loaded")

  component <- list(
    V_ot = 1.0,
    C_ot = c(0.05, 0.04),
    avg_target_est = 0.50,
    avg_source_est = c(0.40, 0.70),
    n_t = 50
  )

  expected_lambda_max <- max(abs(2 * (component$C_ot - component$V_ot) / component$n_t) /
                               (component$avg_target_est - component$avg_source_est)^2)
  lambda_max <- .aggregation_lambda_max(component)
  grid <- .aggregation_lambda_grid(component)

  expect_equal(lambda_max, expected_lambda_max)
  expect_equal(unname(grid[1]), expected_lambda_max)
  expect_equal(attr(grid, "lambda_max"), expected_lambda_max)
  expect_true(all(diff(grid) < 0))
})

test_that("aggregation lambda_max warns for unpenalized active coordinates", {
  skip_if_not(exists(".aggregation_lambda_max"), message = "aggregation lambda_max helper not loaded")

  component <- list(
    V_ot = 1.0,
    C_ot = c(0.05, 0.04),
    avg_target_est = 0.50,
    avg_source_est = c(0.50, 0.70),
    n_t = 50
  )

  expect_warning(
    lambda_max <- .aggregation_lambda_max(component),
    "zero aggregation penalty weight"
  )
  expect_true(is.finite(lambda_max))
})

test_that("fold partitioning fails when stratified assignment is impossible", {
  data <- list(
    W_outcome = matrix(rnorm(20), nrow = 10L),
    Z_site = matrix(rnorm(20), nrow = 10L),
    A = c(rep(1L, 9L), 0L),
    Y = rnorm(10L),
    n = 10L
  )

  expect_error(
    partition_into_folds(data, n_folds = 3L, seed = 1L),
    "too few treated or control units"
  )
})

test_that("cross-fit aggregation calls the inner validation selector", {
  path <- test_path("../../R/cross_fitting_aggregation.R")
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")

  expect_true(grepl("select_aggregation_lambda_inner_cv", txt, fixed = TRUE))
  expect_true(grepl(".validation_aggregation_objective", txt, fixed = TRUE))
  expect_true(grepl("weight optimization failed for lambda index", txt, fixed = TRUE))
  expect_false(grepl("error = function(e) rep(NA_real_, K)", txt, fixed = TRUE))
})
