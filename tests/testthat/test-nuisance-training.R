test_that("exact training-subset reuse preserves fits and the caller RNG", {
  set.seed(182)
  x <- matrix(rnorm(600), 200, 3)
  arguments <- list(W_outcome = x, Y = x[, 1] + rnorm(200),
                    A = rep(1, 200), A_val = 1L, lambda = NULL,
                    nlambda = 10L, family = "gaussian", lambda_rule = "min")
  cache <- new.env(parent = emptyenv())
  state <- .Random.seed
  fit <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "t", c(2L, 4L), cache)
  expect_identical(.Random.seed, state)
  cached <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "t", c(4L, 2L), cache)
  uncached <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "t", c(2L, 4L))
  expect_identical(cached, fit)
  expect_identical(uncached, fit)

  # A matching fold label is insufficient when its actual training data change.
  arguments$Y <- arguments$Y + .5
  changed <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "t", c(2L, 4L), cache)
  changed_uncached <- .fit_nuisance_training_subset("fit_initial_outcome", arguments, "t", c(2L, 4L))
  expect_identical(changed, changed_uncached)
  expect_gt(abs(changed[1] - fit[1]), .4)
  expect_identical(.Random.seed, state)
})

test_that("unsupported RCAL combinations and invalid penalties fail explicitly", {
  expect_error(fit_site_aipw(list(), use_rcal = TRUE, use_crossfit = TRUE),
               "not supported", fixed = TRUE)
  expect_error(fit_initial_outcome(matrix(0, 20, 2), rep(0, 20), rep(1, 20),
                                  lambda = NA_real_), "finite numeric scalar", fixed = TRUE)
  expect_error(.validate_nuisance_cache_flag(NA, "test"), "TRUE or FALSE", fixed = TRUE)
})

test_that("failed nuisance fits retain their site, arm, fold and seed context", {
  failure <- tryCatch(.fit_nuisance_training_subset("fit_initial_density_ratio", list(
    Z_site = matrix(0, 20, 2), A = rep(1, 20), mean_phi = c(1, 0, 0),
    A_val = 1L, lambda = NA_real_), "s1", c(3L, 2L)), nuisance_fit_error = identity)
  expect_s3_class(failure, "nuisance_fit_error")
  expect_identical(failure$nuisance_context$site, "s1")
  expect_identical(failure$nuisance_context$training_folds, 2:3)
  expect_s3_class(failure$parent, "error")
})

test_that("inner initial nuisances exclude their evaluation fold in both protocols and arms", {
  skip_on_cran()
  set.seed(314)
  generated <- generate_simulation_data(
    n_total = 2000L, K = 1L, p = 4L, config = "C1", dgp_type = "face",
    outcome_type = "continuous", n_target = 1000L, n_source_sizes = 1000L,
    ate_deviation = 0, n_deviated_sites = 0L, warn_ignored = FALSE
  )
  data <- split_data_by_site(generated)
  folds <- build_crossfit_folds(data, 3L)
  changed <- data
  for (site in names(data)) {
    site_folds <- if (site == "t") folds$target_folds else folds$source_folds[[site]]
    evaluation_rows <- site_folds[[3L]]$original_idx
    changed[[site]]$Y[evaluation_rows] <- rev(changed[[site]]$Y[evaluation_rows])
    changed[[site]]$W_outcome[evaluation_rows, 1L] <-
      changed[[site]]$W_outcome[evaluation_rows, 1L] + .2
    changed[[site]]$Z_site[evaluation_rows, 1L] <-
      changed[[site]]$Z_site[evaluation_rows, 1L] + .2
  }
  changed_folds <- build_crossfit_folds(changed, 3L)
  for (mode in c("one_round", "two_round")) for (arm in 0:1) {
    fit <- function(split, partition, use_cache) {
      suppressWarnings(run_crossfit(
        split, n_folds = 3L, communication_mode = mode, A_val = arm,
        nlambda_init = 8L, family = "gaussian", verbose = FALSE, n_cores = 1L,
        precomputed_folds = partition, use_lambda_cache = use_cache
      ))
    }
    original_fit <- fit(data, folds, TRUE)
    uncached_fit <- fit(data, folds, FALSE)
    changed_fit <- fit(changed, changed_folds, TRUE)
    original_source <- original_fit$fold_results[[1L]]$source_results$s1
    changed_source <- changed_fit$fold_results[[1L]]$source_results$s1
    for (field in c("per_k2_alpha", "per_k2_gamma")) {
      before <- original_source[[field]]$k2_3
      after <- changed_source[[field]]$k2_3
      expect_identical(as.numeric(after), as.numeric(before))
      expect_identical(attr(after, "lambda_used"), attr(before, "lambda_used"))
      expect_identical(attr(after, "cv_seed"), attr(before, "cv_seed"))
    }
    expect_identical(original_fit$estimate, uncached_fit$estimate)
    expect_identical(original_fit$se, uncached_fit$se)
    expect_identical(original_fit$use_lambda_cache, TRUE)
    expect_identical(original_fit$nuisance_training_policy, NUISANCE_TRAINING_POLICY)
  }
})

test_that("rho reuse retains cache policy and rejects legacy training metadata", {
  skip_on_cran()
  data <- get_smoke_data_split()
  fit <- suppressWarnings(run_tate_crossfit(
    data, n_folds = 3L, communication_mode = "one_round", nlambda_init = 5L,
    family = "gaussian", verbose = FALSE, n_cores = 1L, use_lambda_cache = FALSE
  ))
  changed <- data
  changed$s1$Y[changed$s1$A == 1] <- changed$s1$Y[changed$s1$A == 1] + .1
  reused <- suppressWarnings(.reuse_one_round_tate_across_rho(data, changed, fit, "s1", 1))
  independent <- suppressWarnings(run_tate_crossfit(
    changed, n_folds = 3L, communication_mode = "one_round", nlambda_init = 5L,
    family = "gaussian", verbose = FALSE, n_cores = 1L, use_lambda_cache = FALSE
  ))
  expect_identical(reused$use_lambda_cache, FALSE)
  expect_equal(reused$estimate, independent$estimate, tolerance = 1e-12)
  expect_equal(reused$se, independent$se, tolerance = 1e-12)
  legacy_arm <- fit$arm_results$mu1
  legacy_arm$nuisance_training_policy <- NULL
  expect_error(.refit_one_round_crossfit_sources(data, changed, legacy_arm, "s1"),
               "nuisance_training_policy", fixed = TRUE)
})
