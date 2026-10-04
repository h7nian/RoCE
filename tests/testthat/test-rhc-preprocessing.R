rhc_preprocessing_fixture <- function() {
  raw <- load_rhc_raw()
  arguments <- list(raw = raw, outcome = "death30", site_var = "ninsclas",
    site_recode = c("No insurance" = NA_character_, "Medicare & Medicaid" = NA_character_),
    covariate_profile = "historical")
  cohort <- do.call(build_rhc_cohort, arguments)
  split_arguments <- list(K = 4L, target_site = "Private", phi = base::identity, seed = 42L)
  data <- do.call(build_rhc_data_split, c(list(cohort = cohort), split_arguments))
  list(raw = raw, arguments = arguments, cohort = cohort, split_arguments = split_arguments,
    data = data, folds = build_crossfit_folds(data, 10L))
}

test_that("RHC reference rows control imputation and reject invalid boundaries", {
  expect_equal(.rhc_impute_continuous(c(1, 3, NA, 100), training_rows = 1:3), c(1, 3, 2, 100))
  expect_equal(.rhc_impute_mode(c("a", "b", NA, "b", "b"), training_rows = 1:3),
    c("a", "b", "a", "b", "b"))
  for (rows in list(integer(), c(1, 1), c(0, 2), c(1, 5), c(1, NA), c(1, 1.5))) {
    expect_error(.validate_rhc_preprocessing_rows(rows, 4), "training rows")
  }
})

test_that("all-row preprocessing reproduces historical RHC values and RNG", {
  withr::local_collate("C.UTF-8")
  x <- rhc_preprocessing_fixture()
  expect_identical(digest::digest(x$data, algo = "sha256"),
    "a6d8135c01cd4c1113e134b024b2b440a2e31c3e21e2e0a41b9b4d14d6c5aa1a")
  set.seed(711)
  state <- .Random.seed
  transformed <- .rhc_training_preprocessed_data(seq_len(nrow(x$cohort)), x$arguments, x$split_arguments)
  expect_identical(.Random.seed, state)
  for (site in names(x$data)) expect_identical(transformed[[site]], x$data[[site]])
  expect_identical(attr(transformed, "feature_center"), attr(x$data, "feature_center"))
  expect_identical(attr(transformed, "feature_scale"), attr(x$data, "feature_scale"))
})

test_that("held-out RHC covariates cannot change outer training transformations", {
  withr::local_collate("C.UTF-8")
  x <- rhc_preprocessing_fixture()
  first <- .prepare_rhc_outer_preprocessing(x$cohort, x$data, x$folds, x$arguments, x$split_arguments)
  rows <- .rhc_cohort_rows_by_site(x$cohort, x$data)
  site_folds <- c(list(t = x$folds$target_folds), x$folds$source_folds)
  excluded <- unlist(lapply(names(rows), function(site) rows[[site]][site_folds[[site]][[1L]]$original_idx]))
  retained_raw <- which(!x$raw$ninsclas %in% names(x$arguments$site_recode))
  changed <- x$arguments
  changed$raw$urin1[retained_raw[excluded]] <- 1000000
  changed$raw$race[retained_raw[excluded]] <- "other"
  changed_cohort <- do.call(build_rhc_cohort, changed)
  changed_data <- do.call(build_rhc_data_split, c(list(cohort = changed_cohort), x$split_arguments))
  second <- .prepare_rhc_outer_preprocessing(changed_cohort, changed_data, x$folds, changed, x$split_arguments)
  expect_identical(attr(first$folds, "preprocessing")$parameters[[1]],
    attr(second$folds, "preprocessing")$parameters[[1]])
  for (site in names(rows)) {
    get <- function(result) if (site == "t") result$folds$target_folds else result$folds$source_folds[[site]]
    a <- attr(.outer_fold_views(get(first), 1), ".data_ref")
    b <- attr(.outer_fold_views(get(second), 1), ".data_ref")
    training <- setdiff(seq_len(a$n), site_folds[[site]][[1L]]$original_idx)
    expect_identical(a$W_outcome[training, ], b$W_outcome[training, ])
    expect_identical(a$W_outcome, a$Z_site)
    expect_false(identical(a$W_outcome, b$W_outcome))
    expect_identical(lapply(get(first), `[[`, "original_idx"), lapply(get(second), `[[`, "original_idx"))
  }
  callback <- attr(first$baseline_data$t, "training_preprocessor")
  state <- .Random.seed
  callback(setdiff(seq_len(x$data$t$n), x$folds$target_folds[[1]]$original_idx))
  expect_identical(.Random.seed, state)
})

test_that("outer fold views preserve score coordinates and training independence", {
  set.seed(716L)
  data <- split_data_by_site(generate_bounded_data(n_target = 600L,
    n_source_sizes = 600L, K = 1L, p = 4L))
  folds <- build_crossfit_folds(data, 3L)
  run <- function(views) run_tate_crossfit(data, n_folds = 3L, precomputed_folds = views,
    communication_mode = "one_round", target_nuisance_method = "hou_calibrated",
    source_validation_method = "outer_fit", crossfit_layers = 2L, aggregation_mode = "joint_tate",
    calibration_control = list(recipe = "score_derivative", target_radius = 5),
    M_tau = 3, M_tau_inference = 3, nlambda_init = 10L, nuisance_tol = 1e-10,
    nuisance_solver = "proximal_newton", n_cores = 1L, verbose = FALSE)
  original <- run(folds)
  no_op <- folds
  attr(no_op$target_folds, ".outer_preprocessed_views") <- rep(list(folds$target_folds), 3)
  attr(no_op$source_folds$s1, ".outer_preprocessed_views") <- rep(list(folds$source_folds$s1), 3)
  identical_fit <- run(no_op)
  expect_equal(identical_fit$all_phi_tau, original$all_phi_tau, tolerance = 1e-12)
  expect_equal(identical_fit$se, original$se, tolerance = 1e-12)
  expect_equal(identical_fit$fold_weights, original$fold_weights, tolerance = 1e-12)
  changed <- no_op
  changed_target <- data$t
  selected <- folds$target_folds[[1]]$original_idx
  changed_target$W_outcome[selected, 1] <- changed_target$W_outcome[selected, 1] + 1
  changed_target$Z_site <- changed_target$W_outcome
  changed_target$Y[selected] <- 1L
  attr(changed$target_folds, ".outer_preprocessed_views")[[1]] <-
    .rhc_replace_fold_data(folds$target_folds, changed_target)
  perturbed <- run(changed)
  for (arm in c("mu1", "mu0")) {
    a <- original$arm_results[[arm]]$fold_results[[1]]
    b <- perturbed$arm_results[[arm]]$fold_results[[1]]
    expect_equal(as.numeric(a$source_results$s1$alpha_ts), as.numeric(b$source_results$s1$alpha_ts))
    expect_equal(as.numeric(a$source_results$s1$gamma_s), as.numeric(b$source_results$s1$gamma_s))
  }
  expect_equal(perturbed$fold_weights[1, ], original$fold_weights[1, ], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(perturbed$estimate, original$estimate, tolerance = 1e-8)))
  expect_error(reaggregate_tate_crossfit(data, list(preprocessing = list(mode = "outer_fold")), 3),
    "recorded outer preprocessing")
})
