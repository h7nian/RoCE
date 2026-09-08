make_grouped_crossfit_fixture <- function() {
  set.seed(7103)
  sites <- lapply(seq_len(3L), function(j) {
    origin <- rep(seq_len(300L), each = 2L)
    x <- matrix(rnorm(600), ncol = 2L)
    treatment <- rep(0:1, length.out = 300L)
    y <- 0.5 * treatment + 0.3 * x[, 1L] - 0.2 * x[, 2L] + rnorm(300)
    list(n = 600L, X = x[origin, ], X_dagger = x[origin, ],
         W_outcome = x[origin, ], Z_site = x[origin, ],
         W_outcome_true = x[origin, ], Z_site_true = x[origin, ],
         A = treatment[origin], Y = y[origin], cv_group_id = origin)
  })
  names(sites) <- c("t", "s1", "s2")
  views <- lapply(sites, function(site) {
    fold_id <- (site$cv_group_id - 1L) %% 3L + 1L
    folds <- lapply(1:3, function(k) {
      idx <- which(fold_id == k)
      list(original_idx = idx, n = length(idx))
    })
    attr(folds, ".data_ref") <- site
    folds
  })
  list(data = sites,
       folds = list(target_folds = views$t, source_folds = views[c("s1", "s2")]))
}

test_that("fold views preserve nuisance group IDs, including empty combinations", {
  fixture <- make_grouped_crossfit_fixture()
  views <- fixture$folds$target_folds
  one <- materialize_fold(views, 2L)
  expect_identical(one$cv_group_id, fixture$data$t$cv_group_id[one$original_idx])
  combined <- combine_folds(views, c(3L, 1L))
  expect_identical(combined$cv_group_id,
                   fixture$data$t$cv_group_id[combined$original_idx])
  expect_identical(combine_folds(views, integer(0))$cv_group_id, integer(0))
  expect_silent(.validate_nuisance_cv_fold_views(fixture$data$t, views, "test"))
  changed <- fixture$data$t
  changed$cv_group_id <- seq_len(changed$n)
  expect_error(.validate_nuisance_cv_fold_views(changed, views, "test"), "different")
  bad <- views
  bad[[2L]]$original_idx[1L] <- bad[[1L]]$original_idx[1L]
  expect_error(.validate_nuisance_cv_fold_views(fixture$data$t, bad, "test"), "partition")
  changed <- fixture$data$t
  changed$cv_group_id[views[[2L]]$original_idx[1L]] <- changed$cv_group_id[views[[1L]]$original_idx[1L]]
  attr(views, ".data_ref") <- changed
  expect_error(.validate_nuisance_cv_fold_views(changed, views, "test"), "crosses outer")
  expect_error(.validate_nuisance_cv_group_id(matrix(1:6, ncol = 1L), 6L, "test"), "vector")
})

test_that("all target and source nuisance tuning stages receive intact origin groups", {
  skip_on_cran()
  fixture <- make_grouped_crossfit_fixture()
  captured <- list()
  original <- .make_nuisance_cv_fold_id
  local_mocked_bindings(.make_nuisance_cv_fold_id = function(cv_group_id, n_folds, caller) {
    result <- original(cv_group_id, n_folds, caller)
    captured[[length(captured) + 1L]] <<- list(groups = cv_group_id,
                                             folds = result, caller = caller)
    result
  }, .package = "RoCE")
  for (mode in c("one_round", "two_round")) {
    fitted <- run_tate_crossfit(fixture$data, n_folds = 3L,
      communication_mode = mode, family = "gaussian", nlambda_init = 8L,
      n_cores = 1L, verbose = FALSE, precomputed_folds = fixture$folds)
    expect_true(is.finite(fitted$estimate) && is.finite(fitted$se) && fitted$se > 0)
    expect_identical(dim(fitted$fold_weights), c(3L, 2L))
  }
  expect_true(length(captured) > 0L)
  expect_setequal(unique(vapply(captured, `[[`, character(1L), "caller")),
    c("estimate_complement_fold_aipw PS", "estimate_complement_fold_aipw OR",
      "fit_initial_outcome", "fit_initial_density_ratio",
      "fit_unified_density_ratio", "fit_unified_outcome"))
  for (call in captured) {
    expect_false(is.null(call$groups))
    expect_false(is.null(call$folds))
    expect_true(all(vapply(split(call$folds, call$groups),
                          function(x) length(unique(x)) == 1L, logical(1L))))
  }
})

test_that("all-unique group metadata leaves the complete default TATE path unchanged", {
  skip_on_cran()
  fixture <- make_grouped_crossfit_fixture()
  make_variant <- function(unique_groups) {
    data <- lapply(fixture$data, function(site) {
      site$cv_group_id <- if (unique_groups) seq_len(site$n) else NULL
      site
    })
    folds <- fixture$folds
    attr(folds$target_folds, ".data_ref") <- data$t
    for (site in names(folds$source_folds)) attr(folds$source_folds[[site]], ".data_ref") <- data[[site]]
    run_tate_crossfit(data, n_folds = 3L, communication_mode = "one_round",
      family = "gaussian", nlambda_init = 8L, n_cores = 1L,
      verbose = FALSE, precomputed_folds = folds)
  }
  ordinary <- make_variant(FALSE)
  unique <- make_variant(TRUE)
  expect_equal(unique$estimate, ordinary$estimate, tolerance = 1e-12)
  expect_equal(unique$se, ordinary$se, tolerance = 1e-12)
  expect_equal(unique$fold_weights, ordinary$fold_weights, tolerance = 1e-12)
  expect_equal(unique$all_phi_tau, ordinary$all_phi_tau, tolerance = 1e-12)
})
