# Check interface semantics without fitting nuisance models or selecting effects.
test_that("RHC ablations pass only target-compatible calibration controls", {
  for (source in c("calibrated", "standard")) {
    calibrated <- .rhc_calibration_control("hou_calibrated", source, 5)
    ordinary <- .rhc_calibration_control("lasso", source, 5)
    expect_equal(calibrated$target_radius, 5)
    expect_false("target_radius" %in% names(ordinary))
    expect_identical(ordinary$source_nuisance_method, source)
    expect_silent(.validate_calibration_control(ordinary, "lasso"))
  }
  expect_error(.rhc_calibration_control("hou_calibrated", "calibrated", -1), "radius")
})

test_that("RHC covariate sensitivities preserve patients and historical defaults", {
  raw <- load_rhc_raw()
  excluded <- c("Medicare & Medicaid", "No insurance")
  recode <- setNames(rep(NA_character_, 2), excluded)
  historical <- build_rhc_cohort(raw, site_var = "ninsclas", site_recode = recode)
  reference <- RoCE::build_rhc_cohort(raw, site_var = "ninsclas", site_recode = recode)
  expect_identical(historical, reference)
  log_features <- build_rhc_cohort(raw, site_var = "ninsclas", site_recode = recode,
                                 covariate_profile = "log_missing")
  grouped <- build_rhc_cohort(raw, site_var = "ninsclas", site_recode = recode,
                             covariate_profile = "grouped_log_missing")
  for (cohort in list(log_features, grouped)) {
    expect_identical(cohort[c("A", "Y", "site_var")], historical[c("A", "Y", "site_var")])
    expect_equal(nrow(cohort), 5039L)
    expect_true(all(c("age", "das2d3pc", "aps1", "surv2md1") %in% names(cohort)))
    split <- build_rhc_data_split(cohort, K = 4L, target_site = "Private")
    expect_true(all(vapply(split, function(site) identical(site$W_outcome, site$Z_site), logical(1L))))
  }
  expect_equal(ncol(log_features)-3L, 62L)
  expect_equal(ncol(grouped)-3L, 60L)
  for (variable in c("urin1", "crea1", "bili1", "wblc1")) {
    expect_equal(log_features[[paste0("log1p_", variable)]], log1p(historical[[variable]]))
  }
  expect_equal(sum(log_features$urin1_missing), 2651L)
  expect_equal(grouped$injury_diagnosis, as.integer(historical$trauma == 1 | historical$ortho == 1))
  expect_false(any(c("ortho", "trauma", "cat1Colon.Cancer", "cat1Lung.Cancer") %in% names(grouped)))
  expect_equal(grouped$cat1Solid.Cancer, historical$cat1Colon.Cancer + historical$cat1Lung.Cancer)
  changed <- raw
  changed$dth30 <- ifelse(raw$dth30 == "Yes", "No", "Yes")
  alternative <- build_rhc_cohort(changed, site_var = "ninsclas", site_recode = recode,
                                 covariate_profile = "grouped_log_missing")
  expect_identical(alternative[setdiff(names(alternative), "Y")], grouped[setdiff(names(grouped), "Y")])
  changed <- raw
  changed$crea1[which(!raw$ninsclas %in% excluded)[1L]] <- -1
  expect_error(build_rhc_cohort(changed, site_var = "ninsclas", site_recode = recode,
                                covariate_profile = "log_missing"), "nonnegative")
})

rhc_interface_fixture <- function() {
  scope <- new.env(parent = environment(run_rhc_tate_experiment))
  runner <- run_rhc_tate_experiment
  environment(runner) <- scope
  site <- list(n = 60L, X = matrix(seq_len(120) / 120, 60L, 2L),
               A = rep(0:1, 30), Y = rep(c(0, 0, 1), 20))
  site$X_dagger <- site$Z_site <- site$W_outcome <- site$X
  data <- list(t = site, s1 = site, s2 = site)
  attr(data, "site_mapping") <- c(target = "Private", s1 = "Medicare", s2 = "Medicaid")
  attr(data, "site_components") <- attr(data, "site_mapping")
  scope$build_rhc_cohort <- function(...) NULL
  scope$build_rhc_data_split <- function(...) data
  scope$.target_tate_reference <- function(target_folds, family, nuisance_lambda_rule) {
    scope$reference_folds <- target_folds
    list(estimate = .10, se = .02)
  }
  scope$run_tate_crossfit <- function(...) {
    args <- list(...)
    scope$forwarded <- args
    common <- args$aggregation_mode == "common_tate"
    weights <- list(mu1 = c(s1 = .2, s2 = .3), mu0 = c(s1 = .1, s2 = -.1))
    coordinates <- if (common) c("s1", "s2") else c("mu1:s1", "mu1:s2", "mu0:s1", "mu0:s2")
    wald <- matrix(seq_len(3L * length(coordinates)), 3L, dimnames = list(NULL, coordinates))
    list(estimate = .12, se = .025, target_only = list(estimate = .15, se = .03),
         source_estimates = c(s1 = .11, s2 = .13),
         weights = if (common) weights$mu1 else setNames(unlist(weights, use.names = FALSE), coordinates),
         weights_by_arm = weights, fold_wald_statistics = wald,
         fold_penalty_coefficients = pmax(wald - 5, 0),
         crossfit_levels = if (is.null(args$crossfit_layers)) 2L else args$crossfit_layers,
         source_validation_method = args$source_validation_method,
         calibration_control = args$calibration_control, nuisance_solver = args$nuisance_solver,
         nuisance_cv_certificate = args$nuisance_cv_certificate)
  }
  list(run = runner, scope = scope)
}

test_that("rotating the RHC target preserves each insurance cohort and shared features", {
  # This saved fingerprint includes the historical categorical reference levels.
  withr::local_collate("C.UTF-8")
  excluded <- c("Medicare & Medicaid", "No insurance")
  cohort <- build_rhc_cohort(site_var = "ninsclas", site_recode = setNames(rep(NA_character_, 2L), excluded))
  reference <- build_rhc_data_split(cohort, K = 4L, target_site = "Private")
  mapping <- attr(reference, "site_mapping")
  by_insurance <- setNames(unname(reference), unname(mapping))
  expect_equal(digest::digest(reference, algo = "sha256"),
               "a6d8135c01cd4c1113e134b024b2b440a2e31c3e21e2e0a41b9b4d14d6c5aa1a")
  for (target in unname(mapping)) {
    rotated <- build_rhc_data_split(cohort, K = 4L, target_site = target)
    labels <- attr(rotated, "site_mapping")
    expect_identical(unname(labels["target"]), target)
    expect_equal(sum(vapply(rotated, `[[`, integer(1L), "n")), 5039L)
    expect_identical(rotated$t, by_insurance[[target]])
    for (index in seq_along(rotated)) {
      expect_identical(rotated[[index]], by_insurance[[labels[index]]])
      expect_identical(rotated[[index]]$W_outcome, rotated[[index]]$Z_site)
    }
    expect_identical(attr(rotated, "feature_center"), attr(reference, "feature_center"))
    expect_identical(attr(rotated, "feature_scale"), attr(reference, "feature_scale"))
  }
})

test_that("historical RHC defaults retain a single weight column", {
  fixture <- rhc_interface_fixture()
  result <- fixture$run(K = 3, n_folds = 3, comparison_methods = character(), verbose = FALSE)
  expect_identical(fixture$scope$forwarded$aggregation_mode, "common_tate")
  expect_identical(fixture$scope$forwarded$target_nuisance_method, "lasso")
  expect_equal(result$weights, c(s1 = .2, s2 = .3))
  expect_equal(result$pairwise$weight, c(.2, .3))
  expect_equal(result$metadata$target_anchor_weight, .5)
  expect_equal(as.character(result$methods$method), c("Target-only", "RoCE"))
})

test_that("joint RHC preserves arm weights and separates target references", {
  fixture <- rhc_interface_fixture()
  result <- fixture$run(K = 3, n_folds = 3, comparison_methods = character(), verbose = FALSE,
    aggregation_mode = "joint_tate", target_nuisance_method = "hou_calibrated",
    source_validation_method = "outer_fit", crossfit_layers = 2L,
    calibration_control = list(recipe = "score_derivative", target_radius = 5),
    nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
    nuisance_cv_certificate = TRUE, checkpoint_dir = "checkpoint-for-test", fold_seed = 19L)
  expect_equal(result$pairwise$weight_mu1, c(.2, .3))
  expect_equal(result$pairwise$weight_mu0, c(.1, -.1))
  expect_false("weight" %in% names(result$pairwise))
  expect_equal(result$pairwise$mean_wald_statistic_mu1, c(2, 5))
  expect_equal(result$pairwise$mean_wald_statistic_mu0, c(8, 11))
  expect_equal(result$metadata$target_anchor_weight_by_arm, c(mu1 = .5, mu0 = 1))
  expect_null(result$metadata$target_anchor_weight)
  expect_equal(result$target_only$estimate, .10)
  expect_equal(result$target_anchor$estimate, .15)
  expect_equal(as.character(result$methods$method), c("Target-only", "Calibrated target-only", "RoCE"))
  expect_identical(fixture$scope$forwarded$crossfit_layers, 2L)
  expect_identical(fixture$scope$forwarded$checkpoint_dir, "checkpoint-for-test")
  expect_identical(fixture$scope$forwarded$nuisance_cv_certificate, TRUE)
  expect_equal(fixture$scope$forwarded$nuisance_tol, 1e-10)
  expect_identical(fixture$scope$reference_folds, fixture$scope$forwarded$precomputed_folds$target_folds)
  expect_equal(result$metadata$K, 2)
  expect_equal(result$metadata$n_sites, 3)
})

test_that("RHC fold sensitivity changes folds without changing records", {
  fixture <- rhc_interface_fixture()
  call <- function(seed) fixture$run(K = 3, n_folds = 3, comparison_methods = character(),
                                     verbose = FALSE, fold_seed = seed)
  call(19L)
  first <- fixture$scope$forwarded
  call(20L)
  second <- fixture$scope$forwarded
  expect_identical(first$data_split, second$data_split)
  expect_false(identical(first$precomputed_folds, second$precomputed_folds))
  call(19L)
  expect_identical(first$precomputed_folds, fixture$scope$forwarded$precomputed_folds)
  expect_error(call(0), "positive integer")
  expect_error(call(2.5), "positive integer")
  expect_error(call(1e20), "positive integer")
})
