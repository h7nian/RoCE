test_that("both target anchors validate the complete standard source program", {
  skip_on_cran()
  set.seed(4605)
  data <- split_data_by_site(generate_bounded_data(n_target = 600, n_source_sizes = c(600, 600),
    K = 2, p = 4, config = "C2", n_deviated_sites = 0))
  run <- function(target_method, data_split = data) {
    control <- list(recipe = "score_derivative", source_nuisance_method = "standard")
    if (target_method == "hou_calibrated") {
      control$target_propensity_initialization <- "calibrated"
      control$target_radius <- 12
    }
    run_tate_crossfit(data_split, n_folds = 4, communication_mode = "one_round", nlambda_init = 6,
      target_nuisance_method = target_method, source_validation_method = "complete",
      nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, calibration_layout = "compact",
      calibration_control = control, M_tau = 12, M_tau_inference = 12,
      lambda_selection = .5, aggregation_mode = "joint_tate", verbose = FALSE)
  }
  calibrated_target <- run("hou_calibrated")
  ordinary_target <- run("lasso")
  expect_equal(calibrated_target$source_estimates, ordinary_target$source_estimates, tolerance = 1e-12)
  expect_equal(calibrated_target$crossfit_levels, 3L)
  expect_equal(ordinary_target$crossfit_levels, 2L)
  for (fit in list(calibrated_target, ordinary_target)) {
    expect_true(is.finite(fit$estimate) && is.finite(fit$se) && fit$se > 0)
    expect_identical(fit$source_validation_method, "complete")
    for (arm in fit$arm_results) for (outer in seq_along(arm$fold_results)) {
      for (source in arm$fold_results[[outer]]$source_results) {
        expect_identical(source$source_nuisance_method, "standard")
        expect_equal(source$n_calibrated_folds, 0L)
        expect_false(source$nuisance_fit_diagnostics$final_score_calibration_applied)
        expect_null(source$per_k2_gamma)
        for (key in names(source$inner_fits)) {
          validation <- as.integer(sub("k2_", "", key))
          inner <- source$inner_fits[[key]]
          expect_setequal(inner$training_folds, setdiff(1:4, c(outer, validation)))
          expect_identical(inner$source_nuisance_method, "standard")
        }
      }
    }
  }
  changed <- data
  treated <- which(changed$s1$A == 1L)
  changed$s1$Y[treated[1:20]] <- 1 - changed$s1$Y[treated[1:20]]
  independent <- run("hou_calibrated", changed)
  reused <- .reuse_one_round_tate_across_rho(data, changed, calibrated_target,
    changed_sources = "s1", lambda_selection = .5, n_cores = 1L)
  expect_equal(reused$estimate, independent$estimate, tolerance = 1e-12)
  expect_equal(reused$se, independent$se, tolerance = 1e-12)
  expect_equal(reused$fold_weights, independent$fold_weights, tolerance = 1e-12)
  expect_equal(reused$source_estimates, independent$source_estimates, tolerance = 1e-12)
})
