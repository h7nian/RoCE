# test-simulation-smoke.R - End-to-end smoke test for run_single_simulation
#
# This test sits one level ABOVE test-cross_fitting.R: it does not only verify
# that run_crossfit("one_round"/"two_round") returns valid objects, but that
# the orchestrator run_single_simulation correctly assembles result-rows for
# both modes into the data.frame that run_simulation_study later writes to CSV.
#
# This is the direct regression coverage for the bug where the simulation
# summary CSV was missing rows because one of the source-estimates matrices
# never reached the aggregation step. If run_single_simulation ever silently
# drops a requested method again, this test fails.
#
# Cost: one run_single_simulation invocation with a minimal method list (no
# comparison estimators, no oracle) takes a few tens of seconds. Skipped when
# the C++ backend is not compiled.

library(testthat)

# Package loaded by helper-load.R.

test_that("simulation artifacts require TATE estimation", {
  expect_error(
    run_single_simulation(
      sim_id = 1L,
      estimate_ate = FALSE,
      return_fitted_tate = TRUE,
      verbose = FALSE
    ),
    "requires estimate_ate = TRUE"
  )
})

test_that("simulation entry points validate opt-in weight bootstrap", {
  expect_error(
    run_single_simulation(
      sim_id = 1L, n_weight_bootstrap = 1L, verbose = FALSE
    ),
    "0 \\(disabled\\) or at least 2"
  )
  expect_error(
    run_simulation_study(n_sims = 1L, n_weight_bootstrap = -1L),
    "non-negative integer"
  )
  expect_error(
    run_single_simulation(
      sim_id = 1L, estimate_ate = FALSE,
      n_weight_bootstrap = 2L, verbose = FALSE
    ),
    "requires estimate_ate = TRUE"
  )
})

test_that("run_single_simulation assembles rows for target_only + one_round + two_round", {
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")

  methods <- c("target_only", "one_round_crossfit", "two_round_crossfit")

  # Smallest config the pipeline supports while still keeping each fold above
  # MIN_TREATED_FOR_MODEL: n_total = 600, K = 2 -> 200 per site, 67 per fold.
  # NOTE: p must be >= 4 (transform_covariates indexes X[, 4] unconditionally
  # for transform_type "mild"/"strong").
  res <- run_single_simulation(
    sim_id              = 1L,
    n_total             = 600L,
    K                   = 2L,
    p                   = 4L,
    config              = "C1",
    methods             = methods,
    verbose             = FALSE,
    nlambda_init        = 3L,
    outcome_type        = "binary",
    heterogeneity_type  = "none",
    n_folds             = 3L,
    estimate_ate        = FALSE,
    dgp_type            = "roce"
  )

  expect_s3_class(res, "data.frame")
  expect_true("method" %in% names(res),
              info = "result must carry a 'method' column for summarization.")
  expect_true("nlambda_init" %in% names(res))
  expect_true(all(res$nlambda_init == 3L))
  expect_true("nuisance_lambda_rule" %in% names(res))
  expect_true(all(res$nuisance_lambda_rule == "min"))
  expect_true("n_weight_bootstrap" %in% names(res))
  expect_true(all(res$n_weight_bootstrap == 0L))
  weight_bootstrap_fields <- c(
    "variance_weight_relearn_bootstrap", "se_weight_relearn_bootstrap",
    "se_fixed_weight_bootstrap", "weight_uncertainty_ratio",
    "weight_relearn_n_bootstrap", "weight_bootstrap_seed",
    "weight_bootstrap_failures", "weight_bootstrap_multiplier",
    "weight_bootstrap_relearn_weights",
    "weight_bootstrap_screening_rule"
  )
  expect_true(all(weight_bootstrap_fields %in% names(res)))
  expect_true(all(vapply(
    res[weight_bootstrap_fields], function(value) all(is.na(value)),
    logical(1L)
  )))

  missing_methods <- setdiff(methods, unique(as.character(res$method)))
  expect_equal(length(missing_methods), 0L,
               info = sprintf("run_single_simulation dropped requested methods: %s",
                              paste(missing_methods, collapse = ", ")))

  # Every requested method must have at least one row, and core numeric columns
  # must be finite.
  for (m in methods) {
    rows <- res[as.character(res$method) == m, , drop = FALSE]
    expect_gt(nrow(rows), 0L)
    for (col in c("estimate", "se", "bias", "ci_width")) {
      if (col %in% names(rows)) {
        expect_true(all(is.finite(rows[[col]])),
                    info = sprintf("[%s] column '%s' has non-finite values", m, col))
      }
    }
    if ("coverage" %in% names(rows)) {
      expect_true(all(rows$coverage %in% c(TRUE, FALSE, NA)),
                  info = sprintf("[%s] coverage must be logical", m))
    }
  }
})

test_that("simulation entry points reject an unknown nuisance CV rule", {
  expect_error(
    run_single_simulation(
      sim_id = 1L,
      nuisance_lambda_rule = "unsupported",
      verbose = FALSE
    ),
    "nuisance_lambda_rule must be either 'min' or '1se'",
    fixed = TRUE
  )
  expect_error(
    run_simulation_study(
      n_sims = 1L,
      nuisance_lambda_rule = "unsupported"
    ),
    "nuisance_lambda_rule must be either 'min' or '1se'",
    fixed = TRUE
  )
})

test_that("simulation treatment-arm scheduling does not change TATE results", {
  skip_on_cran()
  skip_on_os("windows")
  skip_if_not(exists("fit_general_glm_cpp"), message = "C++ not compiled")

  common_args <- list(
    sim_id = 917L,
    n_total = 360L,
    K = 1L,
    p = 4L,
    config = "C1",
    methods = c("one_round_crossfit", "target_only"),
    verbose = FALSE,
    n_cores_internal = 1L,
    nlambda_init = 3L,
    outcome_type = "continuous",
    n_folds = 3L,
    estimate_ate = TRUE,
    include_hard_threshold_diagnostic = TRUE,
    dgp_type = "face",
    n_target = 120L,
    n_source_sizes = 240L
  )
  sequential <- do.call(
    run_single_simulation,
    c(common_args, list(parallel_treatment_arms = FALSE))
  )
  concurrent <- do.call(
    run_single_simulation,
    c(common_args, list(parallel_treatment_arms = TRUE))
  )

  scientific_columns <- c(
    "method", "estimate", "se", "bias", "coverage", "ci_width", "truth"
  )
  expect_identical(
    sequential[scientific_columns],
    concurrent[scientific_columns]
  )
  hard <- sequential[
    sequential$method == "one_round_crossfit_ate_hard_threshold",
    , drop = FALSE
  ]
  expect_equal(nrow(hard), 1L)
  quadratic <- sequential[
    sequential$method == "one_round_crossfit_ate_quadratic_bias",
    , drop = FALSE
  ]
  expect_equal(nrow(quadratic), 1L)
  expect_true(all(sequential$misspecification_strength == 0))
  expect_true(all(is.finite(unlist(quadratic[c(
    "estimate", "se", "se_fixed_weights", "weight_layer_indirect_variance"
  )]))))
  expect_identical(quadratic$weight_layer_kink_cells, 0)
  expect_true(all(sequential$quadratic_bias_rule_requested))
  expect_true(all(sequential$deviation_mechanism == "treated_arm"))
  expect_true(all(is.finite(unlist(hard[c(
    "face_initial_dr_nonconverged",
    "face_calibrated_dr_nonconverged",
    "face_calibrated_outcome_nonconverged"
  )]))))
})

test_that("FACE heterogeneity labels follow the generated deviation", {
  expect_identical(RoCE:::.face_heterogeneity_type(0, 0L), "none")
  expect_identical(RoCE:::.face_heterogeneity_type(0, 2L), "none")
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 1L),
    "one_deviated_source"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 2L),
    "multiple_deviated_sources"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(0, 0L, 0.5),
    "source_effect_modification"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 1L, 0.5),
    "deviation_and_effect_modification"
  )
  expect_error(RoCE:::.face_heterogeneity_type(-1, 1L), "non-negative")
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 1L, deviation_mechanism = "both_arms"),
    "one_shared_shift_source"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 2L, deviation_mechanism = "both_arms"),
    "multiple_shared_shift_sources"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(1.5, 1L, 0.5, deviation_mechanism = "both_arms"),
    "shared_shift_and_effect_modification"
  )
  expect_identical(
    RoCE:::.face_heterogeneity_type(0, 1L, deviation_mechanism = "both_arms"),
    "none"
  )
})
test_that("hard-threshold diagnostic flag is strictly logical", {
  expect_error(
    run_single_simulation(
      sim_id = 1L,
      include_hard_threshold_diagnostic = "yes"
    ),
    "must be TRUE or FALSE",
    fixed = TRUE
  )
})

test_that("quadratic-bias rule flag is strictly logical", {
  expect_error(
    run_single_simulation(sim_id = 1L, include_quadratic_bias_rule = NA),
    "include_quadratic_bias_rule must be TRUE or FALSE",
    fixed = TRUE
  )
})
