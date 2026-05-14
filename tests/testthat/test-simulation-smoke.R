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
    dgp_type            = "facec"
  )

  expect_s3_class(res, "data.frame")
  expect_true("method" %in% names(res),
              info = "result must carry a 'method' column for summarization.")

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
