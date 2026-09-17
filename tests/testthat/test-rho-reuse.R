rho_equivalence_columns <- function(result) {
  columns <- unique(c(
    "method", "estimate", "se", "bias", "coverage", "ci_width",
    "truth", "ci_lower", "ci_upper",
    grep(
      paste0(
        "^(face_|direct_tate_|target_only_tate_|",
        "comparison_|or_degenerate_folds$)"
      ),
      names(result), value = TRUE
    )
  ))
  setdiff(columns, grep("seconds$", columns, value = TRUE))
}

test_that("both rho reuse paths retain the fitted aggregation rule", {
  observed <- character()
  testthat::local_mocked_bindings(
    .refit_one_round_crossfit_sources = function(fitted_arm, ...) fitted_arm,
    calculate_tate_crossfit_aggregation = function(screening_rule = "soft_penalty", ...) {
      observed <<- c(observed, screening_rule)
      list(aggregation_screening_rule = screening_rule)
    },
    .package = "RoCE"
  )
  arm <- list(nuisance_fit_diagnostics = data.frame(fold = 1L), timing = list())
  fitted <- list(communication_mode = "one_round", arm_results = list(mu1 = arm, mu0 = arm))
  for (rule in c("soft_penalty", "hard_threshold", "quadratic_bias")) {
    fitted$aggregation_screening_rule <- rule
    for (both_arms in c(FALSE, TRUE)) {
      reused <- RoCE:::.reuse_one_round_tate_across_rho(
        list(), list(), fitted, "s1", 1, refit_control_arm = both_arms
      )
      expect_identical(reused$aggregation_screening_rule, rule)
    }
  }
  fitted$aggregation_screening_rule <- NULL
  legacy <- RoCE:::.reuse_one_round_tate_across_rho(list(), list(), fitted, "s1", 1)
  expect_identical(legacy$aggregation_screening_rule, "soft_penalty")
  expect_identical(observed, c(rep(c("soft_penalty", "hard_threshold", "quadratic_bias"), each = 2), "soft_penalty"))
})

test_that("mixed integer and fractional rho artifacts retain exact keys", {
  make_artifact <- function(rho) {
    list(data_split = list(marker = rho), direct_tate_results = list(
      one_round_crossfit = list(estimate = rho)
    ))
  }
  testthat::local_mocked_bindings(
    run_single_simulation = function(...) {
      result <- data.frame(estimate = 0)
      attr(result, "roce_simulation_artifacts") <- make_artifact(0)
      result
    },
    .run_face_positive_rho_update = function(
        rho, simulation_args, reuse_reference, keep_artifact = FALSE) {
      list(result = data.frame(estimate = rho, rho = rho),
           artifacts = if (keep_artifact) make_artifact(rho) else NULL)
    },
    .package = "RoCE"
  )
  rhos <- c(0, 0.5, 1, 1.5, 2, 2.5)
  args <- list(methods = "one_round_crossfit", estimate_ate = TRUE,
               dgp_type = "face")
  for (requested in list(rhos, c(2.5, 0, 1), numeric(0))) {
    grouped <- RoCE:::.run_face_rho_group(
      args, rhos, changed_sources = "s1", artifact_rhos = requested
    )
    expect_identical(names(grouped$results),
                     c("0", "0.5", "1", "1.5", "2", "2.5"))
    expect_length(grouped$artifacts, length(requested))
    expect_identical(names(grouped$artifacts),
                     vapply(requested, format, character(1L),
                            scientific = FALSE, trim = TRUE))
    expect_false(any(vapply(grouped$artifacts, is.null, logical(1L))))
    for (rho in requested) {
      artifact <- grouped$artifacts[[format(rho, trim = TRUE)]]
      expect_equal(artifact$data_split$marker, rho)
      expect_equal(artifact$direct_tate_results$one_round_crossfit$estimate,
                   rho)
    }
  }
  expect_error(RoCE:::.run_face_rho_group(
    args, rhos, changed_sources = "s1", artifact_rhos = c(0, 0)
  ), "unique subset")
})

test_that("same-seed FACE rho reuse reproduces an independent TATE fit", {
  skip_on_cran()

  common_args <- list(
    sim_id = 731L,
    n_total = 360L,
    K = 1L,
    p = 3L,
    config = "C3",
    methods = c("one_round_crossfit", "target_only"),
    verbose = FALSE,
    n_cores_internal = 2L,
    nlambda_init = 10L,
    n_folds = 3L,
    n_bootstrap = 20L,
    estimate_ate = TRUE,
    dgp_type = "face"
  )
  grouped <- suppressWarnings(RoCE:::.run_face_rho_group(
    simulation_args = common_args,
    rho_values = c(0, 0.5),
    changed_sources = "s1",
    artifact_rhos = 0.5
  ))
  reference <- grouped$results[[1L]]
  independent <- suppressWarnings(
    do.call(
      run_single_simulation,
      c(common_args, list(ate_deviation = 0.5, n_deviated_sites = 1L))
    )
  )
  reused <- grouped$results[[2L]]

  comparison_columns <- rho_equivalence_columns(independent)
  expect_equal(
    reused[comparison_columns],
    independent[comparison_columns],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_true(all(reused$rho_reuse_enabled))
  expect_true(all(!reference$rho_reuse_enabled))

  expect_identical(grouped$reuse$positive_rho_workers, 1L)
  expect_identical(grouped$reuse$positive_rho_backend, "serial")
  expect_identical(
    grouped$reuse$positive_rho_source_workers_per_fit, 1L
  )
  reused_fit <- grouped$artifacts[["0.5"]]$
    direct_tate_results$one_round_crossfit
  expect_true(isTRUE(reused_fit$rho_reuse$reused))
  expect_true(isTRUE(reused_fit$rho_reuse$reused_control_arm))
  expect_equal(reused_fit$rho_reuse$changed_sources, "s1")
})

test_that("PSOCK rho workers exactly reproduce serial grouped updates", {
  skip_on_cran()
  skip_if_not(
    identical(Sys.getenv("ROCE_TEST_INSTALLED"), "1"),
    "PSOCK workers require the current package to be installed."
  )

  common_args <- list(
    sim_id = 733L,
    n_total = 180L,
    K = 1L,
    p = 4L,
    config = "C3",
    methods = c("one_round_crossfit", "target_only"),
    verbose = FALSE,
    n_cores_internal = 1L,
    nlambda_init = 4L,
    n_folds = 3L,
    n_bootstrap = 2L,
    estimate_ate = TRUE,
    dgp_type = "face"
  )
  run_group <- function(workers) {
    suppressWarnings(RoCE:::.run_face_rho_group(
      simulation_args = common_args,
      rho_values = c(0, 0.25, 0.5),
      changed_sources = "s1",
      artifact_rhos = 0.5,
      positive_rho_workers = workers
    ))
  }
  serial <- run_group(1L)
  psock <- run_group(2L)
  for (index in seq_along(serial$results)) {
    columns <- rho_equivalence_columns(serial$results[[index]])
    expect_equal(
      psock$results[[index]][columns],
      serial$results[[index]][columns],
      tolerance = 1e-12,
      ignore_attr = TRUE
    )
  }
  expect_identical(psock$reuse$positive_rho_workers, 2L)
  expect_identical(psock$reuse$positive_rho_backend, "psock")
  expect_true(is.list(psock$artifacts[["0.5"]]))
})

test_that("positive-rho worker count is validated and capped", {
  base_args <- list(
    simulation_args = list(
      methods = "one_round_crossfit", dgp_type = "face",
      estimate_ate = TRUE
    ),
    rho_values = c(0, 0.5)
  )
  expect_error(
    do.call(RoCE:::.run_face_rho_group, c(
      base_args, list(positive_rho_workers = 0)
    )),
    "positive_rho_workers must be a positive integer",
    fixed = TRUE
  )
  expect_error(
    do.call(RoCE:::.run_face_rho_group, c(
      base_args, list(positive_rho_workers = 1.5)
    )),
    "positive_rho_workers must be a positive integer",
    fixed = TRUE
  )
})

test_that("rho reuse rejects changes outside the declared treated source", {
  skip_on_cran()

  set.seed(911L)
  reference_data <- generate_simulation_data(
    n_total = 360L, K = 1L, p = 3L, config = "C3",
    outcome_type = "binary", dgp_type = "face",
    ate_deviation = 0, n_deviated_sites = 0L,
    warn_ignored = FALSE
  )
  reference_split <- split_data_by_site(reference_data)
  current_split <- reference_split
  current_split$t$Y[[1L]] <- 1 - current_split$t$Y[[1L]]

  expect_error(
    RoCE:::.validate_one_round_rho_reuse_data(
      reference_split, current_split, "s1", 1L
    ),
    "outcome changed at non-refitted site 't'",
    fixed = TRUE
  )

  current_split <- reference_split
  control_index <- which(current_split$s1$A == 0L)[[1L]]
  current_split$s1$Y[[control_index]] <-
    1 - current_split$s1$Y[[control_index]]
  expect_error(
    RoCE:::.validate_one_round_rho_reuse_data(
      reference_split, current_split, "s1", 1L
    ),
    "changed outcomes outside the refitted treatment arm",
    fixed = TRUE
  )
  # Under the shared-shift deviation both arms of the changed source may differ.
  expect_identical(
    RoCE:::.validate_one_round_rho_reuse_data(
      reference_split, current_split, "s1", 1L, changed_arms = c(0L, 1L)
    ),
    "s1"
  )
  expect_error(
    RoCE:::.validate_one_round_rho_reuse_data(
      reference_split, current_split, "s1", 1L, changed_arms = 0L
    ),
    "containing the refitted arm",
    fixed = TRUE
  )
})

test_that("same-seed FACE shared-shift reuse refits both arms and reproduces an independent fit", {
  skip_on_cran()

  # K = 2 so that the informative source s2 is reused in both arms while s1
  # is refitted in both arms.
  common_args <- list(
    sim_id = 733L,
    n_total = 360L,
    K = 2L,
    p = 3L,
    config = "C1",
    methods = c("one_round_crossfit", "target_only"),
    verbose = FALSE,
    n_cores_internal = 1L,
    nlambda_init = 10L,
    n_folds = 3L,
    n_bootstrap = 20L,
    estimate_ate = TRUE,
    dgp_type = "face",
    deviation_mechanism = "both_arms"
  )
  grouped <- suppressWarnings(RoCE:::.run_face_rho_group(
    simulation_args = common_args,
    rho_values = c(0, 1),
    changed_sources = "s1",
    artifact_rhos = 1
  ))
  independent <- suppressWarnings(
    do.call(
      run_single_simulation,
      c(common_args, list(ate_deviation = 1, n_deviated_sites = 1L))
    )
  )
  reused <- grouped$results[[2L]]
  comparison_columns <- rho_equivalence_columns(independent)
  expect_equal(
    reused[comparison_columns],
    independent[comparison_columns],
    tolerance = 1e-12,
    ignore_attr = TRUE
  )
  expect_true(all(reused$heterogeneity_type == "one_shared_shift_source"))
  expect_true(all(reused$deviation_mechanism == "both_arms"))
  reused_fit <- grouped$artifacts[["1"]]$direct_tate_results$one_round_crossfit
  expect_false(reused_fit$rho_reuse$reused_control_arm)
  expect_identical(reused_fit$rho_reuse$changed_sources, "s1")
  expect_identical(reused_fit$arm_results$mu0$rho_reuse$reused_sources, "s2")
  # The rho = 0 reference is the same dataset under either mechanism.
  reference <- grouped$results[[1L]]
  treated_reference <- suppressWarnings(do.call(
    run_single_simulation,
    c(common_args[names(common_args) != "deviation_mechanism"],
      list(ate_deviation = 0, n_deviated_sites = 0L))
  ))
  expect_equal(
    reference[comparison_columns], treated_reference[comparison_columns],
    tolerance = 1e-12, ignore_attr = TRUE
  )
})
