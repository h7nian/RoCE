simulation_qc_policy_path <- function() {
  testthat::test_path(
    "..", "..", "scripts", "slurm", "simulation_qc_policy.R"
  )
}

test_that("simulation QC separates recovered line searches from hard boundaries", {
  policy_path <- simulation_qc_policy_path()
  skip_if_not(
    file.exists(policy_path),
    "Repository-only Slurm QC policy is absent from the built package."
  )
  source(policy_path, local = TRUE)
  diagnostics <- data.frame(
    diagnostic_status = c(
      "nuisance_line_search_failures_detected",
      "nuisance_nonconvergence_detected",
      "ok"
    ),
    face_max_initial_dr_abs_coefficient_observed = c(2, 3, 90),
    face_max_calibrated_dr_abs_coefficient_observed = c(4, 5, 1),
    face_max_calibrated_outcome_abs_coefficient_observed = c(2, 95, 1),
    stringsAsFactors = FALSE
  )

  augmented <- roce_augment_simulation_qc(
    diagnostics, parameter_bound = 100
  )
  classification <- roce_classify_simulation_qc(augmented)

  expect_identical(
    augmented$nuisance_density_ratio_boundary_detected,
    c(FALSE, FALSE, TRUE)
  )
  expect_identical(
    augmented$nuisance_outcome_boundary_detected,
    c(FALSE, TRUE, FALSE)
  )
  expect_false(classification$implementation_failed[[1L]])
  expect_true(classification$statistical_review_required[[1L]])
  expect_true(classification$implementation_failed[[2L]])
  expect_true(classification$implementation_failed[[3L]])
  expect_match(
    augmented$diagnostic_status[[3L]],
    "nuisance_density_ratio_boundary_detected"
  )
  expect_match(
    augmented$diagnostic_status[[2L]],
    "nuisance_outcome_boundary_detected"
  )
  expect_identical(
    roce_augment_simulation_qc(
      augmented, parameter_bound = 100
    )$diagnostic_status,
    augmented$diagnostic_status
  )
})

test_that("simulation QC validates its schema and thresholds", {
  policy_path <- simulation_qc_policy_path()
  skip_if_not(
    file.exists(policy_path),
    "Repository-only Slurm QC policy is absent from the built package."
  )
  source(policy_path, local = TRUE)
  incomplete <- data.frame(diagnostic_status = "ok")

  expect_error(
    roce_augment_simulation_qc(incomplete, parameter_bound = 100),
    "missing QC columns"
  )
  expect_error(
    roce_augment_simulation_qc(
      transform(
        incomplete,
        face_max_initial_dr_abs_coefficient_observed = 1,
        face_max_calibrated_dr_abs_coefficient_observed = 1,
        face_max_calibrated_outcome_abs_coefficient_observed = 1
      ),
      parameter_bound = 100,
      boundary_fraction = 1
    ),
    "strictly between"
  )
})

test_that("simulation QC classifies design and density-ratio diagnostics", {
  policy_path <- simulation_qc_policy_path()
  skip_if_not(file.exists(policy_path), "Repository QC policy is absent.")
  source(policy_path, local = TRUE)
  diagnostics <- data.frame(
    diagnostic_status = c(
      "invalid_cell_diagnostics",
      "invalid_dr_weight_diagnostics",
      "sparse_binary_cells_detected",
      "density_ratio_clipping_detected"
    ),
    stringsAsFactors = FALSE
  )

  classification <- roce_classify_simulation_qc(diagnostics)

  expect_identical(
    classification$implementation_failed,
    c(TRUE, TRUE, FALSE, FALSE)
  )
  expect_identical(
    classification$statistical_review_required,
    c(FALSE, FALSE, TRUE, TRUE)
  )
})

test_that("simulation QC fails weight-relearning bootstrap violations", {
  policy_path <- simulation_qc_policy_path()
  skip_if_not(file.exists(policy_path), "Repository QC policy is absent.")
  source(policy_path, local = TRUE)
  diagnostics <- data.frame(
    diagnostic_status = c(
      "invalid_weight_relearn_bootstrap_diagnostics",
      "weight_relearn_bootstrap_failures_detected",
      paste(
        "coverage_below_mc_band",
        "invalid_weight_relearn_bootstrap_diagnostics",
        sep = ";"
      )
    ),
    stringsAsFactors = FALSE
  )

  classification <- roce_classify_simulation_qc(diagnostics)

  expect_identical(
    classification$implementation_failed,
    c(TRUE, TRUE, TRUE)
  )
  expect_identical(
    classification$statistical_review_required,
    c(FALSE, FALSE, TRUE)
  )
})
