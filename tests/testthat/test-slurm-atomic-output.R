test_that("Slurm atomic-directory helper commits success and cleans failure", {
  helper_path <- file.path(
    repo_root, "scripts", "slurm", "atomic_output.R"
  )
  skip_if_not(
    file.exists(helper_path),
    "Repository-only Slurm helper is absent from the built package."
  )
  helper_environment <- new.env(parent = baseenv())
  sys.source(
    helper_path,
    envir = helper_environment
  )
  write_atomic <- helper_environment$roce_write_atomic_directory

  test_root <- tempfile("roce_atomic_output_")
  dir.create(test_root)
  on.exit(unlink(test_root, recursive = TRUE, force = TRUE), add = TRUE)

  successful_output <- file.path(test_root, "successful")
  expect_invisible(write_atomic(
    successful_output,
    writer = function(staging_directory) {
      writeLines("complete", file.path(staging_directory, "value.txt"))
    },
    caller = "unit test"
  ))
  expect_identical(
    readLines(file.path(successful_output, "value.txt")), "complete"
  )

  failed_output <- file.path(test_root, "failed")
  expect_error(
    write_atomic(
      failed_output,
      writer = function(staging_directory) {
        writeLines("partial", file.path(staging_directory, "partial.txt"))
        stop("intentional writer failure")
      },
      caller = "unit test"
    ),
    "intentional writer failure"
  )
  expect_false(file.exists(failed_output))
  expect_setequal(
    list.files(test_root, all.files = TRUE, no.. = TRUE),
    "successful"
  )

  expect_error(
    write_atomic(
      successful_output,
      writer = function(staging_directory) {
        writeLines("overwrite", file.path(staging_directory, "value.txt"))
      },
      caller = "unit test"
    ),
    "output already exists"
  )

})

test_that("TATE task helpers share and validate sensitivity defaults", {
  helper_path <- file.path(
    repo_root, "scripts", "slurm", "direct_tate_task_helpers.R"
  )
  skip_if_not(
    file.exists(helper_path),
    "Repository-only Slurm helper is absent from the built package."
  )
  helper_environment <- new.env(parent = baseenv())
  sys.source(helper_path, envir = helper_environment)

  bundle_root <- tempfile("roce_task_helper_")
  on.exit(unlink(bundle_root, recursive = TRUE, force = TRUE), add = TRUE)
  lock_path <- file.path(bundle_root, "locks", "task.lock")
  release_lock <- helper_environment$roce_claim_task_lock(
    lock_path, c("task=1", "job=local")
  )
  expect_identical(
    readLines(file.path(lock_path, "claim.txt")),
    c("task=1", "job=local")
  )
  expect_error(
    helper_environment$roce_claim_task_lock(lock_path, "task=2"),
    "already claimed"
  )
  expect_invisible(release_lock())
  expect_false(file.exists(lock_path))

  expect_true(helper_environment$roce_parse_boolean("TRUE", "flag"))
  expect_false(helper_environment$roce_parse_boolean("0", "flag"))
  expect_error(
    helper_environment$roce_parse_boolean("yes", "flag"),
    "flag must be TRUE/FALSE or 1/0",
    fixed = TRUE
  )
  scheduler <- helper_environment$roce_scheduler_provenance()
  expect_setequal(
    names(scheduler),
    c(
      "slurm_job_id", "slurm_array_job_id", "slurm_array_task_id",
      "slurm_cluster_name", "slurm_job_partition", "slurm_node_list",
      "compute_hostname"
    )
  )
  expect_true(all(vapply(scheduler, function(value) {
    is.character(value) && length(value) == 1L && nzchar(value)
  }, logical(1L))))
  scheduler_frame <- as.data.frame(scheduler, stringsAsFactors = FALSE)
  expect_invisible(
    helper_environment$roce_validate_scheduler_provenance(scheduler_frame)
  )
  invalid_scheduler <- scheduler_frame
  invalid_scheduler$slurm_job_id <- ""
  expect_error(
    helper_environment$roce_validate_scheduler_provenance(invalid_scheduler),
    "nonmissing and nonempty"
  )

  sensitivity_grid <- helper_environment$roce_default_sensitivity_grid()
  reference_cutoff <- 1 / RoCE:::AGG_WALD_LAMBDA
  expect_equal(nrow(sensitivity_grid), 8L)
  expect_equal(
    sensitivity_grid$cutoff,
    c(1, 1.5, 2, 2.5, 3, rep(reference_cutoff, 3L))
  )
  expect_equal(sensitivity_grid$M_tau_inference, c(5, 5, 5, 5, 5, 4, 6, Inf))

  task <- data.frame(
    task_id = 11L, sim_id = 7L, experiment = "negative_transfer",
    config = "C3", p = 100L, K = 4L, rho = 2.5,
    cutoff = reference_cutoff, n_site = 1000L, n_folds = 5L,
    M_tau = 5, M_tau_inference = 5, deviation_mechanism = "treated_arm"
  )
  expect_true(helper_environment$roce_is_reused_sensitivity_task(task))
  expect_false(helper_environment$roce_is_reused_sensitivity_task(
    transform(task, K = 2L)
  ))
  # Sidecar sensitivities exist only for the treated-arm family.
  expect_false(helper_environment$roce_is_reused_sensitivity_task(
    transform(task, deviation_mechanism = "both_arms")
  ))
  expect_identical(helper_environment$roce_task_deviation_mechanism(task), "treated_arm")
  expect_error(
    helper_environment$roce_task_deviation_mechanism(transform(task, deviation_mechanism = "shared")),
    "deviation_mechanism"
  )

  resource_plan <- list(
    allocated_cores = 40L, nuisance_cv_threads = 5L,
    source_workers = 4L, parallel_treatment_arms = TRUE,
    fully_parallel_cores = 40L
  )
  provenance <- list(
    library = "/tmp/roce-audited-library",
    fingerprint = paste(rep("a", 64L), collapse = "")
  )
  workflow_fingerprint <- paste(rep("b", 64L), collapse = "")
  manifest_fingerprint <- paste(rep("c", 64L), collapse = "")
  primary_row <- data.frame(
    sim_id = 7L, p = 100L, K = 4L, config = "C3",
    method = "one_round_crossfit_ate", deviation_mechanism = "treated_arm"
  )
  annotated_primary <- helper_environment$roce_annotate_direct_tate_rows(
    primary_row, task, resource_plan, 12.5, provenance,
    workflow_fingerprint, manifest_fingerprint,
    metadata_policy = "from_task"
  )
  expect_identical(annotated_primary$experiment, "negative_transfer")
  expect_equal(annotated_primary$cutoff, reference_cutoff)

  sidecar_row <- transform(
    primary_row,
    task_id = 11L,
    experiment = "c3_reused_sensitivity",
    primary_experiment = "negative_transfer",
    rho = 2.5,
    cutoff = 1.5,
    aggregation_cutoff = 1.5,
    aggregation_lambda = 1 / 1.5,
    primary_cutoff = reference_cutoff,
    n_site = 1000L,
    n_folds = 5L
  )
  annotated_sidecar <- helper_environment$roce_annotate_direct_tate_rows(
    sidecar_row, task, resource_plan, 12.5, provenance,
    workflow_fingerprint, manifest_fingerprint,
    metadata_policy = "preserve"
  )
  expect_identical(annotated_sidecar$experiment, "c3_reused_sensitivity")
  expect_equal(annotated_sidecar$cutoff, 1.5)
  expect_equal(annotated_sidecar$aggregation_cutoff, 1.5)
  expect_equal(annotated_sidecar$primary_cutoff, reference_cutoff)
  expect_error(
    helper_environment$roce_annotate_direct_tate_rows(
      transform(sidecar_row, experiment = NA_character_), task,
      resource_plan, 12.5, provenance, workflow_fingerprint,
      manifest_fingerprint, metadata_policy = "preserve"
    ),
    "internally inconsistent"
  )
  expect_error(
    helper_environment$roce_annotate_direct_tate_rows(
      transform(primary_row, p = 99L), task, resource_plan, 12.5,
      provenance, workflow_fingerprint, manifest_fingerprint,
      metadata_policy = "from_task"
    ),
    "does not match task"
  )

  primary <- data.frame(
    method = c("target_only_ate", "one_round_crossfit_ate"),
    estimate = c(0.1, 0.2),
    se = c(0.03, 0.02),
    bias = c(0.01, 0.02),
    coverage = c(1, 1),
    ci_width = c(0.12, 0.08)
  )
  sensitivity <- transform(
    primary,
    sensitivity_kind = "reference_identity"
  )
  expect_invisible(
    helper_environment$roce_validate_sensitivity_identity(
      sensitivity, primary
    )
  )
  sensitivity$estimate[[1L]] <- sensitivity$estimate[[1L]] + 1e-4
  expect_error(
    helper_environment$roce_validate_sensitivity_identity(
      sensitivity, primary
    ),
    "did not reproduce"
  )

  csv_bundle_root <- tempfile("roce_csv_bundle_")
  on.exit(unlink(csv_bundle_root, recursive = TRUE, force = TRUE), add = TRUE)
  output_paths <- file.path(
    csv_bundle_root, c("sidecar/a.csv", "primary/b.csv")
  )
  sentinel_path <- file.path(csv_bundle_root, "commits/group.txt")
  expect_invisible(helper_environment$roce_commit_csv_bundle(
    values = list(data.frame(value = 1), data.frame(value = 2)),
    output_paths = output_paths,
    sentinel_path = sentinel_path,
    sentinel_lines = "commit=complete"
  ))
  expect_equal(read.csv(output_paths[[1L]])$value, 1)
  expect_equal(read.csv(output_paths[[2L]])$value, 2)
  expect_identical(readLines(sentinel_path), "commit=complete")
  expect_error(
    helper_environment$roce_commit_csv_bundle(
      values = list(data.frame(value = 3)),
      output_paths = output_paths[[1L]]
    ),
    "refusing to overwrite"
  )
  expect_error(
    helper_environment$roce_commit_csv_bundle(
      values = list(1),
      output_paths = file.path(csv_bundle_root, "not_a_frame.csv")
    ),
    "data frame"
  )
  expect_error(
    helper_environment$roce_commit_csv_bundle(
      values = list(data.frame(value = 4)),
      output_paths = file.path(csv_bundle_root, "duplicate.csv"),
      sentinel_path = file.path(csv_bundle_root, "duplicate.csv"),
      sentinel_lines = "commit=complete"
    ),
    "distinct"
  )
})

test_that("reused sensitivity rows preserve their scientific setting", {
  sensitivity_path <- file.path(
    repo_root, "scripts", "slurm", "direct_tate_sensitivity_rows.R"
  )
  skip_if_not(
    file.exists(sensitivity_path),
    "Repository-only sensitivity helper is absent from the built package."
  )
  sensitivity_environment <- new.env(parent = asNamespace("RoCE"))
  sys.source(sensitivity_path, envir = sensitivity_environment)

  task <- data.frame(
    sim_id = 17L,
    task_id = 23L,
    experiment = "negative_transfer",
    rho = 1.5,
    K = 2L,
    p = 100L,
    config = "C1",
    n_site = 1000L,
    n_folds = 5L,
    deviation_mechanism = "treated_arm"
  )
  row <- sensitivity_environment$roce_make_tate_result_row(
    task = task,
    fit = list(estimate = 0.19, se = 0.05, M_tau = 5),
    method = "one_round_crossfit_ate",
    truth = 0.1,
    cutoff = 1.5,
    M_tau_inference = 5,
    primary_cutoff = 1,
    n_total = 3000L,
    nlambda_init = 100L,
    n_bootstrap = 5000L,
    outcome_family = "binomial"
  )

  expect_equal(row$bias, 0.09)
  expect_true(row$coverage)
  expect_identical(row$heterogeneity_type, "one_deviated_source")
  expect_identical(row$outcome_family, "binomial")
  expect_identical(row$experiment, "c3_reused_sensitivity")
  expect_identical(row$primary_experiment, "negative_transfer")
  # Every column roce_annotate_direct_tate_rows() asserts against the task must
  # be present, or the sidecar fails only after the fit has already run.
  expect_identical(row$deviation_mechanism, "treated_arm")
  for (field in c("sim_id", "p", "K", "config", "deviation_mechanism")) {
    expect_true(field %in% names(row))
  }
  expect_equal(row$cutoff, 1.5)
  expect_equal(row$aggregation_cutoff, 1.5)
  expect_equal(row$primary_cutoff, 1)
})

test_that("sensitivity aggregation follows the recorded primary cutoff", {
  aggregate_path <- file.path(
    repo_root, "scripts", "slurm", "aggregate_reused_sensitivities.R"
  )
  skip_if_not(file.exists(aggregate_path), "Sensitivity aggregator is absent.")
  aggregate_driver <- readLines(aggregate_path, warn = FALSE)

  expect_true(any(grepl(
    '"M_tau_inference", "primary_cutoff"',
    aggregate_driver, fixed = TRUE
  )))
  expect_true(any(grepl(
    "package_primary_cutoff", aggregate_driver, fixed = TRUE
  )))
  expect_false(any(grepl(
    "sidecar$cutoff == 2", aggregate_driver, fixed = TRUE
  )))
  expect_false(any(grepl(
    "fitting$cutoff == 2", aggregate_driver, fixed = TRUE
  )))
})

test_that("RHC artifact driver scopes cleanup inside main", {
  driver_path <- file.path(
    repo_root, "scripts", "slurm", "run_rhc_direct_tate.R"
  )
  skip_if_not(file.exists(driver_path), "Repository RHC driver is absent.")
  driver <- readLines(driver_path, warn = FALSE)

  expect_true(any(grepl("^main <- function\\(\\) \\{$", driver)))
  expect_true(any(grepl("^on[.]exit\\(\\{$", trimws(driver))))
  expect_identical(tail(driver[nzchar(trimws(driver))], 1L), "main()")
  expect_true(any(grepl("committed_files", driver, fixed = TRUE)))
  expect_true(any(grepl("commit_complete", driver, fixed = TRUE)))
  expect_true(any(grepl("AGG_WALD_LAMBDA", driver, fixed = TRUE)))
  expect_false(any(grepl(
    "aggregation_lambda = 0.5", driver, fixed = TRUE
  )))
  audit_path <- file.path(
    repo_root, "scripts", "slurm", "audit_rhc_direct_tate.R"
  )
  audit <- readLines(audit_path, warn = FALSE)
  expect_true(any(grepl(
    "comparison_density_ratio_clipping_diagnostics", audit, fixed = TRUE
  )))
  expect_true(any(grepl(
    "expected_primary_cutoff", audit, fixed = TRUE
  )))
  expect_true(any(grepl("^main <- function", audit)))
  expect_identical(
    tail(audit[nzchar(trimws(audit))], 1L),
    "main()"
  )
  expect_true(any(grepl(
    "roce_write_atomic_directory", audit, fixed = TRUE
  )))
  expect_true(any(grepl(
    'file.path(output_root, "audits", audit_label)', audit, fixed = TRUE
  )))
  expect_true(any(grepl(
    "rhc_direct_tate_audit_failed.txt", audit, fixed = TRUE
  )))
})

test_that("seed-level task drivers scope locks and commit atomically", {
  driver_paths <- file.path(
    repo_root, "scripts", "slurm",
    c(
      "run_direct_tate_task.R", "run_direct_tate_cutoff_task.R",
      "run_grouped_cutoff_diagnostic_task.R"
    )
  )
  skip_if_not(all(file.exists(driver_paths)), "Repository task drivers absent.")

  for (driver_path in driver_paths) {
    driver <- readLines(driver_path, warn = FALSE)
    expect_true(any(grepl("^main <- function", driver)))
    expect_identical(
      tail(driver[nzchar(trimws(driver))], 1L),
      "main()"
    )
    expect_true(any(grepl("roce_claim_task_lock", driver, fixed = TRUE)))
    expect_true(any(grepl("on.exit(release_lock()", driver, fixed = TRUE)))
    expect_true(any(grepl("roce_commit_csv_bundle", driver, fixed = TRUE)))
  }
})

test_that("bootstrap benchmark preserves provenance and refuses overwrite", {
  benchmark_r_path <- file.path(
    repo_root, "scripts", "slurm", "benchmark_bootstrap_replicates.R"
  )
  benchmark_sh_path <- file.path(
    repo_root, "scripts", "slurm", "benchmark_bootstrap_replicates.sh"
  )
  skip_if_not(
    file.exists(benchmark_r_path) && file.exists(benchmark_sh_path),
    "Repository bootstrap benchmark drivers are absent."
  )
  benchmark_r <- readLines(benchmark_r_path, warn = FALSE)
  benchmark_sh <- readLines(benchmark_sh_path, warn = FALSE)

  expect_true(any(grepl(
    "refusing to overwrite existing bootstrap benchmark",
    benchmark_r, fixed = TRUE
  )))
  expect_true(any(grepl(
    "roce_commit_csv_bundle", benchmark_r, fixed = TRUE
  )))
  expect_true(any(grepl(
    "package_fingerprint", benchmark_r, fixed = TRUE
  )))
  expect_true(any(grepl(
    "workflow_fingerprint", benchmark_r, fixed = TRUE
  )))
  expect_true(any(grepl(
    "ROCE_BOOTSTRAP_WORKFLOW_FINGERPRINT", benchmark_sh, fixed = TRUE
  )))
})

test_that("C3 remainder diagnostic is scoped and workflow-provenanced", {
  diagnostic_r_path <- file.path(
    repo_root, "scripts", "slurm", "diagnose_c3_target_remainder.R"
  )
  diagnostic_sh_path <- file.path(
    repo_root, "scripts", "slurm", "diagnose_c3_target_remainder.sh"
  )
  skip_if_not(
    file.exists(diagnostic_r_path) && file.exists(diagnostic_sh_path),
    "Repository C3 diagnostic drivers are absent."
  )
  diagnostic_r <- readLines(diagnostic_r_path, warn = FALSE)
  diagnostic_sh <- readLines(diagnostic_sh_path, warn = FALSE)

  expect_true(any(grepl("^main <- function", diagnostic_r)))
  expect_identical(tail(diagnostic_r[nzchar(trimws(diagnostic_r))], 1L), "main()")
  expect_true(any(grepl("outcome_family", diagnostic_r, fixed = TRUE)))
  expect_true(any(grepl("workflow_fingerprint", diagnostic_r, fixed = TRUE)))
  expect_true(any(grepl(
    "c3_target_remainder_contrast.csv", diagnostic_r, fixed = TRUE
  )))
  expect_true(any(grepl(
    "ROCE_C3_DIAGNOSTIC_WORKFLOW_FINGERPRINT",
    diagnostic_sh, fixed = TRUE
  )))
  comparison_path <- file.path(
    repo_root, "scripts", "slurm", "compare_c3_nuisance_rules.R"
  )
  expect_true(file.exists(comparison_path))
  comparison <- readLines(comparison_path, warn = FALSE)
  expect_true(any(grepl("^main <- function", comparison)))
  expect_identical(
    tail(comparison[nzchar(trimws(comparison))], 1L), "main()"
  )
})

test_that("cutoff selection is disjoint and does not tune on coverage", {
  selector_path <- file.path(
    repo_root, "scripts", "slurm", "select_grouped_cutoff.R"
  )
  summary_path <- file.path(
    repo_root, "scripts", "slurm", "summarize_grouped_cutoff_diagnostic.R"
  )
  skip_if_not(
    file.exists(selector_path) && file.exists(summary_path),
    "Repository cutoff selector or summarizer is absent."
  )
  selector <- readLines(selector_path, warn = FALSE)
  summarizer <- readLines(summary_path, warn = FALSE)

  expect_true(any(grepl("simulation_ids.*501--510", selector)))
  expect_true(any(grepl("zero$rmse <= rho_zero_tolerance", selector,
                        fixed = TRUE)))
  expect_true(any(grepl("worst_nontransport_rmse", selector, fixed = TRUE)))
  expect_true(any(grepl("coverage_not_used", selector, fixed = TRUE)))
  expect_false(any(grepl("which.max.*coverage", selector)))
  expect_true(any(grepl(
    "cutoff_summary_sha256", summarizer, fixed = TRUE
  )))
  expect_true(any(grepl(
    "observed_summary_fingerprint", selector, fixed = TRUE
  )))
  expect_true(any(grepl(
    "cutoff_decision_fingerprint", selector, fixed = TRUE
  )))
})

test_that("manifest builders preserve existing experiment definitions", {
  builder_paths <- file.path(
    repo_root, "scripts", "slurm",
    c(
      "build_direct_tate_manifest.R",
      "build_direct_tate_cutoff_manifest.R",
      "build_direct_tate_rho_group_manifest.R",
      "split_direct_tate_manifest_by_k.R",
      "build_grouped_cutoff_diagnostic_manifest.R"
    )
  )
  skip_if_not(all(file.exists(builder_paths)), "Manifest builders are absent.")
  builders <- lapply(builder_paths, readLines, warn = FALSE)

  expect_true(all(vapply(builders, function(builder) {
    any(grepl("refusing to overwrite existing", builder, fixed = TRUE)) ||
      any(grepl("refusing to overwrite existing K-specific", builder,
                fixed = TRUE))
  }, logical(1L))))
  expect_true(all(vapply(builders, function(builder) {
    any(grepl("file.rename", builder, fixed = TRUE))
  }, logical(1L))))
  expect_true(any(grepl(
    "ROCE_PRIMARY_CUTOFF", builders[[1L]], fixed = TRUE
  )))
  expect_true(any(grepl(
    "ROCE_PRIMARY_CUTOFF", builders[[5L]], fixed = TRUE
  )))
  expect_false(any(grepl(
    "cutoffs = 2", builders[[1L]], fixed = TRUE
  )))
})

test_that("production submission is cutoff-gated and batch-bounded", {
  driver_paths <- file.path(
    repo_root, "scripts", "slurm",
    c(
      "prepare_mc500_manifests.sh",
      "submit_main_direct_tate.sh",
      "submit_rho_group_direct_tate.sh",
      "run_rho_group_checkpoint_audit.sh",
      "submit_direct_tate.sh",
      "submit_cutoff_diagnostic.sh"
    )
  )
  skip_if_not(all(file.exists(driver_paths)), "Production drivers are absent.")
  drivers <- lapply(driver_paths, readLines, warn = FALSE)

  expect_true(all(vapply(drivers[seq_len(4L)], function(driver) {
    any(grepl("PRIMARY_CUTOFF", driver, fixed = TRUE))
  }, logical(1L))))
  expect_true(any(grepl(
    "CUTOFF_SELECTION_FINGERPRINT", drivers[[1L]], fixed = TRUE
  )))
  expect_true(any(grepl(
    "CUTOFF_SELECTION_GATE", drivers[[2L]], fixed = TRUE
  )))
  expect_true(any(grepl(
    "CUTOFF_SELECTION_GATE", drivers[[3L]], fixed = TRUE
  )))
  # HISTORY #0010 raised the grouped caps so the 18 setting blocks can run
  # concurrently, each scoped by ROCE_SETTING and bounded by its own ladder.
  expect_true(any(grepl(
    '"${BATCH_SIZE}" -gt 500', drivers[[3L]], fixed = TRUE
  )))
  expect_true(any(grepl(
    '"${MAX_CONCURRENT}" -gt 50', drivers[[3L]], fixed = TRUE
  )))
  expect_true(any(grepl("ROCE_SETTING", drivers[[3L]], fixed = TRUE)))
  expect_true(any(grepl(
    "enforce_roce_array_safety_cap 5 2", drivers[[5L]], fixed = TRUE
  )))
  expect_true(any(grepl(
    "enforce_roce_array_safety_cap 5 2", drivers[[6L]], fixed = TRUE
  )))
})

test_that("installed-package provenance covers both lazy-load artifacts", {
  helper_path <- file.path(
    repo_root, "scripts", "slurm", "package_library_utils.sh"
  )
  skip_if_not(file.exists(helper_path), "Package-library helper is absent.")
  helper <- readLines(helper_path, warn = FALSE)
  expect_true(any(grepl("R/RoCE.rdb", helper, fixed = TRUE)))
  expect_true(any(grepl("R/RoCE.rdx", helper, fixed = TRUE)))
})

test_that("production resource policy is queue-aware and internally exact", {
  resource_path <- file.path(
    repo_root, "scripts", "slurm", "resource_topology.R"
  )
  skip_if_not(file.exists(resource_path), "Repository resource helper absent.")
  resource_environment <- new.env(parent = baseenv())
  sys.source(resource_path, envir = resource_environment)

  expect_identical(
    resource_environment$roce_expected_production_cv_threads(c(2L, 4L, 8L)),
    c(5L, 5L, 2L)
  )
  rows <- data.frame(
    K = c(2L, 4L, 8L),
    n_folds = 5L,
    allocated_cores = c(20L, 40L, 32L),
    nuisance_cv_threads = c(5L, 5L, 2L),
    source_workers_per_arm = c(2L, 4L, 8L),
    parallel_treatment_arms = TRUE,
    fully_parallel_cores = c(20L, 40L, 32L)
  )
  expect_invisible(
    resource_environment$roce_validate_production_resource_metadata(rows)
  )
  rows$nuisance_cv_threads[[3L]] <- 5L
  expect_error(
    resource_environment$roce_validate_production_resource_metadata(rows),
    "resource metadata|production resource metadata"
  )
})
