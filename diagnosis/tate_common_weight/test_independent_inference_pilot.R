#!/usr/bin/env Rscript

library(testthat)
source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")
source("scripts/slurm/atomic_output.R")
source("scripts/slurm/result_provenance.R")

pilot_fixture <- function() {
  estimate <- c(0, 0.2, 0.4, -0.2, 1, 2)
  analytic_se <- rep(0.1, 6L)
  truth <- rep(0.15, 6L)
  lower <- estimate - qnorm(0.975) * analytic_se
  upper <- estimate + qnorm(0.975) * analytic_se
  data.frame(
    sim_id = 10001L, rho = .inference_pilot_rhos(),
    method = "one_round_crossfit_ate", estimate = estimate, truth = truth,
    se = analytic_se, ci_lower = lower, ci_upper = upper,
    coverage = truth >= lower & truth <= upper,
    se_fixed_weight_bootstrap = 0.1,
    se_weight_relearn_bootstrap = 0.2,
    variance_weight_relearn_bootstrap = 0.04,
    weight_uncertainty_ratio = 2,
    weight_relearn_n_bootstrap = 1000L,
    weight_bootstrap_failures = 0L
  )
}

test_that("manifest fixes independent IDs and every scientific setting", {
  manifest <- roce_inference_pilot_manifest()
  expect_identical(manifest$task_id, 1:100)
  expect_identical(manifest$sim_id, 10001:10100)
  expect_length(intersect(manifest$sim_id, 1:500), 0L)
  expect_equal(manifest$n_weight_bootstrap, rep(1000L, 100L))
  expect_equal(manifest$n_bootstrap, rep(5000L, 100L))
  expect_equal(manifest$rho_values, rep("0;0.5;1;1.5;2;2.5", 100L))
  expect_equal(roce_validate_inference_pilot_manifest(manifest, "1")$sim_id,
               10001L)
  expect_equal(roce_validate_inference_pilot_manifest(manifest, 100)$sim_id,
               10100L)

  path <- tempfile(fileext = ".csv")
  on.exit(unlink(path), add = TRUE)
  write.csv(manifest, path, row.names = FALSE)
  roundtrip <- read.csv(path, stringsAsFactors = FALSE)
  expect_equal(roce_validate_inference_pilot_manifest(roundtrip, 1),
               manifest[1L, , drop = FALSE])

  for (field in names(manifest)) {
    broken <- manifest
    broken[[field]][1L] <- NA
    expect_error(roce_validate_inference_pilot_manifest(broken, 1))
  }
  for (field in c("p", "K", "cutoff", "n_site", "n_folds",
                  "nlambda_init", "n_bootstrap", "n_weight_bootstrap",
                  "M_tau", "M_tau_inference", "sim_id", "task_id")) {
    broken <- manifest
    broken[[field]][1L] <- broken[[field]][1L] + 1
    expect_error(roce_validate_inference_pilot_manifest(broken, 1))
  }
  expect_error(roce_validate_inference_pilot_manifest(manifest[-1L, ], 1))
  expect_error(roce_validate_inference_pilot_manifest(manifest[100:1, ], 1))
  expect_error(roce_validate_inference_pilot_manifest(manifest[, -1L], 1))
  expect_error(roce_validate_inference_pilot_manifest(
    cbind(manifest, extra = 1), 1))
  broken <- manifest
  broken$p <- matrix(broken$p, ncol = 1L)
  expect_error(roce_validate_inference_pilot_manifest(broken, 1))
  for (task in list(0, 101, 1.5, NA_real_, Inf, "x", c(1, 2),
                    TRUE, 1 + 1i, matrix(1, 1, 1), as.Date("2000-01-01"))) {
    expect_error(roce_validate_inference_pilot_manifest(manifest, task))
  }
})

test_that("diagnostic CI is separate and preserves all analytic inputs", {
  input <- pilot_fixture()
  before <- input
  observed <- roce_build_inference_audit(input)
  expect_identical(input, before)
  expect_equal(observed$analytic_se, input$se)
  expect_equal(observed$analytic_ci_lower, input$ci_lower)
  expect_equal(observed$analytic_ci_upper, input$ci_upper)
  expect_identical(observed$analytic_coverage, input$coverage)
  expect_equal(observed$weight_relearn_ci_lower,
               input$estimate - qnorm(0.975) * 0.2)
  expect_equal(observed$weight_relearn_ci_upper,
               input$estimate + qnorm(0.975) * 0.2)
  expect_identical(observed$weight_relearn_coverage,
    input$truth >= observed$weight_relearn_ci_lower &
      input$truth <= observed$weight_relearn_ci_upper)
  expect_true(all(observed$inference_status == "diagnostic_only"))
  expect_true(all(!observed$analytic_coverage | observed$weight_relearn_coverage))
  expect_equal(roce_build_inference_audit(input[6:1, ]), observed)

  other <- input
  other$method <- "target_only_ate"
  other$se_weight_relearn_bootstrap <- NA_real_
  expect_equal(roce_build_inference_audit(rbind(other, input)), observed)
  input$coverage <- as.integer(input$coverage)
  expect_equal(roce_build_inference_audit(input), observed)
})

test_that("invalid derived-CI inputs fail closed", {
  input <- pilot_fixture()
  expect_error(roce_build_inference_audit(input[-1L, ]))
  expect_error(roce_build_inference_audit(rbind(input, input[1L, ])))
  expect_error(roce_build_inference_audit(input[, -1L]))
  mutations <- list(
    function(x) { x$rho[1L] <- x$rho[2L]; x },
    function(x) { x$rho[1L] <- Inf; x },
    function(x) { x$sim_id[1L] <- 10002L; x },
    function(x) { x$sim_id[] <- 1L; x },
    function(x) { x$sim_id[] <- 10001.5; x },
    function(x) { x$sim_id[1L] <- NA_real_; x },
    function(x) { x$coverage[1L] <- NA; x },
    function(x) { x$coverage[1L] <- !x$coverage[1L]; x },
    function(x) { x$coverage <- rep(2, 6L); x },
    function(x) { x$coverage <- rep("yes", 6L); x },
    function(x) { x$ci_lower[1L] <- x$ci_lower[1L] - 0.1; x },
    function(x) { x$ci_upper[1L] <- x$ci_upper[1L] + 0.1; x },
    function(x) { x$se[1L] <- -1; x },
    function(x) { x$se_fixed_weight_bootstrap[1L] <- 0; x },
    function(x) { x$se_weight_relearn_bootstrap[1L] <- Inf; x },
    function(x) { x$variance_weight_relearn_bootstrap[1L] <- 0.01; x },
    function(x) { x$weight_uncertainty_ratio[1L] <- 1; x },
    function(x) { x$weight_relearn_n_bootstrap[1L] <- 500L; x },
    function(x) { x$weight_bootstrap_failures[1L] <- 1L; x }
  )
  for (mutate in mutations) {
    expect_error(roce_build_inference_audit(mutate(input)))
  }
})

test_that("scientific failure is retained without a success payload", {
  parent <- tempfile("inference_failure_test_")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  output <- file.path(parent, "failed_seed")
  task <- roce_inference_pilot_manifest()[1L, ]
  hashes <- list(package = paste(rep("a", 64L), collapse = ""),
                 workflow = paste(rep("b", 64L), collapse = ""),
                 manifest = paste(rep("c", 64L), collapse = ""))
  .inference_pilot_write_failure(output, task, "synthetic fit failure",
    "simulation", hashes, roce_sha256_file, roce_write_atomic_directory)
  expect_setequal(list.files(output),
    c("attempt_failure.csv", "metadata.txt", "sha256.txt"))
  failure <- read.csv(file.path(output, "attempt_failure.csv"))
  expect_identical(failure$status, "failed")
  expect_false(failure$retry_authorized)
  expect_identical(failure$failure_message, "synthetic fit failure")
  lines <- readLines(file.path(output, "sha256.txt"))
  expect_equal(vapply(file.path(output, substring(lines, 67L)),
                     roce_sha256_file, character(1L), USE.NAMES = FALSE),
               substr(lines, 1L, 64L))
  expect_error(.inference_pilot_write_failure(output, task, "overwrite",
    "simulation", hashes, roce_sha256_file, roce_write_atomic_directory))
  expect_identical(read.csv(file.path(output, "attempt_failure.csv")), failure)
})

test_that("scheduler mapping and output basename fail before fitting", {
  old <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("SLURM_ARRAY_TASK_ID") else
    Sys.setenv(SLURM_ARRAY_TASK_ID = old), add = TRUE)
  Sys.unsetenv("SLURM_ARRAY_TASK_ID")
  expect_true(.inference_pilot_validate_scheduler_task(1L))
  Sys.setenv(SLURM_ARRAY_TASK_ID = "1")
  expect_true(.inference_pilot_validate_scheduler_task(1L))
  expect_error(.inference_pilot_validate_scheduler_task(2L), "TASK_ID")
  for (value in c("0", "01", "1.5", "NaN", "-1", "abc")) {
    Sys.setenv(SLURM_ARRAY_TASK_ID = value)
    expect_error(.inference_pilot_validate_scheduler_task(1L), "TASK_ID")
  }
  Sys.unsetenv("SLURM_ARRAY_TASK_ID")
  parent <- tempfile("inference_path_test_")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  manifest <- file.path(parent, "manifest.csv")
  write.csv(roce_inference_pilot_manifest(), manifest, row.names = FALSE)
  expect_error(run_independent_inference_pilot_main(
    c(manifest, "1", file.path(parent, "wrong_seed"))), "seed_010001")
  expect_false(dir.exists(file.path(parent, "wrong_seed")))
})

test_that("publication failure remains a failure and never overwrites", {
  parent <- tempfile("inference_publication_test_")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE), add = TRUE)
  task <- roce_inference_pilot_manifest()[1L, ]
  hashes <- list(package = "test_package", workflow = "test_workflow",
                 manifest = "test_manifest")
  failure_output <- file.path(parent, "publication_failure")
  expect_error(roce_publish_inference_pilot(
    failure_output, task, hashes,
    function(directory) stop("synthetic publication failure"),
    roce_sha256_file, roce_write_atomic_directory),
    "synthetic publication failure")
  failed <- read.csv(file.path(failure_output, "attempt_failure.csv"))
  expect_identical(failed$status, "failed")
  expect_identical(failed$failure_stage, "publication")
  expect_false(failed$retry_authorized)
  expect_false(file.exists(file.path(failure_output, "results.csv")))

  success_output <- file.path(parent, "success")
  writer <- function(directory) writeLines("original",
    file.path(directory, "payload.txt"))
  expect_identical(roce_publish_inference_pilot(
    success_output, task, hashes, writer, roce_sha256_file,
    roce_write_atomic_directory), success_output)
  expect_error(roce_publish_inference_pilot(
    success_output, task, hashes, writer, roce_sha256_file,
    roce_write_atomic_directory))
  expect_identical(list.files(success_output), "payload.txt")
  expect_identical(readLines(file.path(success_output, "payload.txt")), "original")
  expect_false(file.exists(file.path(success_output, "attempt_failure.csv")))

  unavailable <- function(...) stop("synthetic unavailable publication")
  expect_error(roce_publish_inference_pilot(
    file.path(parent, "unavailable"), task, hashes, writer,
    roce_sha256_file, unavailable), "synthetic unavailable publication")
  expect_false(dir.exists(file.path(parent, "unavailable")))
})
