test_that("audit refuses an existing library before invoking R", {
  driver <- file.path(repo_root, "scripts/slurm/run_package_audit_tests.sh")
  skip_if_not(file.exists(driver), "Repository-only audit driver")
  root <- tempfile("audit_existing_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  marker <- file.path(root, "frozen.txt")
  writeLines("unchanged", marker)
  output <- suppressWarnings(system2("bash", shQuote(driver), stdout = TRUE, stderr = TRUE,
    env = c(paste0("ROCE_AUDIT_LIB=", shQuote(root)),
            paste0("ROCE_PROJECT_ROOT=", shQuote(repo_root)))))
  expect_identical(attr(output, "status"), 2L)
  expect_true(any(grepl("Refusing to reuse an existing audit library", output, fixed = TRUE)))
  expect_identical(list.files(root), "frozen.txt")
  expect_identical(readLines(marker), "unchanged")
})

test_that("ladder retains submitter errors and reports aggregate failure", {
  driver <- file.path(repo_root, "scripts/slurm/advance_production_ladder.sh")
  skip_if_not(file.exists(driver), "Repository-only ladder driver")
  root <- tempfile("ladder_test_")
  dir.create(file.path(root, "scripts/slurm"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  submitter <- file.path(root, "scripts/slurm/submit_rho_group_direct_tate.sh")
  writeLines(c("#!/bin/bash", 'echo "unfiltered error for ${ROCE_SETTING}" >&2', "exit 7"), submitter)
  Sys.chmod(submitter, "0755")
  environment <- paste0("ROCE_PROJECT_ROOT=", shQuote(root))
  output <- suppressWarnings(system2("bash", shQuote(driver), stdout = TRUE, stderr = TRUE,
                                    env = environment))
  expect_identical(attr(output, "status"), 1L)
  expect_equal(sum(grepl("unfiltered error", output, fixed = TRUE)), 12L)
  expect_equal(sum(grepl("exit=7", output, fixed = TRUE)), 12L)
  writeLines(c("#!/bin/bash", "echo already committed", "exit 0"), submitter)
  success <- system2("bash", shQuote(driver), stdout = TRUE, stderr = TRUE, env = environment)
  expect_null(attr(success, "status"))
  expect_equal(sum(grepl("already committed", success, fixed = TRUE)), 12L)
})

test_that("rule diagnostic rejects changed provenance and incomplete caches", {
  driver <- file.path(repo_root, "diagnosis/weight_rule_precision/weight_rule_precision.R")
  skip_if_not(file.exists(driver), "Repository-only diagnostic")
  old_directory <- setwd(repo_root)
  on.exit(setwd(old_directory), add = TRUE)
  diagnostic <- new.env(parent = globalenv())
  sys.source(driver, envir = diagnostic)
  settings <- list(config = "C1", package_fingerprint = strrep("a", 64),
                   workflow_fingerprint = strrep("b", 64))
  rows <- data.frame(sim_id = 1L, method = names(diagnostic$RULE_LABELS),
                     estimate = 0, truth = 0, se = 0.1, coverage = TRUE)
  for (field in names(settings)) rows[[field]] <- settings[[field]]
  expect_identical(diagnostic$.validate_rule_rows(rows, 1L, settings), rows)
  expect_error(diagnostic$.validate_rule_rows(rows[-1, ], 1L, settings), "incomplete")
  expect_error(diagnostic$.validate_rule_rows(rows, 2L, settings), "mismatched")
  changed <- settings
  changed$package_fingerprint <- strrep("c", 64)
  expect_error(diagnostic$.validate_rule_rows(rows, 1L, changed), "package_fingerprint")
  changed <- settings
  changed$workflow_fingerprint <- strrep("c", 64)
  expect_error(diagnostic$.validate_rule_rows(rows, 1L, changed), "workflow_fingerprint")
  rows$coverage[1] <- FALSE
  expect_error(diagnostic$.validate_rule_rows(rows, 1L, settings), "numerical")
})
