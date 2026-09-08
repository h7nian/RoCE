#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(testthat))
source("scripts/slurm/result_provenance.R")
source("diagnosis/tate_common_weight/audit_independent_inference_pilot.R")
source("diagnosis/tate_common_weight/run_independent_inference_pilot.R")

write_manifest <- function(directory, files) {
  hashes <- vapply(file.path(directory, files), roce_sha256_file, character(1L))
  writeLines(paste(hashes, files, sep = "  "), file.path(directory, "sha256.txt"))
}

make_success <- function() {
  d <- tempfile("iia_success_"); dir.create(d)
  files <- c("results.csv", "inference_audit.csv", "diagnostic_qc.csv",
             "artifacts.rds", "metadata.txt")
  for (name in files) writeLines(paste("payload", name), file.path(d, name))
  write_manifest(d, files); d
}

test_that("exact successful payload is accepted and tampering fails closed", {
  d <- make_success()
  expect_identical(.iia_classify_bundle(d, roce_sha256_file)$status, "complete")
  writeLines("tampered", file.path(d, "results.csv"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "checksum mismatch")
})

test_that("malformed, extra, duplicate, and invalid hashes fail closed", {
  d <- make_success(); writeLines("stray",file.path(d,"extra.txt"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "extra or missing")
  d <- make_success(); lines <- readLines(file.path(d, "sha256.txt"))
  writeLines(c(lines, lines[[1L]]), file.path(d, "sha256.txt"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "wrong exact payload")
  d <- make_success(); lines <- readLines(file.path(d, "sha256.txt"))
  writeLines(c(lines[-1L], sub("metadata.txt$", "extra.txt", lines[[length(lines)]])),
             file.path(d, "sha256.txt"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "wrong exact payload")
  d <- make_success(); lines <- readLines(file.path(d, "sha256.txt"))
  substr(lines[[1L]], 1L, 1L) <- "z"; writeLines(lines, file.path(d, "sha256.txt"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "wrong exact payload")
})

test_that("failed attempt preserves explicit reason and is not success", {
  d <- tempfile("iia_failed_"); dir.create(d)
  write.csv(data.frame(status="failed",failure_stage="scientific_or_validation",
    failure_message="deliberate fixture failure"),file.path(d,"attempt_failure.csv"),row.names=FALSE)
  writeLines("independent_inference_pilot=failed",file.path(d,"metadata.txt"))
  write_manifest(d,c("attempt_failure.csv","metadata.txt"))
  x <- .iia_classify_bundle(d, roce_sha256_file)
  expect_identical(x$status,"failed")
  expect_identical(x$stage,"scientific_or_validation")
  expect_identical(x$reason,"deliberate fixture failure")
  write.csv(data.frame(status="failed",failure_stage="",failure_message="why"),
            file.path(d,"attempt_failure.csv"),row.names=FALSE)
  write_manifest(d,c("attempt_failure.csv","metadata.txt"))
  expect_error(.iia_classify_bundle(d, roce_sha256_file), "explicit stage/reason")
})

test_that("metadata requires parseable unique keys", {
  p <- tempfile(); writeLines(c("a=1","b=two"),p)
  expect_identical(.iia_metadata(p),c(a="1",b="two"))
  writeLines(c("a=1","a=2"),p); expect_error(.iia_metadata(p),"duplicate")
  writeLines(c("a=1","missing"),p); expect_error(.iia_metadata(p),"malformed")
})

test_that("ancestor manifest discovery is exact and absence fails", {
  root <- tempfile("iia_manifest_"); dir.create(root); dir.create(file.path(root,"a","b"),recursive=TRUE)
  write.csv(data.frame(x=1),file.path(root,"manifest.csv"),row.names=FALSE)
  expect_identical(.iia_manifest_path(file.path(root,"a","b")),normalizePath(file.path(root,"manifest.csv")))
  unlink(file.path(root,"manifest.csv")); expect_error(.iia_manifest_path(file.path(root,"a","b")),"no ancestor")
})

test_that("canonical manifest validates task and sim IDs", {
  manifest <- roce_inference_pilot_manifest()
  row <- roce_validate_inference_pilot_manifest(manifest,1L)
  expect_equal(row$task_id,1L); expect_equal(row$sim_id,10001L)
  bad <- manifest; bad$sim_id[[1L]] <- 99L
  expect_error(roce_validate_inference_pilot_manifest(bad,1L),"mismatch")
  expect_error(roce_validate_inference_pilot_manifest(manifest,101L),"TASK_ID")
})

test_that("strict integer parsing rejects fractional metadata IDs", {
  expect_identical(.iia_strict_integer("1", "task", 1L, 100L),1L)
  expect_error(.iia_strict_integer("1.5","task",1L,100L),"strict integer")
  expect_error(.iia_strict_integer(10001.5,"sim",10001L,10100L),"strict integer")
  for (bad in list(TRUE, 1 + 1i, matrix(1, 1, 1), factor("1"))) {
    expect_error(.iia_strict_integer(bad, "id"), "strict integer")
  }
})

test_that("every saved task field agrees without truncating IDs", {
  expected <- roce_inference_pilot_manifest()[1L, , drop = FALSE]
  expect_invisible(.iia_validate_saved_task(expected, expected))
  for (name in names(expected)) {
    bad <- expected
    bad[[name]][1L] <- NA
    expect_error(.iia_validate_saved_task(bad, expected), "saved task")
  }
  for (name in c("task_id", "sim_id", "n_site", "n_weight_bootstrap")) {
    bad <- expected
    bad[[name]][1L] <- bad[[name]][1L] + 0.5
    expect_error(.iia_validate_saved_task(bad, expected), "saved task")
  }
  bad <- expected
  bad$methods <- "target_only"
  expect_error(.iia_validate_saved_task(bad, expected), "saved task")
  expect_error(.iia_validate_saved_task(expected[, -1L], expected), "saved task")
})

test_that("hard nuisance fields reject unknown and invalid counts", {
  fields <- c("initial_outcome_degenerate","target_only_outcome_degenerate",
    "initial_dr_nonconverged","calibrated_dr_nonconverged",
    "calibrated_outcome_nonconverged","initial_dr_line_search_failures",
    "calibrated_dr_line_search_failures","calibrated_outcome_line_search_failures",
    "initial_dr_support_floor_applied","calibrated_dr_support_floor_applied")
  fit <- list(nuisance_fit_diagnostics=as.data.frame(setNames(rep(list(c(0,1)),length(fields)),fields)))
  expect_equal(nrow(.iia_nuisance(fit,0)),length(fields))
  for (bad in list(NA_real_,Inf,-1,0.5)) {
    x <- fit; x$nuisance_fit_diagnostics[[fields[[1L]]]][[1L]] <- bad
    expect_error(.iia_nuisance(x,0),"known finite nonnegative integers")
  }
})

test_that("numerical errors have fixed finite nonnegative schema", {
  names <- c("a","b"); expect_equal(.iia_numeric_errors(c(a=0,b=1),names),c(a=0,b=1))
  expect_error(.iia_numeric_errors(c(a=0),names),"fixed")
  expect_error(.iia_numeric_errors(c(a=0,b=NA),names),"finite")
  expect_error(.iia_numeric_errors(c(a=0,b=Inf),names),"finite")
  expect_error(.iia_numeric_errors(c(a=0,b=-1),names),"nonnegative")
})

test_that("data and fit scalar schema cannot silently omit variance or influence", {
  site <- list(n=1000L,X=matrix(0,1000,100),A=rep(0L,1000),Y=rep(0L,1000))
  entry <- list(data_split=list(t=site,s1=site,s2=site))
  fit <- list(variance=1,all_phi_agg=rep(0,3000),N_all=3000,
              n_folds=5L,n_sites=3L)
  expect_invisible(.iia_validate_data_fit(entry,fit))
  bad <- fit; bad$variance <- NULL
  expect_error(.iia_validate_data_fit(entry,bad),"variance/influence")
  bad <- fit; bad$all_phi_agg[[1L]] <- NA_real_
  expect_error(.iia_validate_data_fit(entry,bad),"variance/influence")
  badentry <- entry; badentry$data_split$s1$X <- matrix(0,999,100)
  expect_error(.iia_validate_data_fit(badentry,fit),"site schema")
  for (field in c("A", "Y")) {
    for (value in c(Inf, 0.5, 2, NA_real_)) {
      badentry <- entry
      badentry$data_split$s1[[field]][1L] <- value
      expect_error(.iia_validate_data_fit(badentry,fit), "site schema")
    }
  }
  bad <- fit; bad$N_all <- NULL
  expect_error(.iia_validate_data_fit(entry,bad), "strict integer")
  bad <- fit; bad$all_phi_agg <- matrix(bad$all_phi_agg, ncol = 1L)
  expect_error(.iia_validate_data_fit(entry,bad), "variance/influence")
})

test_that("existing audit output is rejected before package work", {
  output <- tempfile("iia_existing_"); dir.create(output)
  expect_error(audit_independent_inference_pilot_main(c(tempdir(),tempdir(),output)),
               "already exists")
})

cat("audit_independent_inference_pilot targeted tests: PASS\n")
