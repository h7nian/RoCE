#!/usr/bin/env Rscript
# Test the current RHC interface against the immutable, checked estimator library.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
library_path <- normalizePath(arguments[1L], mustWork = TRUE)
repository <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
interface <- new.env(parent = asNamespace("RoCE"))
sys.source(file.path(repository, "R/real_data_rhc.R"), envir = interface)
dir.create(output, recursive = TRUE)
test_path <- file.path(output, "test-rhc-current-interface.R")
stopifnot(file.copy(file.path(repository, "tests/testthat/test-rhc-current-interface.R"), test_path))
Sys.setenv(ROCE_TEST_INSTALLED = "1")
results <- testthat::test_file(test_path,
                              env = interface, reporter = "summary", stop_on_failure = TRUE)
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
saveRDS(results, file.path(output, "interface_tests.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
checked_files <- file.path(repository, c("R/real_data_rhc.R",
  "tests/testthat/test-rhc-current-interface.R", "diagnosis/repair/check_rhc_current.R"))
hashes <- vapply(checked_files, function(path) digest::digest(path, algo = "sha256", file = TRUE), character(1L))
jsonlite::write_json(as.list(hashes), file.path(output, "checked_files.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines("RHC_CURRENT_INTERFACE_CHECKS_PASSED", file.path(output, "INTERFACE_CHECKS_PASSED"))
