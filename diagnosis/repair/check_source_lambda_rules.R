#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/"))
library_path <- file.path(root, "Rlib")
.libPaths(c(library_path, .libPaths()))
Sys.setenv(ROCE_TEST_INSTALLED = "1", NOT_CRAN = "true")
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
stopifnot(identical(normalizePath(find.package("RoCE")), normalizePath(file.path(library_path, "RoCE"))))
selected <- paste(c("source-lambda-rules", "calibration-training", "nested-calibration",
  "crossfit-layers", "nuisance-training", "nuisance-cv-certificate", "source-calibration-ablation",
  "source-nuisance-program", "joint-tate", "layer-arm-variance", "score-calibration"), collapse = "|")
results <- testthat::test_dir(file.path(root, "source/tests/testthat"), filter = selected,
  reporter = "summary", stop_on_failure = TRUE)
saveRDS(results, file.path(root, "test_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(root, "session_info.txt"))
for (marker in c("TESTS_PASSED", "CROSSFIT_LAYER_CHECKS_PASSED", "CV_CERTIFICATE_CHECKS_PASSED",
                 "SOURCE_STANDARD_CHECKS_PASSED", "SOURCE_LAMBDA_RULES_CHECKS_PASSED")) {
  writeLines(c(marker, "Selected installed-package regressions and source-stage isolation checks passed."),
             file.path(root, marker))
}
