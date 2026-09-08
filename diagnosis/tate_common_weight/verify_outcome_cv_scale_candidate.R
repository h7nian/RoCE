#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=2L)stop("usage: verify_outcome_cv_scale_candidate.R LIBRARY OUTPUT")
  lib<-normalizePath(args[1],mustWork=TRUE)
  .libPaths(c(lib,.libPaths()));Sys.setenv(ROCE_TEST_INSTALLED="1",NOT_CRAN="true")
  suppressPackageStartupMessages(library(RoCE,lib.loc=lib))
  stopifnot(normalizePath(find.package("RoCE"))==file.path(lib,"RoCE"))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  fingerprint<-.weight_calibration_installed_package_fingerprint(lib,roce_sha256_file)
  files<-list.files("tests/testthat",pattern="[.]R$",full.names=TRUE)
  hashes<-vapply(files,roce_sha256_file,"")
  result<-testthat::test_dir("tests/testthat",reporter="summary",stop_on_failure=FALSE)
  summary<-as.data.frame(result)
  targeted<-testthat::test_file("diagnosis/tate_common_weight/test_outcome_cv_training_scale.R",reporter="summary")
  targeted_summary<-as.data.frame(targeted)
  stopifnot(identical(hashes,vapply(files,roce_sha256_file,"")),
    fingerprint==.weight_calibration_installed_package_fingerprint(lib,roce_sha256_file))
  roce_write_atomic_directory(args[2],function(stage) {
    saveRDS(list(repository=result,targeted=targeted),file.path(stage,"test_results.rds"))
    write.csv(summary[,!vapply(summary,is.list,logical(1L)),drop=FALSE],file.path(stage,"repository_test_summary.csv"),row.names=FALSE)
    write.csv(data.frame(path=files,sha256=hashes),file.path(stage,"test_hashes.csv"),row.names=FALSE)
    writeLines(c(paste0("package_fingerprint=",fingerprint),paste0("passed_assertions=",sum(summary$passed)),
      paste0("failed_assertions=",sum(summary$failed)),paste0("error_tests=",sum(summary$error)),
      paste0("skipped_tests=",sum(summary$skipped)),paste0("targeted_passed=",sum(targeted_summary$passed)),
      "NOT_CRAN=true","inference_validated=FALSE"),file.path(stage,"metadata.txt"))
    payload<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(payload,roce_sha256_file,""),basename(payload),sep="  "),file.path(stage,"sha256.txt"))
  },caller="installed outcome CV candidate tests")
  stopifnot(sum(summary$failed)==0,!any(summary$error),sum(targeted_summary$failed)==0,!any(targeted_summary$error))
}
if(sys.nframe()==0L)main()
