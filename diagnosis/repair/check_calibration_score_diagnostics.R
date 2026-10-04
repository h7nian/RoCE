#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=4L) stop("Usage: checked_R_library repository execution_checks new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
repository <- normalizePath(arguments[2L],mustWork=TRUE)
fixtures <- normalizePath(arguments[3L],mustWork=TRUE)
output <- arguments[4L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/repair/layer_comparison_summary.R"))
source(file.path(repository,"diagnosis/repair/calibration_score_diagnostics.R"))
library(testthat)
packets <- models <- list(); rows <- list()
for(variant in c("full","standard_source","ordinary_target")) {
  saved <- readRDS(file.path(fixtures,paste0(variant,".rds")))
  artifacts <- attr(saved,"roce_simulation_artifacts")
  fit <- artifacts$direct_tate_results$one_round_crossfit
  packet <- layer_inference_packet(fit,artifacts$data_split,artifacts$arm_truth)
  packet$data_hash <- digest::digest(artifacts$data_split,algo="sha256")
  packets[[variant]] <- packet
  model <- fit_packet_joint_weights(packet)
  models[[variant]] <- model
  evaluated <- evaluate_packet_weights(packet,model$weights,model)
  reference <- reaggregate_tate_crossfit(artifacts$data_split,fit,M_tau_inference=12,
    aggregation_mode="joint_tate",verbose=FALSE)
  detail <- summarize_layer_arms(reference,artifacts$data_split)
  test_that(paste(variant,"packet reproduces the independently assembled estimator"), {
    expect_equal(evaluated$estimates[["tate"]],reference$estimate,tolerance=1e-12)
    expect_equal(evaluated$tate_variance,reference$variance,tolerance=1e-12)
    expect_equal(evaluated$tate_variance_fixed_weights,reference$variance_fixed_weights,tolerance=1e-12)
    expect_equal(unname(evaluated$covariance),unname(detail$covariance),tolerance=1e-12)
    for(arm in c("mu1","mu0")) expect_equal(unname(model$weights[[arm]]),
      unname(reference$fold_weights_by_arm[[arm]]),tolerance=1e-12)
    candidates <- summarize_packet_candidates(packet)
    tate <- candidates[candidates$quantity=="tate",]
    expect_equal(tate$estimate[1L],fit$target_only$estimate,tolerance=1e-12)
    expect_equal(tate$estimate[-1L],as.numeric(fit$source_estimates),tolerance=1e-12)
    expect_true(all(candidates$score_variance>0))
  })
  rows[[variant]] <- cbind(variant=variant,summarize_packet_candidates(packet))
}
test_that("shared-weight comparisons reject different data and evaluation folds", {
  expect_silent(validate_paired_score_packets(packets$full,packets$standard_source))
  changed <- packets$standard_source
  changed$data_hash <- paste(rep("0",64),collapse="")
  expect_error(evaluate_packet_weights(changed,models$full$weights,models$full),"identical data")
  changed <- packets$standard_source
  for(arm in c("mu1","mu0")) changed$arms[[arm]]$fold_info[c(1L,2L)] <- changed$arms[[arm]]$fold_info[c(2L,1L)]
  expect_error(evaluate_packet_weights(changed,models$full$weights,models$full),"identical outer folds")
  changed <- packets$full
  changed$arms$mu1$fold_info[[1L]]$target_idx[1L] <- changed$arms$mu1$fold_info[[1L]]$target_idx[2L]
  expect_error(fit_packet_joint_weights(changed),"exactly once")
  expect_error(fit_packet_joint_weights(packets$full,0),"positive")
})
transferred <- evaluate_packet_weights(packets$standard_source,models$full$weights,models$full)
own <- evaluate_packet_weights(packets$standard_source,models$standard_source$weights,models$standard_source)
test_that("cutoff and weight-transfer diagnostics preserve the fixed-weight component", {
  fixed <- evaluate_packet_weights(packets$standard_source,models$full$weights)
  expect_identical(fixed$estimates,transferred$estimates)
  expect_identical(fixed$covariance_fixed_weights,transferred$covariance_fixed_weights)
  expect_gt(max(abs(own$estimates-transferred$estimates)),1e-10)
  adjusted <- fit_packet_joint_weights(packets$full,1.5)
  expect_equal(adjusted$fits[[1L]]$lambda,1/1.5)
  expect_true(all(is.finite(evaluate_packet_weights(packets$full,adjusted$weights,adjusted)$covariance)))
})
dir.create(output,recursive=TRUE)
write.csv(do.call(rbind,rows),file.path(output,"candidate_scores.csv"),row.names=FALSE)
saveRDS(list(shared_eta=transferred,own_eta=own,
  source_hashes=tools::md5sum(file.path(repository,"diagnosis/repair",
    c("calibration_score_diagnostics.R","check_calibration_score_diagnostics.R")))),file.path(output,"checks.rds"))
writeLines(c("CALIBRATION_SCORE_DIAGNOSTICS_CHECKS_PASSED",
  "K8 rho1.5: candidate means, full joint point/variance reconstruction, data/fold guards and shared-eta calculations checked.",
  "This validates the diagnostic interface, not a finite-sample calibration advantage or a full nuisance-refit variance."),
  file.path(output,"checks.txt"))
