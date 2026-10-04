#!/usr/bin/env Rscript
# Full-grid K8 execution checks, separate from the prescribed MC seeds.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: checked_R_library repository_root new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
repository <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/repair/layer_comparison_summary.R"))
dir.create(output,recursive=TRUE)
base_arguments <- list(sim_id=701L,n_target=1000L,n_source_sizes=rep(1000L,8L),K=8L,p=10L,
  config="C2",n_folds=10L,nlambda_init=100L,dgp_type="bounded",outcome_type="binary",
  estimand_type="superpopulation",ate_deviation=1.5,n_deviated_sites=1L,deviation_mechanism="treated_arm",
  methods=c("one_round_crossfit","target_only"),estimate_ate=TRUE,return_fitted_tate=TRUE,
  include_quadratic_bias_rule=FALSE,n_bootstrap=5000L,n_cores_internal=1L,parallel_treatment_arms=TRUE,
  nuisance_solver="proximal_newton",nuisance_tol=1e-10,nuisance_cv_certificate=TRUE,
  crossfit_layers=2L,source_validation_method="outer_fit",calibration_layout="compact",
  additional_aggregation_modes=c("separate_arms","joint_tate"),M_tau=12,M_tau_inference=12,
  aggregation_lambda=.5,checkpoint_dir=file.path(output,"checkpoints"),verbose=FALSE)
results <- list(); references <- list(); source_estimates <- list()
for(variant in c("full","standard_source","ordinary_target")) {
  inputs <- base_arguments
  inputs$target_nuisance_method <- if(variant=="ordinary_target") "lasso" else "hou_calibrated"
  inputs$calibration_control <- list(recipe="score_derivative",
    target_propensity_initialization=if(variant=="ordinary_target") "logistic" else "calibrated",
    source_nuisance_method=if(variant=="standard_source") "standard" else "calibrated")
  if(variant!="ordinary_target") inputs$calibration_control$target_radius <- 12
  started <- proc.time()[["elapsed"]]
  result <- do.call(run_single_simulation,inputs)
  artifacts <- attr(result,"roce_simulation_artifacts")
  data <- artifacts$data_split
  stopifnot(length(data)==9L,all(vapply(data,function(site) site$n==1000L,logical(1L))),
    all(is.finite(result$estimate)),all(is.finite(result$se)),all(result$se>0),
    all(result$roce_crossfit_levels==2L))
  fit <- artifacts$direct_tate_results$one_round_crossfit
  for(mode in c("common_tate","separate_arms","joint_tate")) {
    selected <- if(mode=="common_tate") fit else reaggregate_tate_crossfit(data,fit,
      M_tau_inference=12,aggregation_mode=mode,verbose=FALSE)
    detail <- summarize_layer_arms(selected,data)
    stopifnot(all(is.finite(detail$covariance)))
  }
  reference <- result[result$method=="target_only_ate",c("estimate","se","truth")]
  references[[variant]] <- as.numeric(reference[1L,])
  source_estimates[[variant]] <- fit$source_estimates
  data_hash <- digest::digest(data,algo="sha256")
  if(variant=="full") {
    expected_hash <- data_hash
    expected_reference <- references[[variant]]
  } else {
    stopifnot(identical(data_hash,expected_hash),
              max(abs(references[[variant]]-expected_reference))<1e-12)
  }
  saveRDS(result,file.path(output,paste0(variant,".rds")))
  attr(result,"roce_simulation_artifacts") <- NULL
  result$variant <- variant
  result$check_elapsed_seconds <- proc.time()[["elapsed"]]-started
  results[[variant]] <- result
  cat("Checked",variant,"K8, p10, C2, rho1.5, two layers, ten folds, grid100\n")
}
stopifnot(max(abs(source_estimates$full-source_estimates$ordinary_target))<1e-12)
write.csv(do.call(rbind,results),file.path(output,"results.csv"),row.names=FALSE)
writeLines(c("K8_TWO_LAYER_EXECUTION_CHECKS_PASSED",
  "1000 observations at each of9 sites; p10 C2 rho1.5; seed701 outside the main1:200 panel.",
  "Complete source-calibrated, standard-source and ordinary-target variants; all three aggregations checked.",
  "Data hashes, ordinary-target references, source estimates across target anchors, and arm/TATE variance identities checked.",
  "Execution checks only; no coverage or calibration-benefit inference."),file.path(output,"checks.txt"))
