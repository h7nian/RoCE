#!/usr/bin/env Rscript
# Conditional population moments for the saved p4, rho0 fitted-score fixtures.
# Quadrature is a diagnostic oracle, not an input to any feasible estimator.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: checked_R_library repository_root new_scratch_output")
library_path <- normalizePath(arguments[[1L]],mustWork=TRUE)
repository <- normalizePath(arguments[[2L]],mustWork=TRUE)
output <- arguments[[3L]]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/next_paper/candidate_score_summary.R"))
source(file.path(repository,"diagnosis/next_paper/honest_candidate_scores.R"))
dir.create(output,recursive=TRUE)

source(file.path(repository,"diagnosis/next_paper/honest_population_moments.R"))

fixtures <- data.frame(scenario=c("C2","C3","C1"),sources=c(1L,2L,2L),version=c(1L,3L,4L))
root <- "/scratch.global/zhan9381/FACE-HD/implementation/next_paper/v9"
rows <- list(); outputs <- list(); errors <- list()
for(index in seq_len(nrow(fixtures))) {
  case <- fixtures[index,]
  path <- file.path(root,paste0("honest_candidate_checks_v",case$version),"checked_candidates.rds")
  checked <- readRDS(path); fitted <- checked$fitted; packet <- checked$packet
  set.seed(49231)
  dgp <- generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,case$sources),p=4L,config=case$scenario)
  data <- split_data_by_site(dgp)
  # The saved fit must be from this exact DGP draw, not merely a matching label.
  invisible(evaluate_honest_candidates(fitted,data))
  target_rows <- which(dgp$R=="t")
  stopifnot(max(abs(plogis(drop(dgp$Z_site_true[target_rows,1:4]%*%c(.35,-.25,.125,-.0625)))-
    dgp$p_treat_true[target_rows]))<1e-14)
  coarse <- honest_population_moments(fitted,packet,dgp,16L)
  fine <- honest_population_moments(fitted,packet,dgp,24L)
  difference <- max(abs(coarse$means-fine$means),abs(coarse$covariance-fine$covariance))
  stopifnot(difference<1e-8)
  oracle_variance_error <- NA_real_
  if(case$scenario=="C1") {
    oracle_fit <- fitted; oracle_packet <- packet
    oracle_packet$training_shifts[,] <- 0
    for(arm in c("mu1","mu0")) {
      value <- as.integer(arm=="mu1")
      alpha <- if(value==1L) dgp$alpha1_true else dgp$alpha0_true
      oracle_fit$models[[arm]]$target$outcome <- alpha
      oracle_fit$models[[arm]]$target$weight <- c(0,(2*value-1)*c(.35,-.25,.125,-.0625))
      oracle_fit$models[[arm]]$common$coefficients <- alpha
      for(site in names(oracle_fit$models[[arm]]$source)) {
        oracle_fit$models[[arm]]$source[[site]]$outcome <- alpha
        oracle_fit$models[[arm]]$source[[site]]$weight <- dgp$gamma_params[[paste0(site,"_",value)]]
      }
    }
    oracle <- honest_population_moments(oracle_fit,oracle_packet,dgp,24L)
    stopifnot(max(abs(sweep(oracle$means,2L,oracle$truth,"-")))<1e-12)
    quad <- RoCE:::.bounded_quadrature(24L)
    design <- cbind(1,quad$X)
    means <- cbind(plogis(drop(design%*%dgp$alpha1_true)),plogis(drop(design%*%dgp$alpha0_true)))
    propensity <- plogis(drop(quad$X%*%c(.35,-.25,.125,-.0625)))
    target_effect <- oracle$truth[["mu1"]]-oracle$truth[["mu0"]]
    common_variance <- sum(quad$weights*(means[,1L]-means[,2L]-target_effect)^2)/packet$evaluation_sizes[["t"]]
    anchor_variance <- common_variance+sum(quad$weights*(means[,1L]*(1-means[,1L])/propensity+
      means[,2L]*(1-means[,2L])/(1-propensity)))/packet$evaluation_sizes[["t"]]
    true_variances <- c(target_anchor=anchor_variance)
    for(site in names(oracle_fit$models$mu1$source)) {
      weights <- sapply(1:0,function(value) exp(-drop(design%*%dgp$gamma_params[[paste0(site,"_",value)]])))
      true_variances[site] <- common_variance+sum(quad$weights*rowSums(weights*means*(1-means)))/packet$evaluation_sizes[[site]]
    }
    width <- nrow(oracle$means)
    calculated <- vapply(seq_len(width),function(site) {
      contrast <- numeric(2L*width);contrast[c(site,width+site)] <- c(1,-1)
      drop(crossprod(contrast,oracle$covariance%*%contrast))
    },numeric(1L))
    oracle_variance_error <- max(abs(calculated-true_variances))
    stopifnot(oracle_variance_error<1e-12)
  }
  # A paired research alternative: use target IPW OR directly in every source
  # residual, retaining the identical fitted merged weights. No extra fitting
  # or training-mean translation is needed for this evaluation comparison.
  common_fit <- fitted
  for(arm in c("mu1","mu0")) for(site in names(common_fit$models[[arm]]$source)) {
    common_fit$models[[arm]]$source[[site]]$outcome <- common_fit$models[[arm]]$common$coefficients
  }
  common_packet <- evaluate_honest_candidates(common_fit,data)
  stopifnot(max(abs(common_packet$training_shifts))==0)
  common_population <- honest_population_moments(common_fit,common_packet,dgp,24L)
  common_coarse <- honest_population_moments(common_fit,common_packet,dgp,16L)
  common_difference <- max(abs(common_population$means-common_coarse$means),
    abs(common_population$covariance-common_coarse$covariance))
  stopifnot(common_difference<1e-8)
  variants <- list(translated_source=list(population=fine,packet=packet),
                   common_outcome=list(population=common_population,packet=common_packet))
  width <- nrow(fine$means)
  for(construction in names(variants)) for(site in seq_len(width)) for(quantity in c("mu1","mu0","tate")) {
    population <- variants[[construction]]$population
    evaluated <- variants[[construction]]$packet
    contrast <- numeric(2L*width)
    if(quantity!="mu0") contrast[site] <- 1
    if(quantity!="mu1") contrast[width+site] <- if(quantity=="tate") -1 else 1
    arm_contrast <- if(quantity=="mu1") c(1,0) else if(quantity=="mu0") c(0,1) else c(1,-1)
    expected <- sum(population$means[site,]*arm_contrast)
    truth <- sum(population$truth*arm_contrast)
    variance <- drop(crossprod(contrast,population$covariance %*% contrast))
    empirical_variance <- drop(crossprod(contrast,evaluated$covariance %*% contrast))
    rows[[length(rows)+1L]] <- data.frame(config=case$scenario,K=case$sources,construction=construction,
      protocol=fitted$configuration$communication_mode,site=rownames(fine$means)[site],quantity=quantity,
      conditional_mean=expected,truth=truth,conditional_bias=expected-truth,
      original_conditional_bias=sum(fine$original_means[site,]*arm_contrast)-truth,
      conditional_sd=sqrt(variance),bias_over_conditional_sd=(expected-truth)/sqrt(variance),
      estimated_to_conditional_variance=empirical_variance/variance,
      evaluation_z=(sum(evaluated$means[site,]*arm_contrast)-expected)/sqrt(variance))
  }
  outputs[[case$scenario]] <- variants
  errors[[case$scenario]] <- data.frame(config=case$scenario,quadrature_difference=difference,
    replacement_identity_error=fine$identity_error,oracle_variance_error=oracle_variance_error,
    common_outcome_quadrature_difference=common_difference)
}
write.csv(do.call(rbind,rows),file.path(output,"conditional_moments.csv"),row.names=FALSE)
write.csv(do.call(rbind,errors),file.path(output,"oracle_checks.csv"),row.names=FALSE)
saveRDS(outputs,file.path(output,"population_packets.rds"))
writeLines(c("CONDITIONAL_POPULATION_AUDIT_PASSED",
  "Three saved training datasets, p4,rho0,n1000/site. 16/24 quadrature agreement checked.",
  "Conditional nuisance bias and variance diagnostic; no repeated-training coverage or growing-K claim."),
  file.path(output,"checks.txt"))
cat("CONDITIONAL_POPULATION_AUDIT_PASSED\n")
