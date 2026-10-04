#!/usr/bin/env Rscript
# Verify the paper's two-arm covariance formulation in an independent basis.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: checked_R_library repository new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
repository <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/repair/calibration_score_diagnostics.R"))
assert_close <- function(first,second,label,tolerance=1e-12) {
  if(length(first)!=length(second) || any(!is.finite(first)) || any(!is.finite(second)) ||
     max(abs(first-second))>tolerance) stop(label)
}
base <- "/scratch.global/zhan9381/FACE-HD/implementation/r11/current_paper_mc200_v1"
records <- list()
for(K in c(2L,8L)) {
  root <- file.path(base,paste0("main_method_p100_k",K,"_s1_200"))
  manifest <- read.csv(file.path(root,"manifest.csv"),stringsAsFactors=FALSE)
  selected <- subset(manifest,sim_id==1L & rho %in% c(0,1.5) & config %in% c("C1","C2","C3"))
  stopifnot(nrow(selected)==6L)
  for(task_index in seq_len(nrow(selected))) {
    task <- selected[task_index,]
    directory <- file.path(root,"tasks",task$task_id)
    stopifnot(file.exists(file.path(directory,"COMPLETE")))
    packet <- readRDS(file.path(directory,"score_derivative_layer_records.rds"))
    packet$data_hash <- readLines(file.path(directory,"data_sha256.txt"),warn=FALSE)
    learning <- fit_packet_joint_weights(packet,aggregation_cutoff=2)
    m <- K+1L
    first <- seq_len(m); second <- m+first
    max_covariance_error <- max_objective_error <- 0
    for(fold in seq_len(packet$n_folds)) {
      arms <- lapply(packet$arms,function(arm)
        RoCE:::.inner_fold_records(arm$inner_fold_info[[fold]]))
      sites <- list(Map(cbind,arms$mu1$target,arms$mu0$target))
      for(j in seq_len(K)) sites[[j+1L]] <- Map(function(u1,u0) {
        values <- matrix(0,length(u1),2L*m)
        values[,j+1L] <- u1;values[,m+j+1L] <- u0
        values
      },arms$mu1$source[[j]],arms$mu0$source[[j]])
      sigma <- matrix(0,2L*m,2L*m)
      means <- numeric(2L*m)
      for(blocks in sites) {
        count <- sum(vapply(blocks,nrow,integer(1L)))
        means <- means+colSums(do.call(rbind,blocks))/count
        # The production training moments center each calibration fold before pooling.
        for(block in blocks) {
          centered <- sweep(block,2L,colMeans(block),"-")
          sigma <- sigma+crossprod(centered)/count^2
        }
      }
      transform <- matrix(0,2L*m,1L+2L*K)
      transform[1L,1L] <- 1;transform[m+1L,1L] <- -1
      for(j in seq_len(K)) {
        transform[1L,1L+j] <- -1;transform[1L+j,1L+j] <- 1
        transform[m+1L,1L+K+j] <- 1;transform[m+1L+j,1L+K+j] <- -1
      }
      projected <- crossprod(transform,sigma%*%transform)
      moments <- learning$moments[[fold]]
      reference <- rbind(c(moments$constant,moments$l),cbind(moments$l,moments$Q))
      assert_close(projected,reference,"Candidate covariance does not reproduce increment covariance")
      max_covariance_error <- max(max_covariance_error,max(abs(projected-reference)))
      discrepancies <- c(means[first][-1L]-means[first][1L],means[second][-1L]-means[second][1L])
      assert_close(c(discrepancies[seq_len(K)],-discrepancies[K+seq_len(K)]),
                   moments$discrepancy,"Control-arm discrepancy sign mismatch")
      for(eta in list(learning$fits[[fold]]$weights,seq(-.2,.4,length.out=2L*K))) {
        w1 <- c(1-sum(eta[seq_len(K)]),eta[seq_len(K)])
        w0 <- c(1-sum(eta[K+seq_len(K)]),eta[K+seq_len(K)])
        variance <- drop(crossprod(w1,sigma[first,first]%*%w1)+
                         crossprod(w0,sigma[second,second]%*%w0)-
                         2*crossprod(w1,sigma[first,second]%*%w0))
        quadratic <- moments$constant+2*sum(moments$l*eta)+drop(crossprod(eta,moments$Q%*%eta))
        assert_close(variance,quadratic,"Two-arm variance does not reproduce joint quadratic")
        scale <- c(diag(sigma[first,first])[-1L]+sigma[1L,1L]-2*sigma[1L,first[-1L]],
                   diag(sigma[second,second])[-1L]+sigma[m+1L,m+1L]-2*sigma[m+1L,second[-1L]])
        stopifnot(all(scale>RoCE:::VARIANCE_MIN))
        penalty <- pmax(abs(discrepancies)/sqrt(scale)/2-1,0)
        assert_close(penalty,learning$fits[[fold]]$penalty,"Arm-specific penalty mismatch")
        if(identical(eta,learning$fits[[fold]]$weights)) {
          objective <- moments$N_all*variance+sum(penalty*abs(eta))
          assert_close(objective,learning$fits[[fold]]$objective,"Full manuscript objective differs",1e-10)
          max_objective_error <- max(max_objective_error,abs(objective-learning$fits[[fold]]$objective))
        }
      }
    }
    site_scores <- lapply(packet$site_sizes,function(size) matrix(NA_real_,size,2L))
    for(arm_index in 1:2) {
      arm <- c("mu1","mu0")[arm_index]
      for(fold in seq_len(packet$n_folds)) {
        info <- packet$arms[[arm]]$fold_info[[fold]]
        eta <- learning$weights[[arm]][fold,]
        target <- (1-sum(eta))*(info$varphi_ot+info$fold_target_estimate)
        for(j in seq_len(K)) {
          target <- target+eta[j]*(info$zeta_components[[j]]+info$mu_pred_ts[j])
          site_scores[[j+1L]][info$source_idx[[j]],arm_index] <- eta[j]*(info$xi_components[[j]]+info$delta_ts[j])
        }
        site_scores[[1L]][info$target_idx,arm_index] <- target
      }
    }
    stopifnot(all(vapply(site_scores,function(value) all(is.finite(value)),logical(1L))))
    point <- Reduce(`+`,lapply(site_scores,colMeans))
    covariance <- Reduce(`+`,lapply(site_scores,function(value)
      crossprod(sweep(value,2L,colMeans(value),"-"))/nrow(value)^2))
    reconstructed <- evaluate_packet_weights(packet,learning$weights,learning)
    assert_close(point,reconstructed$estimates[c("mu1","mu0")],"Data-level arm formula mismatch")
    assert_close(covariance,reconstructed$covariance_fixed_weights,"Data-level covariance mismatch")
    values <- read.csv(file.path(directory,"results.csv"),stringsAsFactors=FALSE)
    result <- subset(values,method=="one_round_crossfit_ate_joint_tate")
    stopifnot(nrow(result)==1L)
    assert_close(point[1L]-point[2L],result$estimate,"Full TATE expression mismatch")
    assert_close(reconstructed$tate_variance,result$se^2,"Reported eta-adjusted variance mismatch")
    records[[length(records)+1L]] <- data.frame(config=task$config,K=K,rho=task$rho,sim_id=task$sim_id,
      folds=packet$n_folds,max_covariance_error=max_covariance_error,max_objective_error=max_objective_error,
      covariance_mu1_mu0=covariance[1L,2L],point_estimate=point[1L]-point[2L])
  }
}
# A positive-definite known-covariance counterexample to separate-arm equivalence.
within <- diag(c(1,2));cross <- matrix(c(0,.1,.8,0),2L,2L)
sigma <- rbind(cbind(within,cross),cbind(t(cross),within))
stopifnot(min(eigen(sigma,symmetric=TRUE,only.values=TRUE)$values)>0)
transform <- cbind(c(1,0,-1,0),c(-1,1,0,0),c(0,0,1,-1))
projected <- crossprod(transform,sigma%*%transform)
joint <- drop(solve(projected[-1L,-1L],-projected[-1L,1L]))
separate <- rep(1/3,2L)
variance <- function(eta) drop(crossprod(c(1,eta),projected%*%c(1,eta)))
assert_close(joint,c(8/39,7/13),"Independent example joint optimum mismatch")
stopifnot(variance(joint)<variance(separate),max(abs(joint-separate))>.1)
quadrature <- RoCE:::.bounded_quadrature()
oracle <- lapply(c("C1","C2","C3"),function(config) {
  features <- RoCE:::.bounded_active_features(quadrature$X,if(config=="C2") .75 else 0)
  predictor <- drop(features%*%c(.25,.20,.10,.05))
  treated <- plogis(.5+predictor);control <- plogis(-.5+predictor)
  mean1 <- sum(quadrature$weights*treated);mean0 <- sum(quadrature$weights*control)
  covariance <- sum(quadrature$weights*(treated-mean1)*(control-mean0))
  stopifnot(covariance>0)
  data.frame(config=config,mu1=mean1,mu0=mean0,tate=mean1-mean0,
    covariance_predictions=covariance,oracle_aipw_covariance_n1000=covariance/1000)
})
dir.create(output,recursive=TRUE)
write.csv(do.call(rbind,records),file.path(output,"fit_checks.csv"),row.names=FALSE)
write.csv(do.call(rbind,oracle),file.path(output,"oracle_cross_arm_covariance.csv"),row.names=FALSE)
write.csv(data.frame(method=c("joint","separate"),eta1=c(joint[1L],separate[1L]),
  eta0=c(joint[2L],separate[2L]),variance=c(variance(joint),variance(separate))),
  file.path(output,"independent_covariance_example.csv"),row.names=FALSE)
writeLines(c("JOINT_TATE_MANUSCRIPT_CHECKS_PASSED",
  "Twelve saved fits, 120 outer training problems at rho0/1.5: covariance-basis transformation, objective, arm penalty, full site-weighted estimate, fixed-weight covariance and reported variance match.",
  "A positive-definite example gives different jointly/independently optimized weights.",
  "Known-nuisance target AIPW cross-arm covariance is positive under all three bounded DGP scenarios despite iid patients."),file.path(output,"checks.txt"))
print(do.call(rbind,records));print(do.call(rbind,oracle))
