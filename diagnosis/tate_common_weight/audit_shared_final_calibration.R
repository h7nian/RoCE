#!/usr/bin/env Rscript
.shared_final_losses <- function(U,A,Y,outcome_weight,tilt_weight,target_summary,alpha,gamma) {
  eta<-drop(U%*%alpha); log_weight<-drop(U%*%gamma)
  n<-nrow(U)
  list(outcome_loss=mean(A*outcome_weight*(pmax(eta,0)+log1p(exp(-abs(eta)))-Y*eta)),
    tilting_loss=sum(target_summary*gamma)+mean(A*tilt_weight*exp(-log_weight)),
    outcome_gradient=drop(crossprod(U,A*outcome_weight*(plogis(eta)-Y)))/n,
    tilting_gradient=target_summary-drop(crossprod(U,A*tilt_weight*exp(-log_weight)))/n)
}

main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=5L)stop("usage: audit_shared_final_calibration.R LIBRARY SEED PILOT_ROOT CONFIG OUTPUT")
  config<-args[4]
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/shared_calibration_basis.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));expected<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  root<-file.path(args[3],paste0("shared_final_",config,"_v1"))
  saved<-checked(root,"shared_final_states.rds")
  original<-checked(file.path(args[3],paste0("source_basis_",config,"_outer_1_v1")),"source_validation_states.rds")
  artifact<-checked(args[2],"artifacts.rds")$group_result$artifacts[["0"]]
  data<-.source_working_basis(artifact$data_split,config)
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  stopifnot(saved$working_config==config,identical(saved$target_ids,lapply(info,`[[`,"target_idx")),
    identical(saved$source_ids,lapply(info,function(x)x$source_idx[[1L]])))
  target_U<-.apply_shared_calibration_basis(saved$basis_spec,data$t$W_outcome,data$t$Z_site)
  source_U<-.apply_shared_calibration_basis(saved$basis_spec,data$s1$W_outcome,data$s1$Z_site)
  target_ids<-saved$target_ids;source_ids<-saved$source_ids
  train<-unlist(source_ids[2:5],use.names=FALSE);U<-cbind(1,source_U[train,,drop=FALSE])
  reports<-states<-list()
  clip<-function(x)pmax(-5,pmin(5,x))
  for(arm in 0:1) {
    fitted<-saved$results[[paste0("mu",arm)]];old<-original$results[[paste0("mu",arm)]]$result
    stopifnot(!fitted$failed,identical(fitted$initial_alpha,old$per_k2_alpha),
      identical(fitted$initial_gamma,old$per_k2_gamma))
    alpha<-fitted$final_fit$alpha;gamma<-fitted$final_fit$gamma
    source_outcome_weight<-source_tilt_weight<-numeric();summaries<-list()
    initial_clipping<-c(outcome=0L,tilting=0L)
    for(k in 2:5) {
      key<-paste0("k2_",k);s<-source_ids[[k]];t<-target_ids[[k]]
      alpha_initial<-old$per_k2_alpha[[key]];gamma_initial<-old$per_k2_gamma[[key]]
      a_source<-drop(cbind(1,data$s1$W_outcome[s,,drop=FALSE])%*%alpha_initial)
      a_target<-drop(cbind(1,data$t$W_outcome[t,,drop=FALSE])%*%alpha_initial)
      stopifnot(max(abs(a_target))<RoCE:::LOGISTIC_CLIP-1)
      g_source<-drop(cbind(1,data$s1$Z_site[s,,drop=FALSE])%*%gamma_initial)
      m_source<-plogis(clip(a_source));m_target<-plogis(a_target)
      summaries[[key]]<-colMeans(cbind(1,target_U[t,,drop=FALSE])*(m_target*(1-m_target)))
      source_outcome_weight<-c(source_outcome_weight,exp(-clip(g_source)))
      source_tilt_weight<-c(source_tilt_weight,m_source*(1-m_source))
      initial_clipping<-initial_clipping+c(sum(abs(a_source)>5),sum(abs(g_source)>5))
    }
    target_summary<-Reduce(`+`,summaries)/4
    stopifnot(max(abs(target_summary-Reduce(`+`,fitted$mean_gradients)/4))<1e-12)
    evaluate<-function(a,g).shared_final_losses(U,as.numeric(data$s1$A[train]==arm),data$s1$Y[train],
      source_outcome_weight,source_tilt_weight,target_summary,a,g)
    value<-evaluate(alpha,gamma)
    kkt<-function(parameter,gradient) {
      lambda<-attr(parameter,"lambda_used");stopifnot(length(lambda)==1L,is.finite(lambda))
      penalty<-c(0,rep(lambda,length(parameter)-1L));active<-parameter!=0
      residual<-pmax(abs(gradient)-penalty,0)
      residual[active]<-abs(gradient[active]+penalty[active]*sign(parameter[active]))
      max(residual)
    }
    kkt_error<-max(kkt(alpha,value$outcome_gradient),kkt(gamma,value$tilting_gradient))
    derivative_error<-0
    for(k in 1:6) {
      direction<-sin(seq_along(alpha)*k);direction<-direction/sqrt(sum(direction^2));eps<-1e-6
      derivative_error<-max(derivative_error,
        abs((evaluate(alpha+eps*direction,gamma)$outcome_loss-evaluate(alpha-eps*direction,gamma)$outcome_loss)/(2*eps)-sum(value$outcome_gradient*direction)),
        abs((evaluate(alpha,gamma+eps*direction)$tilting_loss-evaluate(alpha,gamma-eps*direction)$tilting_loss)/(2*eps)-sum(value$tilting_gradient*direction)))
    }
    s<-source_ids[[1L]];t<-target_ids[[1L]]
    native<-RoCE:::calculate_correction_term_cpp(source_U[s,,drop=FALSE],data$s1$A[s],data$s1$Y[s],
      gamma,alpha,source_U[s,,drop=FALSE],5,1L,1L,arm)
    target_score<-plogis(drop(cbind(1,target_U[t,,drop=FALSE])%*%alpha))
    source_score<-as.numeric(data$s1$A[s]==arm)*exp(-clip(drop(cbind(1,source_U[s,,drop=FALSE])%*%gamma)))*
      (data$s1$Y[s]-plogis(drop(cbind(1,source_U[s,,drop=FALSE])%*%alpha)))
    point_error<-abs(mean(source_score)-native$delta_ts)
    stopifnot(kkt_error<1e-4,derivative_error<1e-7,point_error<1e-12,!native$clip_diagnostics$any_safety_clipped)
    reports[[paste0("mu",arm)]]<-data.frame(arm,kkt_error,derivative_error,point_error,
      initial_outcome_clipped=initial_clipping[1],initial_tilting_clipped=initial_clipping[2],
      mean_estimate=mean(target_score)+mean(source_score),warning_count=length(fitted$warnings))
    states[[paste0("mu",arm)]]<-list(losses=value,target_score=target_score,source_score=source_score,
      inference_clipping=native$clip_diagnostics)
  }
  reports<-do.call(rbind,reports)
  roce_write_atomic_directory(args[5],function(stage) {
    write.csv(reports,file.path(stage,"calibration_checks.csv"),row.names=FALSE)
    saveRDS(states,file.path(stage,"audited_site_scores.rds"))
    writeLines(c("inference_validated=FALSE","original_initial_coefficients_unchanged=TRUE",
      paste0("fit_manifest_sha256=",roce_sha256_file(file.path(root,"sha256.txt"))),
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/audit_shared_final_calibration.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="shared final calibration equation audit")
  print(reports,row.names=FALSE)
}
if(sys.nframe()==0L)main()
