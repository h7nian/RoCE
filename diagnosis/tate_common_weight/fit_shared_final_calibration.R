#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=5L)stop("usage: fit_shared_final_calibration.R LIBRARY SEED PILOT_ROOT CONFIG OUTPUT")
  config<-args[4];if(!config %in% c("C2","C3"))stop("CONFIG must be C2 or C3")
  if(file.exists(args[5]))stop("output already exists")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/shared_calibration_basis.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  code<-file.path("diagnosis/tate_common_weight",c("fit_shared_final_calibration.R",
    "shared_calibration_basis.R","source_working_basis.R","SHARED_FINAL_CALIBRATION_PROTOCOL.md"))
  hashes<-vapply(code,roce_sha256_file,"")
  stopifnot(.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));expected<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  artifact<-checked(args[2],"artifacts.rds")$group_result$artifacts[["0"]]
  original_root<-file.path(args[3],paste0("source_basis_",config,"_outer_1_v1"))
  original<-checked(original_root,"source_validation_states.rds")
  stopifnot(identical(original$global_fold_map,1:5),identical(original$working_config,config))
  data<-.source_working_basis(artifact$data_split,config)
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  target_ids<-lapply(info,`[[`,"target_idx");source_ids<-lapply(info,function(x)x$source_idx[[1L]])
  training_ids<-unlist(target_ids[2:5],use.names=FALSE)
  spec<-.fit_shared_calibration_basis(data$t$W_outcome[training_ids,,drop=FALSE],data$t$Z_site[training_ids,,drop=FALSE],
    if(config=="C2")"site" else "outcome")
  target_U<-.apply_shared_calibration_basis(spec,data$t$W_outcome,data$t$Z_site)
  source_U<-.apply_shared_calibration_basis(spec,data$s1$W_outcome,data$s1$Z_site)
  source_pieces<-lapply(source_ids[2:5],function(index)source_U[index,,drop=FALSE])
  stacked_ids<-unlist(source_ids[2:5],use.names=FALSE)
  stack<-source_U[stacked_ids,,drop=FALSE]
  plugin<-RoCE:::.make_plugin_block_design(source_pieces,caller="shared final calibration candidate")
  results<-list()
  for(arm in 0:1) {
    old<-original$results[[paste0("mu",arm)]];stopifnot(!old$failed)
    fit<-old$result; keys<-paste0("k2_",2:5)
    stopifnot(identical(names(fit$per_k2_alpha),keys),identical(names(fit$per_k2_gamma),keys))
    alpha_plugin<-lapply(fit$per_k2_alpha,function(a)drop(spec$outcome_map%*%a))
    gamma_plugin<-lapply(fit$per_k2_gamma,function(g)drop(spec$site_map%*%g))
    mean_gradients<-lapply(2:5,function(k)RoCE:::.mean_glm_gradient_site_basis(
      data$t$W_outcome[target_ids[[k]],,drop=FALSE],target_U[target_ids[[k]],,drop=FALSE],
      fit$per_k2_alpha[[paste0("k2_",k)]],1L,1L))
    prediction_error<-max(vapply(2:5,function(k) {
      key<-paste0("k2_",k);idx<-source_ids[[k]]
      max(abs(drop(cbind(1,source_U[idx,,drop=FALSE])%*%alpha_plugin[[key]])-
        drop(cbind(1,data$s1$W_outcome[idx,,drop=FALSE])%*%fit$per_k2_alpha[[key]])),
        abs(drop(cbind(1,source_U[idx,,drop=FALSE])%*%gamma_plugin[[key]])-
        drop(cbind(1,data$s1$Z_site[idx,,drop=FALSE])%*%fit$per_k2_gamma[[key]])))
    },numeric(1L)))
    stopifnot(prediction_error<1e-10)
    warnings<-character()
    answer<-tryCatch(RoCE:::with_seed(883101L+arm,withCallingHandlers({
      gamma<-RoCE:::fit_unified_density_ratio(Z_site=stack,A=data$s1$A[stacked_ids],
        mean_grad_psi=Reduce(`+`,mean_gradients)/4,alpha_init=c(0,unlist(alpha_plugin,use.names=FALSE)),
        lambda=NULL,calibrated=TRUE,M_tau=5,W_outcome=plugin,A_val=arm,family_int=1L,link_int=1L,
        warm_start=Reduce(`+`,gamma_plugin)/4,nlambda=100L,lambda_rule="min")
      alpha<-RoCE:::fit_unified_outcome(W_outcome=stack,Y=data$s1$Y[stacked_ids],A=data$s1$A[stacked_ids],
        A_val=arm,gamma_s=c(0,unlist(gamma_plugin,use.names=FALSE)),lambda=NULL,family_int=1L,link_int=1L,
        calibrated=TRUE,M_tau=5,Z_site=plugin,warm_start=Reduce(`+`,alpha_plugin)/4,nlambda=100L,lambda_rule="min")
      list(alpha=alpha,gamma=gamma)
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})),error=identity)
    failed<-inherits(answer,"error")
    results[[paste0("mu",arm)]]<-list(failed=failed,failure=if(failed)conditionMessage(answer)else NULL,
      warnings=warnings,final_fit=if(failed)NULL else answer,initial_alpha=fit$per_k2_alpha,initial_gamma=fit$per_k2_gamma,
      alpha_plugin=alpha_plugin,gamma_plugin=gamma_plugin,mean_gradients=mean_gradients,prediction_error=prediction_error)
  }
  stopifnot(identical(hashes,vapply(code,roce_sha256_file,"")))
  roce_write_atomic_directory(args[5],function(stage) {
    saveRDS(list(results=results,basis_spec=spec,target_ids=target_ids,source_ids=source_ids,working_config=config),
      file.path(stage,"shared_final_states.rds"))
    write.csv(data.frame(arm=0:1,failed=vapply(results,`[[`,logical(1L),"failed"),
      warning_count=vapply(results,function(x)length(x$warnings),integer(1L))),file.path(stage,"fit_status.csv"),row.names=FALSE)
    write.csv(data.frame(path=code,sha256=hashes),file.path(stage,"code_hashes.csv"),row.names=FALSE)
    writeLines(c("production_method_changed=FALSE","inference_validated=FALSE","initial_nuisance_refits=0",
      paste0("original_fit_manifest_sha256=",roce_sha256_file(file.path(original_root,"sha256.txt")))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="shared final basis native calibration")
  if(any(vapply(results,`[[`,logical(1L),"failed")))stop("failed final fits retained; review before evaluation")
}
if(sys.nframe()==0L)main()
