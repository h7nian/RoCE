#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=4L)stop("usage: run_matched_outcome_cv.R LIBRARY VERSION PILOT_ROOT OUTPUT")
  version<-args[2];stopifnot(version %in% c("v19","corrected"))
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  source("diagnosis/tate_common_weight/shared_calibration_basis.R")
  fingerprint<-.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)
  expected<-if(version=="v19")"2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7" else
    "6a715f3a41f2f42ef58bc08310d4ecca061524b0ab555d52a6221101a63a0a93"
  stopifnot(fingerprint==expected)
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));h<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(h)==1L,roce_sha256_file(file.path(root,name))==h);readRDS(file.path(root,name))
  }
  artifact<-checked(file.path(args[3],"seed_010013"),"artifacts.rds")$group_result$artifacts[["0"]]
  saved<-checked(file.path(args[3],"shared_final_C3_v1"),"shared_final_states.rds")
  data<-.source_working_basis(artifact$data_split,"C3")
  U<-.apply_shared_calibration_basis(saved$basis_spec,data$s1$W_outcome,data$s1$Z_site)
  train<-unlist(saved$source_ids[2:5],use.names=FALSE)
  X<-U[train,,drop=FALSE];A<-data$s1$A[train];Y<-data$s1$Y[train]
  plugin<-RoCE:::.make_plugin_block_design(lapply(saved$source_ids[2:5],function(ids)U[ids,,drop=FALSE]),caller="matched outcome CV")
  initial<-saved$results$mu1
  gamma<-c(0,unlist(initial$gamma_plugin,use.names=FALSE))
  n_cv<-RoCE:::.nuisance_cv_fold_count(A,1L,"matched outcome CV")
  folds<-RoCE:::with_seed(994102L,sample(rep(seq_len(n_cv),length.out=sum(A==1L))))
  lambda_max<-RoCE:::compute_lambda_max_outcome(X,Y,A,gamma,A_val=1L,
    family_int=1L,link_int=1L,Z_site=plugin,calibrated=TRUE,M_tau=5)
  grid<-RoCE:::build_lambda_grid(lambda_max=lambda_max,
    lambda_min_ratio=RoCE:::LAMBDA_MIN_RATIO_LOW_DIM,nlambda=100L)
  inputs<-list(X=X,A=A,Y=Y,plugin=plugin,gamma=gamma,lambda_grid=grid,fold_id=folds,
    source_ids=train,M_fit=5,arm=1L,seed=994102L)
  warnings<-character()
  answer<-tryCatch(withCallingHandlers({
    cv<-RoCE:::select_lambda_cv_calibrated_outcome_cpp(X,Y,A,gamma,grid,n_cv,
      RoCE:::MAX_ITER_DEFAULT,RoCE:::TOL_DEFAULT,1L,5,plugin,1L,1L,folds)
    fit<-RoCE:::fit_unified_outcome_cpp(X,Y,A,gamma,1L,1L,cv$lambda_min,
      RoCE:::MAX_ITER_DEFAULT,RoCE:::TOL_DEFAULT,1L,TRUE,5,plugin,Reduce(`+`,initial$alpha_plugin)/4)
    list(cv=cv,fit=fit)
  },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),error=identity)
  failed<-inherits(answer,"error")
  roce_write_atomic_directory(args[4],function(stage) {
    saveRDS(inputs,file.path(stage,"inputs.rds"))
    saveRDS(list(result=answer,warnings=warnings),file.path(stage,"cv_result.rds"))
    if(!failed)write.csv(data.frame(lambda=grid,cv_score=as.numeric(answer$cv$cv_scores),
      cv_se=as.numeric(answer$cv$cv_se),n_valid=as.integer(answer$cv$n_valid_folds)),file.path(stage,"cv_path.csv"),row.names=FALSE)
    writeLines(c(paste0("version=",version),paste0("package_fingerprint=",fingerprint),
      paste0("failed=",failed),paste0("warning_count=",length(warnings)),
      paste0("cv_threads=",Sys.getenv("ROCE_NUISANCE_CV_THREADS","1")),
      "inference_validated=FALSE","initial_nuisance_refits=0","bootstrap_draws=0",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/run_matched_outcome_cv.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="matched high-dimensional outcome CV")
  if(failed)stop(conditionMessage(answer))
  stopifnot(isTRUE(answer$fit$converged))
  cat(version,"lambda_min",answer$cv$lambda_min,"active",sum(answer$fit$alpha[-1]!=0),"\n")
}
if(sys.nframe()==0L)main()
