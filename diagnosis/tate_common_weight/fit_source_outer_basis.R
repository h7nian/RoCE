#!/usr/bin/env Rscript
# Full outer-training nuisance refit for a paired C2/C3 basis comparison.
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=4L) stop("usage: fit_source_outer_basis.R LIBRARY SEED CONFIG OUTPUT")
  config<-args[3]
  if(!config %in% c("C2","C3")) stop("CONFIG must be C2 or C3")
  if(file.exists(args[4])) stop("output already exists")
  .libPaths(c(args[1],.libPaths())); suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R"); source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  stopifnot(.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(args[2],"sha256.txt"))==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  lines<-readLines(file.path(args[2],"sha256.txt"))
  expected<-substr(lines[substring(lines,67L)=="artifacts.rds"],1L,64L)
  stopifnot(length(expected)==1L,roce_sha256_file(file.path(args[2],"artifacts.rds"))==expected)
  artifact<-readRDS(file.path(args[2],"artifacts.rds"))$group_result$artifacts[["0"]]
  data<-.source_working_basis(artifact$data_split,config)
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  ids<-list(target=lapply(info,`[[`,"target_idx"),source=lapply(info,function(x)x$source_idx[[1L]]))
  views<-function(site_data,indices) {
    folds<-lapply(indices,function(index)list(original_idx=index,n=length(index)))
    attr(folds,".data_ref")<-site_data
    folds
  }
  target_folds<-views(data$t,ids$target); source_folds<-list(s1=views(data$s1,ids$source))
  splits<-lapply(2:5,function(k) {
    record<-list(local_calibration_fold=k,global_calibration_fold=k)
    for(site in names(ids)) {
      record[[paste0(site,"_training_ids")]]<-unlist(ids[[site]][setdiff(2:5,k)],use.names=FALSE)
      record[[paste0(site,"_calibration_ids")]]<-ids[[site]][[k]]
      stopifnot(!any(c(record[[paste0(site,"_training_ids")]],ids[[site]][[k]]) %in% ids[[site]][[1L]]))
    }
    record
  })
  results<-list()
  for(arm in 0:1) {
    warnings<-character(); inputs<-initial<-list()
    answer<-tryCatch(RoCE:::with_seed(881101L+arm,withCallingHandlers({
      cached_lambda<-NULL
      for(k in 2:5) {
        training<-RoCE:::combine_folds(target_folds,setdiff(2:5,k))
        calibration<-RoCE:::materialize_fold(target_folds,k)
        alpha<-RoCE:::fit_initial_outcome(training$W_outcome,training$Y,training$A,arm,
          lambda=cached_lambda,nlambda=100L,family="binomial",lambda_rule="min",cv_group_id=training$cv_group_id)
        lambda<-attr(alpha,"lambda_used")
        if(is.null(cached_lambda) && length(lambda)==1L && is.finite(lambda)) cached_lambda<-lambda
        key<-paste0("k2_",k); initial[[key]]<-alpha
        inputs[[key]]<-list(alpha_init=alpha,mean_phi=c(1,colMeans(training$Z_site)),
          mean_grad_psi_init=RoCE:::.mean_glm_gradient_site_basis(calibration$W_outcome,
            calibration$Z_site,alpha,1L,1L))
      }
      RoCE:::process_source_site(s="s1",source_folds=source_folds,target_folds=target_folds,
        k1=1L,n_folds=5L,A_val=arm,M_tau=5,M_tau_inference=5,data_split=data,
        get_fold_inputs=function(site,k)inputs[[paste0("k2_",k)]],family_int=1L,link_int=1L,
        use_lambda_cache=TRUE,nuisance_nlambda=100L,nuisance_lambda_rule="min")
    },warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")})),error=identity)
    failed<-inherits(answer,"error")
    results[[paste0("mu",arm)]]<-list(failed=failed,warnings=warnings,
      failure=if(failed)conditionMessage(answer)else NULL,result=if(failed)NULL else answer,
      initial_outcomes=initial,target_summaries=inputs,seed=881101L+arm)
  }
  roce_write_atomic_directory(args[4],function(stage) {
    saveRDS(list(results=results,splits=splits,global_fold_map=1:5,working_config=config,
      validation_ids=lapply(ids,`[[`,1L),excluded_outer_ids=lapply(ids,`[[`,1L)),
      file.path(stage,"source_validation_states.rds"))
    write.csv(data.frame(arm=0:1,failed=vapply(results,`[[`,logical(1L),"failed"),
      warning_count=vapply(results,function(x)length(x$warnings),integer(1L))),file.path(stage,"fit_status.csv"),row.names=FALSE)
    writeLines(c("fit_scope=complete_outer_training","outer_evaluation_fold=1",
      "inference_validated=FALSE","original_saved_fit_identity_claimed=FALSE",
      paste0("working_config=",config),paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/fit_source_outer_basis.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="source full outer-training basis fit")
  if(any(vapply(results,`[[`,logical(1L),"failed")))stop("outer-training nuisance failures retained")
}
if(sys.nframe()==0L) main()
