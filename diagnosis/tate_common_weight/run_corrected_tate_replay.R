#!/usr/bin/env Rscript
# Full original TATE estimator replay after the outcome-CV normalization fix.
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: run_corrected_tate_replay.R LIBRARY SEED_BUNDLE OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  fingerprint<-.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)
  stopifnot(fingerprint=="6a715f3a41f2f42ef58bc08310d4ecca061524b0ab555d52a6221101a63a0a93",
    roce_sha256_file(file.path(args[2],"sha256.txt"))=="f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  manifest<-readLines(file.path(args[2],"sha256.txt"))
  expected<-substr(manifest[substring(manifest,67L)=="artifacts.rds"],1L,64L)
  stopifnot(length(expected)==1L,roce_sha256_file(file.path(args[2],"artifacts.rds"))==expected)
  artifact<-readRDS(file.path(args[2],"artifacts.rds"))$group_result$artifacts[["0"]]
  data<-artifact$data_split;old<-artifact$direct_tate_results$one_round_crossfit
  info<-old$intermediates$fold_info
  views<-function(site,ids) {
    stopifnot(identical(sort(unlist(ids,use.names=FALSE)),seq_len(site$n)))
    result<-lapply(ids,function(index)list(original_idx=index,n=length(index)))
    attr(result,".data_ref")<-site;result
  }
  source_names<-setdiff(names(data),"t")
  folds<-list(target_folds=views(data$t,lapply(info,`[[`,"target_idx")),
    source_folds=setNames(lapply(seq_along(source_names),function(j)
      views(data[[source_names[j]]],lapply(info,function(x)x$source_idx[[j]]))),source_names))
  warnings<-list()
  answer<-tryCatch(withCallingHandlers(RoCE:::run_tate_crossfit(data,n_folds=5L,
    communication_mode="one_round",lambda_selection=1,verbose=FALSE,M_tau=5,M_tau_inference=5,
    n_cores=1L,nlambda_init=100L,family="binomial",use_lambda_cache=TRUE,
    precomputed_folds=folds,nuisance_lambda_rule="min",parallel_arms=FALSE),
    warning=function(w) {
      warnings[[length(warnings)+1L]]<<-list(message=conditionMessage(w),
        call=if(is.null(conditionCall(w)))NA_character_ else paste(deparse(conditionCall(w)),collapse=" "))
      invokeRestart("muffleWarning")
    }),error=identity)
  failed<-inherits(answer,"error")
  checks<-NULL;comparison<-NULL
  if(!failed) {
    phi<-answer$all_phi_tau;sizes<-vapply(data,`[[`,integer(1L),"n")
    centered_ss<-sum(vapply(split(phi,rep(seq_along(sizes),sizes)),function(x)sum((x-mean(x))^2),numeric(1L)))
    checks<-c(point_error=abs(mean(phi)-answer$estimate),variance_error=abs(centered_ss/length(phi)^2-answer$variance),
      target_only_error=abs(answer$target_only$estimate-old$target_only$estimate))
    comparison<-data.frame(version=c("v19","corrected"),estimate=c(old$estimate,answer$estimate),
      se=c(old$se,answer$se),target_only=c(old$target_only$estimate,answer$target_only$estimate),
      truth=artifact$tate_truth,inference_validated=FALSE)
  }
  roce_write_atomic_directory(args[3],function(stage) {
    saveRDS(list(fit=answer,checks=checks,warnings=warnings),file.path(stage,"tate_replay.rds"))
    if(!failed)write.csv(comparison,file.path(stage,"comparison.csv"),row.names=FALSE)
    writeLines(c(paste0("failed=",failed),paste0("package_fingerprint=",fingerprint),
      "sim_id=10013","rho=0","config=C1","data_and_outer_folds_reused=TRUE",
      "bootstrap_draws=0","new_mc_replications=0","inference_validated=FALSE",
      "warning_capture_scope=sequential_R_conditions","warning_capture_complete=FALSE",
      paste0("warning_count=",length(warnings)),paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/run_corrected_tate_replay.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="corrected full TATE replay")
  if(failed)stop(conditionMessage(answer))
  stopifnot(max(checks)<1e-10)
  print(comparison,row.names=FALSE);print(checks)
}
if(sys.nframe()==0L)main()
