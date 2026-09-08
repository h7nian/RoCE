#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: check_corrected_weight_layer.R SEED_BUNDLE REPLAY_BUNDLE OUTPUT")
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  for(file in c("quadratic_bias_weights.R","quadratic_weight_influence.R","balanced_fold_gradient.R"))
    source(file.path("diagnosis/tate_common_weight",file))
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));h<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(h)==1L,roce_sha256_file(file.path(root,name))==h);readRDS(file.path(root,name))
  }
  data<-checked(args[1],"artifacts.rds")$group_result$artifacts[["0"]]$data_split
  fit<-checked(args[2],"tate_replay.rds")$fit
  result<-.quadratic_weight_influence(fit)
  masses<-lapply(result$gradient,function(x)rep(1,length(x)))
  checks<-do.call(rbind,lapply(1:8,function(k) {
    direction<-lapply(seq_along(masses),function(s)sin(seq_along(masses[[s]])*k+s));epsilon<-1e-5
    plus<-Map(function(x,d)x+epsilon*d,masses,direction);minus<-Map(function(x,d)x-epsilon*d,masses,direction)
    numerical<-(.quadratic_weight_functional(fit,plus)-.quadratic_weight_functional(fit,minus))/(2*epsilon)
    data.frame(direction=k,error=abs(numerical-sum(unlist(result$gradient)*unlist(direction))),
      fixed_weight_error=abs(numerical-sum(unlist(result$direct)*unlist(direction))))
  }))
  stopifnot(max(checks$error)<1e-8,abs(.quadratic_weight_functional(fit)-result$estimate)<1e-12)
  balanced<-lapply(seq_along(data),function(site) {
    folds<-integer(data[[site]]$n)
    for(k in seq_along(fit$intermediates$fold_info)) {
      info<-fit$intermediates$fold_info[[k]]
      ids<-if(site==1L)info$target_idx else info$source_idx[[site-1L]]
      folds[ids]<-k
    }
    .balanced_arm_fold_gradient(result$gradient[[site]],data[[site]]$A,folds)
  })
  summary<-data.frame(original_estimate=fit$estimate,quadratic_estimate=result$estimate,
    fixed_weight_se=sqrt(result$fixed_variance),weight_layer_se=sqrt(result$weight_linearized_variance),
    balanced_weight_layer_se=sqrt(sum(vapply(balanced,`[[`,numeric(1L),"variance"))),
    indirect_variance=result$indirect_variance,cross_term=result$direct_indirect_cross_term,
    maximum_derivative_error=max(checks$error),inference_validated=FALSE)
  roce_write_atomic_directory(args[3],function(stage) {
    write.csv(summary,file.path(stage,"weight_layer_summary.csv"),row.names=FALSE)
    write.csv(checks,file.path(stage,"directional_checks.csv"),row.names=FALSE)
    saveRDS(list(weight_layer=result,balanced=balanced),file.path(stage,"weight_gradients.rds"))
    writeLines(c("nuisance_fits_differentiated=FALSE","inference_validated=FALSE","resampling_draws=0",
      paste0("replay_manifest_sha256=",roce_sha256_file(file.path(args[2],"sha256.txt"))),
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_corrected_weight_layer.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="corrected replay analytic weight layer")
  print(summary,row.names=FALSE)
}
if(sys.nframe()==0L)main()
