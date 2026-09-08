#!/usr/bin/env Rscript
# Exact Gaussian-lasso fixture for the shared outcome CV normalization path.
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=2L)stop("usage: check_outcome_cv_scale.R LIBRARY OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  # Keep each +/- pair in the same fold so every training design is centered.
  # This removes intercept/slope iteration error from the normalization test.
  x<-c(as.vector(rbind(-(1:30)/30,(1:30)/30)),seq(-1,1,length.out=40))
  X<-matrix(x,ncol=1);A<-rep(c(1,0),c(60,40))
  Y<-1+2*x+.3*sin(3*x);arm_ids<-which(A==1);fold_id<-rep(rep(1:5,length.out=30),each=2)
  grid<-c(.01,.02,.04);q<-mean(A==1)
  closed_fit<-function(ids,scale,lambda) {
    xc<-x[ids]-mean(x[ids]);yc<-Y[ids]-mean(Y[ids])
    covariance<-scale*mean(xc*yc);variance<-scale*mean(xc^2)
    slope<-sign(covariance)*max(abs(covariance)-lambda,0)/variance
    c(mean(Y[ids])-slope*mean(x[ids]),slope)
  }
  native<-RoCE:::select_lambda_cv_calibrated_outcome_cpp(X,Y,A,c(0,0),grid,5L,10000L,1e-10,
    1L,5,X,0L,0L,fold_id)
  reconstructed<-standardized<-numeric(length(grid));rows<-list()
  for(i in seq_along(grid))for(fold in 1:5) {
    train<-arm_ids[fold_id!=fold];validation<-arm_ids[fold_id==fold]
    fraction<-length(train)/length(arm_ids)
    current<-closed_fit(train,length(train)/length(x),grid[i])
    comparable<-closed_fit(train,q,grid[i])
    equivalent<-closed_fit(train,q,grid[i]/fraction)
    stopifnot(max(abs(current-equivalent))<1e-12)
    loss<-function(beta)sum(.5*(Y[validation]-beta[1]-beta[2]*x[validation])^2)/length(x)
    reconstructed[i]<-reconstructed[i]+loss(current)/5
    standardized[i]<-standardized[i]+loss(comparable)/5
    rows[[length(rows)+1L]]<-data.frame(lambda=grid[i],fold,train_fraction=fraction,
      cv_effective_lambda=grid[i]/fraction,final_refit_lambda=grid[i],
      current_slope=current[2],training_mean_scale_slope=comparable[2])
  }
  error<-max(abs(native$cv_scores-reconstructed))
  print(data.frame(lambda=grid,native=as.numeric(native$cv_scores),reconstructed,standardized))
  cat("observed CV score discrepancy",error,"\n")
  stopifnot(error<1e-7,max(abs(reconstructed-standardized))>1e-5)
  full<-RoCE:::fit_unified_outcome_cpp(X,Y,A,c(0,0),0L,0L,grid[1],10000L,1e-10,1L,TRUE,5,X,numeric())
  stopifnot(max(abs(full$alpha-closed_fit(arm_ids,q,grid[1])))<1e-7)
  scores<-data.frame(lambda=grid,native_cv_score=as.numeric(native$cv_scores),
    reconstructed_cv_score=reconstructed,training_mean_scale_cv_score=standardized)
  roce_write_atomic_directory(args[2],function(stage) {
    write.csv(scores,file.path(stage,"cv_score_reconstruction.csv"),row.names=FALSE)
    write.csv(do.call(rbind,rows),file.path(stage,"fold_penalty_scales.csv"),row.names=FALSE)
    writeLines(c("production_code_changed=FALSE","coverage_explanation_established=FALSE",
      paste0("native_cv_score_error=",error),paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_outcome_cv_scale.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="outcome CV objective scale audit")
  print(scores);cat("native CV reconstruction error",error,"\n")
}
if(sys.nframe()==0L)main()
