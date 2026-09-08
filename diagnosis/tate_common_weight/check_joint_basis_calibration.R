#!/usr/bin/env Rscript
# Structural population check only; no production basis change or SE scaling.
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=2L)stop("usage: check_joint_basis_calibration.R LIBRARY OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  results<-list()
  for(case in c("correct_outcome","correct_tilting")) {
    x<-if(case=="correct_outcome")c(-1,0,1) else c(-1,-.3,.4,1)
    U<-cbind(1,x,x^2);linear<-U[,1:2,drop=FALSE]
    target_mass<-rep(1/length(x),length(x))
    counts<-if(case=="correct_outcome")c(90L,120L,90L)else c(80L,120L,120L,80L)
    treated_counts<-if(case=="correct_outcome")c(24L,84L,36L)else rep(50L,4L)
    source_mass<-treated_counts/sum(counts)
    true_mean<-plogis(-.2+.7*x+1.2*x^2+if(case=="correct_tilting") .8*x^3 else 0)
    W0<-if(case=="correct_outcome")U else linear
    Z0<-if(case=="correct_outcome")linear else U
    alpha0<-optim(rep(0,ncol(W0)),function(a)sum(target_mass*(log1p(exp(drop(W0%*%a)))-true_mean*drop(W0%*%a))),
      function(a)drop(crossprod(W0,target_mass*(plogis(drop(W0%*%a))-true_mean))),method="BFGS",
      control=list(reltol=1e-14,maxit=10000))$par
    gamma0<-optim(rep(0,ncol(Z0)),function(g)sum(target_mass*drop(Z0%*%g)+source_mass*exp(-drop(Z0%*%g))),
      function(g)drop(crossprod(Z0,target_mass-source_mass*exp(-drop(Z0%*%g)))),method="BFGS",
      control=list(reltol=1e-14,maxit=10000))$par
    source_x<-rep(x,counts);source_y<-rep(true_mean,counts)
    source_a<-unlist(lapply(seq_along(x),function(k)c(rep(1,treated_counts[k]),rep(0,counts[k]-treated_counts[k]))))
    source_U<-cbind(source_x,source_x^2)
    source_W0<-if(ncol(W0)==3L)source_U else matrix(source_x,ncol=1)
    source_Z0<-if(ncol(Z0)==3L)source_U else matrix(source_x,ncol=1)
    d0<-plogis(drop(W0%*%alpha0));d0<-d0*(1-d0)
    for(repair in c(FALSE,TRUE)) {
      W<-if(repair)U else W0;Z<-if(repair)U else Z0
      source_W<-if(repair)source_U else source_W0
      source_Z<-if(repair)source_U else source_Z0
      alpha<-RoCE:::fit_unified_outcome_cpp(source_W,source_y,source_a,gamma0,
        1L,1L,0,10000L,1e-12,1L,TRUE,5,source_Z0,numeric())$alpha
      gamma<-RoCE:::fit_unified_density_ratio_cpp(source_Z,source_a,drop(crossprod(Z,target_mass*d0)),
        alpha0,0,10000L,1e-12,TRUE,5,source_W0,1L,1L,1L,numeric())$gamma
      prediction<-plogis(drop(W%*%alpha));weight<-exp(-drop(Z%*%gamma))
      stopifnot(max(abs(drop(source_Z0%*%gamma0[-1]+gamma0[1])))<5,
        max(abs(drop(W0%*%alpha0)))<5,max(abs(drop(Z%*%gamma)))<5)
      da<-drop(crossprod(W,(target_mass-source_mass*weight)*prediction*(1-prediction)))
      dg<- -drop(crossprod(Z,source_mass*weight*(true_mean-prediction)))
      point<-sum(target_mass*prediction)+sum(source_mass*weight*(true_mean-prediction))
      error<-point-sum(target_mass*true_mean)
      if(repair)stopifnot(max(abs(c(da,dg)))<1e-6,abs(error)<1e-7)
      else stopifnot(max(abs(c(da,dg)))>1e-3)
      results[[length(results)+1L]]<-data.frame(case,shared_final_basis=repair,
        outcome_gradient=max(abs(da)),tilting_gradient=max(abs(dg)),population_point_error=error,
        outcome_model_error=max(abs(prediction-true_mean)))
    }
  }
  results<-do.call(rbind,results)
  roce_write_atomic_directory(args[2],function(stage) {
    write.csv(results,file.path(stage,"population_checks.csv"),row.names=FALSE)
    writeLines(c("production_method_changed=FALSE","high_dimensional_inference_validated=FALSE",
      "penalties=0","truncation_inactive=TRUE",paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_joint_basis_calibration.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="shared final basis population check")
  print(results,row.names=FALSE)
}
if(sys.nframe()==0L)main()
