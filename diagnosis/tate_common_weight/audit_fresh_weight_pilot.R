#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: audit_fresh_weight_pilot.R LIBRARY INPUT OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  for(file in c("quadratic_bias_weights.R","quadratic_weight_influence.R","balanced_fold_gradient.R"))
    source(file.path("diagnosis/tate_common_weight",file))
  manifest<-readLines(file.path(args[2],"sha256.txt"))
  for(name in c("artifacts.rds","results.csv","metadata.txt")) {
    expected<-substr(manifest[substring(manifest,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(args[2],name))==expected)
  }
  bundle<-readRDS(file.path(args[2],"artifacts.rds"));rows<-read.csv(file.path(args[2],"results.csv"))
  stopifnot(is.null(bundle$failure),nrow(rows)==24L,length(unique(rows$sim_id))==1L,
    unique(rows$sim_id) %in% 20001:20005,!any(rows$inference_validated))
  methods<-c("target_only","original_common_weights","quadratic_fixed_weights","quadratic_analytic_weights")
  rhos<-c(0,.5,1,1.5,2,2.5)
  expected<-expand.grid(rho=rhos,method=methods,stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(paste(rows$rho,rows$method)),
    setequal(paste(rows$rho,rows$method),paste(expected$rho,expected$method)),
    identical(names(bundle$grouped$artifacts),as.character(rhos)))
  checks<-list()
  for(rho in rhos) {
    artifact<-bundle$grouped$artifacts[[as.character(rho)]];fit<-artifact$direct_tate_results$one_round_crossfit
    data<-artifact$data_split;stored<-bundle$gradients[[as.character(rho)]];q<-stored$weight_layer
    stopifnot(identical(names(data),c("t","s1","s2")))
    sizes<-vapply(data,`[[`,integer(1L),"n");phi<-fit$all_phi_tau
    variance<-sum(vapply(split(phi,rep(seq_along(sizes),sizes)),function(x)sum((x-mean(x))^2),numeric(1L)))/sum(sizes)^2
    balanced<-lapply(seq_along(data),function(site) {
      labels<-integer(sizes[site]);counts<-integer(sizes[site])
      for(k in 1:5) {
        f<-fit$intermediates$fold_info[[k]];ids<-if(site==1L)f$target_idx else f$source_idx[[site-1L]]
        labels[ids]<-k;counts<-counts+tabulate(ids,nbins=sizes[site])
      }
      stopifnot(all(counts==1L))
      .balanced_arm_fold_gradient(q$gradient[[site]],data[[site]]$A,labels)
    })
    vq<-sum(vapply(balanced,`[[`,numeric(1L),"variance"))
    covariance_error<-abs(q$weight_linearized_variance-q$fixed_variance-q$indirect_variance-q$direct_indirect_cross_term)
    masses<-lapply(q$gradient,function(x)rep(1,length(x)));derivative_error<-0
    for(k in 1:2) {
      direction<-lapply(seq_along(masses),function(s)sin(seq_along(masses[[s]])*k+s));eps<-1e-5
      numerical<-(.quadratic_weight_functional(fit,Map(function(m,d)m+eps*d,masses,direction))-
        .quadratic_weight_functional(fit,Map(function(m,d)m-eps*d,masses,direction)))/(2*eps)
      derivative_error<-max(derivative_error,abs(numerical-sum(unlist(q$gradient)*unlist(direction))))
    }
    x<-rows[rows$rho==rho,];x<-x[match(methods,x$method),]
    estimates<-c(fit$target_only$estimate,fit$estimate,q$estimate,q$estimate)
    standard_errors<-c(fit$target_only$se,fit$se,sqrt(q$fixed_variance),sqrt(vq))
    errors<-c(point=abs(mean(phi)-fit$estimate),variance=abs(variance-fit$variance),
      covariance=covariance_error,rows=max(abs(x$estimate-estimates),abs(x$se-standard_errors),
        abs(x$truth-artifact$tate_truth),abs(x$bias-(estimates-artifact$tate_truth))),
      ci=max(abs(x$ci_lower-(estimates-qnorm(.975)*standard_errors)),
        abs(x$ci_upper-(estimates+qnorm(.975)*standard_errors))))
    stopifnot(max(errors)<1e-9,derivative_error<1e-8,
      all(x$coverage==(x$ci_lower<=x$truth & x$ci_upper>=x$truth)))
    diagnostics<-fit$nuisance_fit_diagnostics
    hard<-grep("nonconverged|degenerate",names(diagnostics),value=TRUE)
    stopifnot(all(as.matrix(diagnostics[hard])==0))
    checks[[as.character(rho)]]<-data.frame(rho,maximum_algebra_error=max(errors),
      maximum_weight_derivative_error=derivative_error,
      cv_invalid_fold_fits=sum(as.matrix(diagnostics[grep("cv_invalid_fold_fits",names(diagnostics),value=TRUE)])))
  }
  checks<-do.call(rbind,checks)
  roce_write_atomic_directory(args[3],function(stage) {
    write.csv(checks,file.path(stage,"numerical_checks.csv"),row.names=FALSE)
    writeLines(c("numerical_gate=passed","inference_validated=FALSE","coverage_gate_passed=FALSE",
      paste0("sim_id=",unique(rows$sim_id)),paste0("input_manifest_sha256=",roce_sha256_file(file.path(args[2],"sha256.txt"))),
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/audit_fresh_weight_pilot.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="fresh analytic-weight pilot numerical audit")
  print(checks,row.names=FALSE)
}
if(sys.nframe()==0L)main()
