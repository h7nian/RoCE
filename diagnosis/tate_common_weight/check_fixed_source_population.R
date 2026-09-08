#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: check_fixed_source_population.R LIBRARY PILOT_ROOT OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/source_score_derivative_rows.R")
  stopifnot(.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")
  holdout<-new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R",holdout)
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));expected<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  states<-list();coefficients<-list()
  for(config in c("C2","C3")) {
    probes<-checked(file.path(args[2],paste0("source_basis_",config,"_outer_1_audit_v1")),"validation_equation_states.rds")
    selected<-checked(file.path(args[2],paste0("source_basis_",config,"_stagewise_v1")),"selection_states.rds")
    for(arm in 0:1) {
      key<-paste(config,arm,sep=":");states[[key]]<-probes[[paste0("mu",arm)]]$state
      coefficients[[key]]<-selected$final_states[[paste0("mu",arm)]]$fit$coefficients
    }
  }
  # Check exact conditional averaging on the real score implementation.
  for(key in names(states)) for(site in c("target","source")) {
    state<-states[[key]];a<-coefficients[[key]];block<-state[[site]]
    block$W<-block$W[1:3,,drop=FALSE];block$Z<-block$Z[1:3,,drop=FALSE]
    block$A<-rep(.35,3);block$Y<-rep(.6,3)
    direct<-holdout$.projection_holdout_scores(state,block,site,a)$corrected
    direct_gradient<-.source_score_derivative_rows(state,block,site,a)
    averaged<-numeric(3);averaged_gradient<-matrix(0,3,length(a))
    for(A in 0:1)for(Y in 0:1) {
      weight<-dbinom(A,1,.35)*dbinom(Y,1,.6);b<-block;b$A<-rep(A,3);b$Y<-rep(Y,3)
      averaged<-averaged+weight*holdout$.projection_holdout_scores(state,b,site,a)$corrected
      averaged_gradient<-averaged_gradient+weight*.source_score_derivative_rows(state,b,site,a)
    }
    stopifnot(max(abs(direct-averaged))<1e-12,max(abs(direct_gradient-averaged_gradient))<1e-12)
  }
  ps<-RoCE:::get_face_ps_parameters(100L);outcome<-RoCE:::get_face_outcome_parameters(100L)
  calibration<-RoCE:::get_face_binary_calibration(100L,kappa=RoCE:::FACE_KAPPA)
  scores<-gradients<-list()
  for(batch in 1:25) {
    covariates<-RoCE:::with_seed(918731L+batch,RoCE:::generate_face_covariates(2000L,c(2000L,2000L),100L))
    for(site in c("target","source")) {
      label<-if(site=="target")"t" else "s1"
      X<-covariates$X[covariates$R==label,,drop=FALSE]
      centered<-sweep(X,2L,covariates$kappa,"-")
      propensity<-RoCE:::calculate_face_propensity(X,ps$alpha1,ps$alpha2)
      signal<-RoCE:::face_binary_logit(drop(centered%*%outcome$beta_linear+(X^2)%*%outcome$beta_squared),calibration)
      for(config in c("C2","C3"))for(arm in 0:1) {
        key<-paste(config,arm,sep=":");state<-states[[key]];a<-coefficients[[key]]
        truth<-plogis(signal+arm*RoCE:::FACE_BINARY_ATE_TARGET)
        block<-list(W=cbind(1,if(config=="C2")X else cbind(centered,X^2)),
          Z=cbind(1,if(config=="C3")X else cbind(X,X^2)),
          A=if(arm==1L)propensity else 1-propensity,Y=truth)
        score<-holdout$.projection_holdout_scores(state,block,site,a)
        gradient<-.source_score_derivative_rows(state,block,site,a)
        record_key<-paste(batch,site,key,sep=":")
        gradients[[record_key]]<-colMeans(gradient)
        scores[[record_key]]<-data.frame(batch,site,config,arm,n=nrow(X),
          original_mean=mean(score$original),corrected_mean=mean(score$corrected),conditional_truth=mean(truth))
      }
    }
  }
  scores<-do.call(rbind,scores)
  roce_write_atomic_directory(args[3],function(stage) {
    write.csv(scores,file.path(stage,"batch_scores.csv"),row.names=FALSE)
    saveRDS(gradients,file.path(stage,"batch_gradients.rds"))
    writeLines(c("conditional_expectation_identity_passed=TRUE","n_covariates_per_evaluated_site=50000",
      "nuisance_refits=0","bootstrap_draws=0","parameters_retuned=FALSE","inference_validated=FALSE",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_fixed_source_population.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="fixed source population evaluation")
  print(aggregate(cbind(original_mean,corrected_mean,conditional_truth)~site+config+arm,scores,mean))
}
if(sys.nframe()==0L)main()
