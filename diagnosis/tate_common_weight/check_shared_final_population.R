#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: check_shared_final_population.R LIBRARY PILOT_ROOT OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/shared_calibration_basis.R")
  source("diagnosis/tate_common_weight/calibrated_site_score.R")
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));expected<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  shared<-original<-list()
  for(config in c("C2","C3")) {
    shared[[config]]<-checked(file.path(args[2],paste0("shared_final_",config,"_v1")),"shared_final_states.rds")
    original[[config]]<-checked(file.path(args[2],paste0("source_basis_",config,"_outer_1_v1")),"source_validation_states.rds")
    audit<-read.csv(file.path(args[2],paste0("shared_final_",config,"_audit_v1"),"calibration_checks.csv"))
    stopifnot(nrow(audit)==2L,max(audit$kkt_error)<1e-4,max(audit$point_error)<1e-12)
  }
  ps<-RoCE:::get_face_ps_parameters(100L);outcome<-RoCE:::get_face_outcome_parameters(100L)
  calibration<-RoCE:::get_face_binary_calibration(100L,config="C1",kappa=RoCE:::FACE_KAPPA)
  scores<-gradients<-clipping<-list()
  for(batch in 1:25) {
    covariates<-RoCE:::with_seed(928731L+batch,RoCE:::generate_face_covariates(2000L,c(2000L,2000L),100L))
    for(site in c("target","source")) {
      X<-covariates$X[covariates$R==if(site=="target")"t" else "s1",,drop=FALSE]
      centered<-sweep(X,2L,covariates$kappa,"-")
      propensity<-RoCE:::calculate_face_propensity(X,ps$alpha1,ps$alpha2)
      signal<-RoCE:::face_binary_logit(drop(centered%*%outcome$beta_linear+(X^2)%*%outcome$beta_squared),calibration)
      for(config in c("C2","C3"))for(arm in 0:1) {
        W<-if(config=="C2")X else cbind(centered,X^2)
        Z<-if(config=="C3")X else cbind(X,X^2)
        U<-cbind(1,.apply_shared_calibration_basis(shared[[config]]$basis_spec,W,Z))
        truth<-plogis(signal+arm*RoCE:::FACE_BINARY_ATE_TARGET)
        A<-if(arm==1L)propensity else 1-propensity
        old<-original[[config]]$results[[paste0("mu",arm)]]$result
        new<-shared[[config]]$results[[paste0("mu",arm)]]$final_fit
        baseline<-.calibrated_site_score(cbind(1,W),cbind(1,Z),A,truth,old$alpha_ts,old$gamma_s,site)
        candidate<-.calibrated_site_score(U,U,A,truth,new$alpha,new$gamma,site)
        key<-paste(batch,site,config,arm,sep=":")
        gradients[[key]]<-colMeans(candidate$gradient)
        scores[[key]]<-data.frame(batch,site,config,arm,n=nrow(X),original_mean=mean(baseline$score),
          corrected_mean=mean(candidate$score),conditional_truth=mean(truth))
        clipping[[key]]<-data.frame(batch,site,config,arm,n=nrow(X),
          original_clipped=baseline$n_clipped,shared_final_clipped=candidate$n_clipped)
      }
    }
  }
  roce_write_atomic_directory(args[3],function(stage) {
    write.csv(do.call(rbind,scores),file.path(stage,"batch_scores.csv"),row.names=FALSE)
    saveRDS(gradients,file.path(stage,"batch_gradients.rds"))
    write.csv(do.call(rbind,clipping),file.path(stage,"clipping_counts.csv"),row.names=FALSE)
    writeLines(c("inference_validated=FALSE","corrected_method=shared_final_calibration",
      "n_covariates_per_evaluated_site=50000","seed_base=928731","parameters_retuned=FALSE",
      "initial_or_final_nuisance_refits=0","bootstrap_draws=0",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_shared_final_population.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="shared final calibration independent evaluation")
}
if(sys.nframe()==0L)main()
