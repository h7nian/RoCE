#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L)stop("usage: check_source_score_derivative_rows.R LIBRARY SEED PILOT_ROOT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("diagnosis/tate_common_weight/source_working_basis.R")
  source("diagnosis/tate_common_weight/source_projection_validation_state.R")
  source("diagnosis/tate_common_weight/source_score_derivative_rows.R")
  equations<-new.env(parent=globalenv());holdout<-new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R",equations)
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R",holdout)
  artifact<-readRDS(file.path(args[2],"artifacts.rds"))$group_result$artifacts[["0"]]
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info[[1L]]
  for(config in c("C2","C3")) {
    data<-.source_working_basis(artifact$data_split,config)
    probes<-readRDS(file.path(args[3],paste0("source_basis_",config,"_outer_1_audit_v1"),"validation_equation_states.rds"))
    selected<-readRDS(file.path(args[3],paste0("source_basis_",config,"_stagewise_v1"),"selection_states.rds"))
    for(arm in 0:1) {
      state<-probes[[paste0("mu",arm)]]$state;a<-selected$final_states[[paste0("mu",arm)]]$fit$coefficients
      blocks<-list(target=equations$.saved_projection_block(data$t,info$target_idx,arm),
        source=equations$.saved_projection_block(data$s1,info$source_idx[[1L]],arm))
      rows<-lapply(names(blocks),function(site).source_score_derivative_rows(state,blocks[[site]],site,a));names(rows)<-names(blocks)
      value<-equations$.evaluate_saved_projection(.source_projection_validation_state(state,blocks$target,blocks$source),state$theta,TRUE)
      mean_error<-max(abs(colMeans(rows$target)+colMeans(rows$source)-value$score_gradient+drop(crossprod(value$jacobian,a))))
      error<-0
      for(k in 1:6) for(site in names(blocks)) {
        direction<-sin(seq_along(a)*k);direction<-direction/sqrt(sum(direction^2));epsilon<-1e-6
        plus<-minus<-state;plus$theta<-state$theta+epsilon*direction;minus$theta<-state$theta-epsilon*direction
        numeric_rows<-(holdout$.projection_holdout_scores(plus,blocks[[site]],site,a)$corrected-
          holdout$.projection_holdout_scores(minus,blocks[[site]],site,a)$corrected)/(2*epsilon)
        error<-max(error,abs(numeric_rows-drop(rows[[site]]%*%direction)))
      }
      stopifnot(mean_error<1e-10,error<1e-7)
      noise<-sqrt(apply(rows$target,2,var)/nrow(rows$target)+apply(rows$source,2,var)/nrow(rows$source))
      gradient<-colMeans(rows$target)+colMeans(rows$source)
      usable<-noise>0
      cat(config,"arm",arm,"mean identity",mean_error,"row derivative error",error,
        "max gradient/noise (descriptive only)",max(abs(gradient[usable])/noise[usable]),"\n")
    }
  }
}
if(sys.nframe()==0L)main()
