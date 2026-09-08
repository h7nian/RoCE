#!/usr/bin/env Rscript
# Local fixed-coefficient score derivatives; not a full fitted-estimator IF.
main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=4L) stop("usage: audit_selected_source_sensitivity.R LIBRARY SEED PILOT_ROOT OUTPUT")
  .libPaths(c(args[1],.libPaths())); suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R"); source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  source("diagnosis/tate_common_weight/source_projection_validation_state.R")
  equations<-new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R",equations)
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt")); expected<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  artifact<-checked(args[2],"artifacts.rds")$group_result$artifacts[["0"]]
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info[[1L]]
  blocks<-directions<-list()
  for(config in c("C2","C3")) {
    data<-.source_working_basis(artifact$data_split,config)
    probes<-checked(file.path(args[3],paste0("source_basis_",config,"_outer_1_audit_v1")),"validation_equation_states.rds")
    selected<-checked(file.path(args[3],paste0("source_basis_",config,"_stagewise_v1")),"selection_states.rds")
    for(arm in 0:1) {
      state<-probes[[paste0("mu",arm)]]$state
      a<-selected$final_states[[paste0("mu",arm)]]$fit$coefficients
      target<-equations$.saved_projection_block(data$t,info$target_idx,arm)
      source<-equations$.saved_projection_block(data$s1,info$source_idx[[1L]],arm)
      heldout<-.source_projection_validation_state(state,target,source)
      for(region in c("training","outer_evaluation")) {
        current<-if(region=="training")state else heldout
        value<-equations$.evaluate_saved_projection(current,state$theta,TRUE)
        adjusted_gradient<-value$score_gradient-drop(crossprod(value$jacobian,a))
        for(name in names(state$indices)) {
          pos<-state$indices[[name]]
          blocks[[length(blocks)+1L]]<-data.frame(config,arm,region,block=name,
            coefficient_l1=sum(abs(a[pos])), original_max_abs_gradient=max(abs(value$score_gradient[pos])),
            corrected_max_abs_gradient=max(abs(adjusted_gradient[pos])))
        }
        for(k in 1:6) {
          direction<-sin(seq_along(a)*k); direction<-direction/sqrt(sum(direction^2));epsilon<-1e-6
          plus<-equations$.evaluate_saved_projection(current,state$theta+epsilon*direction)
          minus<-equations$.evaluate_saved_projection(current,state$theta-epsilon*direction)
          numeric_gradient<-((plus$score-sum(plus$moments*a))-(minus$score-sum(minus$moments*a)))/(2*epsilon)
          error<-abs(numeric_gradient-sum(adjusted_gradient*direction))
          stopifnot(error<1e-7)
          directions[[length(directions)+1L]]<-data.frame(config,arm,region,direction=k,error)
        }
      }
    }
  }
  blocks<-do.call(rbind,blocks); directions<-do.call(rbind,directions)
  stopifnot(nrow(blocks)==80L,nrow(directions)==48L)
  roce_write_atomic_directory(args[4],function(stage) {
    write.csv(blocks,file.path(stage,"block_sensitivity.csv"),row.names=FALSE)
    write.csv(directions,file.path(stage,"directional_checks.csv"),row.names=FALSE)
    writeLines(c("inference_validated=FALSE","projection_coefficients_held_fixed=TRUE",
      "outer_evaluation_used_for_selection=FALSE",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/audit_selected_source_sensitivity.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="selected source local sensitivity")
  print(aggregate(cbind(original_max_abs_gradient,corrected_max_abs_gradient)~config+arm+region,blocks,max))
  print(blocks[blocks$config=="C3" & blocks$region=="training",])
}
if(sys.nframe()==0L) main()
