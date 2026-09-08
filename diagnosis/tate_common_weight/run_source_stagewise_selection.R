#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(4L,5L)) stop("usage: run_source_stagewise_selection.R LIBRARY SEED PILOT_ROOT OUTPUT [CONFIG]")
  working_config<-if(length(args)==5L) args[5] else "C1"
  .libPaths(c(args[1], .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  stopifnot(.weight_calibration_installed_package_fingerprint(args[1],roce_sha256_file)==
    "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(args[2],"sha256.txt"))==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  files <- c("sparse_moment_projection.R", "source_final_projection.R", "source_joint_projection.R",
    "source_projection_validation_state.R", "target_projection_validation.R", "source_stagewise_selection.R",
    "source_working_basis.R")
  for (file in files) source(file.path("diagnosis/tate_common_weight", file))
  equations <- new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", equations)
  holdout <- new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R", holdout)
  code_files <- file.path("diagnosis/tate_common_weight", c(files, "run_source_stagewise_selection.R",
    "probe_saved_nuisance_projection.R", "audit_saved_projection_holdout.R", "SOURCE_STAGEWISE_SELECTION_PROTOCOL.md"))
  code_hashes <- vapply(code_files,roce_sha256_file,"")
  checked <- function(root,name) {
    manifest<-readLines(file.path(root,"sha256.txt"))
    expected<-substr(manifest[substring(manifest,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(root,name))==expected)
    readRDS(file.path(root,name))
  }
  bundle<-checked(args[2],"artifacts.rds")
  artifact<-bundle$group_result$artifacts[["0"]]
  artifact$data_split<-.source_working_basis(artifact$data_split,working_config)
  audit_prefix<-if(working_config=="C1") "source_validation" else paste0("source_basis_",working_config)
  outer_probes<-NULL; outer_manifest<-NA_character_
  if(working_config!="C1") {
    outer_root<-file.path(args[3],paste0(audit_prefix,"_outer_1_audit_v1"))
    outer_probes<-checked(outer_root,"validation_equation_states.rds")
    outer_manifest<-roce_sha256_file(file.path(outer_root,"sha256.txt"))
  }
  info<-artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  candidates<-c(null=NA_real_,c0.5=.5,c1=1,c2=2)
  rows<-attempts<-list(); manifests<-character()
  catch_projection<-function(expr) tryCatch(expr,error=function(e) {
    if (!inherits(e,"roce_projection_error")) stop(e)
    e
  })
  block_risk<-function(value,index,a,name) {
    pos<-index[[name]]; other<-setdiff(seq_along(a),pos)
    sign<-if(grepl("^alpha_",name)) -1 else 1
    H<-sign*value$jacobian[pos,pos,drop=FALSE]
    rhs<-sign*(value$score_gradient[pos]-drop(crossprod(value$jacobian[other,pos,drop=FALSE],a[other])))
    .5*sum(a[pos]*drop(H%*%a[pos]))-sum(rhs*a[pos])
  }
  for(fold in 2:5) {
    root<-file.path(args[3],paste0(audit_prefix,"_split_1_",fold,"_audit_v1"))
    probes<-checked(root,"validation_equation_states.rds")
    manifests[as.character(fold)]<-roce_sha256_file(file.path(root,"sha256.txt"))
    for(arm in 0:1) {
      probe<-probes[[paste0("mu",arm)]]; state<-probe$state
      stopifnot(ncol(state$target$W)==ncol(artifact$data_split$t$W_outcome)+1L,
        ncol(state$target$Z)==ncol(artifact$data_split$t$Z_site)+1L)
      target<-equations$.saved_projection_block(artifact$data_split$t,info[[fold]]$target_idx,arm)
      source<-equations$.saved_projection_block(artifact$data_split$s1,info[[fold]]$source_idx[[1L]],arm)
      validation<- .source_projection_validation_state(state,target,source)
      value<-equations$.evaluate_saved_projection(validation,state$theta,TRUE)
      for(final_candidate in names(candidates)) {
        final<-if(final_candidate=="null") NULL else catch_projection(
          .fit_source_final_projection(state,probe$evaluated,candidates[[final_candidate]]))
        final_failed<-inherits(final,"error")
        a_final<-numeric(length(state$theta))
        if(!final_failed && !is.null(final)) for(name in names(final))
          a_final[state$indices[[name]]]<-final[[name]]$coefficients
        final_risk<-if(final_failed) NA_real_ else sum(vapply(c("alpha_final","gamma_final"),
          function(name) block_risk(value,state$indices,a_final,name),numeric(1L)))
        for(initial_candidate in names(candidates)) {
          key<-paste(fold,arm,final_candidate,initial_candidate,sep=":")
          attempt<-if(final_failed) final else if(final_candidate=="null")
            list(coefficients=a_final,null_projection=TRUE) else catch_projection(
              .fit_source_joint_projection(state,probe$evaluated,candidates[[final_candidate]],
                if(initial_candidate=="null") NULL else candidates[[initial_candidate]]))
          initial_failed<-!final_failed && inherits(attempt,"error")
          initial_risk<-if(final_failed || initial_failed) NA_real_ else {
            final_index<-unlist(state$indices[c("alpha_final","gamma_final")],use.names=FALSE)
            stopifnot(identical(attempt$coefficients[final_index],a_final[final_index]))
            sum(vapply(setdiff(names(state$indices),c("alpha_final","gamma_final")),
              function(name) block_risk(value,state$indices,attempt$coefficients,name),numeric(1L)))
          }
          attempts[[key]]<-attempt
          rows[[length(rows)+1L]]<-data.frame(final_candidate,initial_candidate,arm,fold,
            n_validation=nrow(target$W)+nrow(source$W),final_risk,initial_risk,final_failed,initial_failed)
        }
      }
    }
  }
  losses<-do.call(rbind,rows)
  selection<-.select_source_stagewise_candidates(losses)
  final_choice<-selection$final_selection$selected
  initial_choice<-selection$initial_selection$selected
  # No outer evaluation data were used in selection. Refit on outer training.
  final_states<-list(); final_reports<-list()
  for(arm in 0:1) {
    state<-if(is.null(outer_probes)) equations$.saved_projection_state(bundle,arm) else
      outer_probes[[paste0("mu",arm)]]$state
    value<-equations$.evaluate_saved_projection(state,state$theta,TRUE)
    if(!is.null(outer_probes)) stopifnot(identical(value,outer_probes[[paste0("mu",arm)]]$evaluated),
      ncol(state$target$W)==ncol(artifact$data_split$t$W_outcome)+1L,
      ncol(state$target$Z)==ncol(artifact$data_split$t$Z_site)+1L)
    fitted<-if(final_choice=="null") list(coefficients=numeric(length(state$theta)),null_projection=TRUE) else
      catch_projection(.fit_source_joint_projection(state,value,candidates[[final_choice]],
        if(initial_choice=="null") NULL else candidates[[initial_choice]]))
    failed<-inherits(fitted,"error")
    scores<-if(failed) NULL else lapply(c("target","source"),function(site) {
      data<-artifact$data_split[[if(site=="target") "t" else "s1"]]
      ids<-if(site=="target") info[[1L]]$target_idx else info[[1L]]$source_idx[[1L]]
      block<-equations$.saved_projection_block(data,ids,arm)
      stopifnot(!any(ids %in% state[[site]]$index))
      holdout$.projection_holdout_scores(state,block,site,fitted$coefficients)
    })
    final_states[[paste0("mu",arm)]]<-list(fit=fitted,scores=scores)
    final_reports[[paste0("mu",arm)]]<-data.frame(arm,failed,
      original_mean=if(failed) NA_real_ else sum(vapply(scores,function(x)mean(x$original),numeric(1L))),
      corrected_mean=if(failed) NA_real_ else sum(vapply(scores,function(x)mean(x$corrected),numeric(1L))))
  }
  final_reports<-do.call(rbind,final_reports)
  stopifnot(identical(code_hashes,vapply(code_files,roce_sha256_file,"")))
  roce_write_atomic_directory(args[4],function(stage) {
    write.csv(losses,file.path(stage,"candidate_losses.csv"),row.names=FALSE)
    write.csv(final_reports,file.path(stage,"outer_fold_report.csv"),row.names=FALSE)
    saveRDS(list(selection=selection,attempts=attempts,final_states=final_states),file.path(stage,"selection_states.rds"))
    write.csv(data.frame(path=code_files,sha256=code_hashes),file.path(stage,"code_hashes.csv"),row.names=FALSE)
    writeLines(c("inference_validated=FALSE","common_weights_recomputed=FALSE","nuisance_refits=0",
      paste0("working_config=",working_config),paste0("outer_training_manifest_sha256=",outer_manifest),
      paste0("final_candidate=",final_choice),paste0("initial_candidate=",initial_choice),
      paste0("validation_",names(manifests),"_manifest_sha256=",manifests)),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="source stagewise selection")
  print(selection$final_selection$scores); print(selection$initial_selection$scores); print(final_reports)
  if(any(final_reports$failed)) stop("selected outer-training projection failed; retained without switching candidates")
}
if(sys.nframe()==0L) main()
