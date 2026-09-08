#!/usr/bin/env Rscript
# Dimension/derivative fixtures only: reducing coefficient vectors is not refitting C2/C3.
.reduce_source_basis_fixture <- function(state, basis, retained_columns) {
  stopifnot(basis %in% c("W", "Z"), is.numeric(retained_columns), !anyDuplicated(retained_columns))
  kind <- if (basis == "W") "alpha_" else "gamma_"
  values <- lapply(names(state$indices), function(name) {
    x <- state$theta[state$indices[[name]]]
    if (startsWith(name, kind)) x[retained_columns] else x
  })
  names(values) <- names(state$indices)
  reduce <- function(block) { block[[basis]] <- block[[basis]][, retained_columns, drop=FALSE]; block }
  state$target <- reduce(state$target); state$source <- reduce(state$source)
  state$pieces <- lapply(state$pieces, function(piece) lapply(piece, reduce))
  sizes <- lengths(values); ends <- cumsum(sizes)
  state$indices <- setNames(lapply(seq_along(sizes), function(k)
    seq.int(ends[k]-sizes[k]+1L,ends[k])),names(values))
  state$theta <- unlist(values,use.names=FALSE)
  state
}

main <- function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=3L) stop("usage: check_source_unequal_bases.R LIBRARY AUDITED_EQUATIONS OUTPUT")
  .libPaths(c(args[1],.libPaths())); suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R"); source("scripts/slurm/atomic_output.R")
  equations<-new.env(parent=globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R",equations)
  lines<-readLines(file.path(args[2],"sha256.txt"))
  expected<-substr(lines[substring(lines,67L)=="validation_equation_states.rds"],1L,64L)
  stopifnot(length(expected)==1L,roce_sha256_file(file.path(args[2],"validation_equation_states.rds"))==expected)
  probes<-readRDS(file.path(args[2],"validation_equation_states.rds"))
  rows<-list()
  for(arm in 0:1) for(basis in c("W","Z")) {
    state<-.reduce_source_basis_fixture(probes[[paste0("mu",arm)]]$state,basis,1:101)
    value<-equations$.evaluate_saved_projection(state,state$theta,TRUE)
    alpha<-state$theta[state$indices$alpha_final]; gamma<-state$theta[state$indices$gamma_final]
    native<-RoCE:::calculate_correction_term_cpp(state$source$Z[,-1,drop=FALSE],state$source$A,
      state$source$Y,gamma,alpha,state$source$W[,-1,drop=FALSE],5,1L,1L,1L)
    native_mean<-mean(RoCE:::predict_glm_cpp(state$target$W[,-1,drop=FALSE],alpha,1L,1L))+native$delta_ts
    stopifnot(abs(native_mean-value$score)<1e-12,!native$clip_diagnostics$any_safety_clipped)
    for(k in 1:6) {
      direction<-sin(seq_along(state$theta)*k);direction<-direction/sqrt(sum(direction^2));epsilon<-1e-6
      plus<-equations$.evaluate_saved_projection(state,state$theta+epsilon*direction)
      minus<-equations$.evaluate_saved_projection(state,state$theta-epsilon*direction)
      rows[[length(rows)+1L]]<-data.frame(arm,reduced_basis=basis,direction=k,
        outcome_dimension=ncol(state$target$W),site_dimension=ncol(state$target$Z),
        native_score_error=abs(native_mean-value$score),
        jacobian_error=max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(value$jacobian%*%direction))),
        score_gradient_error=abs((plus$score-minus$score)/(2*epsilon)-sum(value$score_gradient*direction)))
    }
  }
  rows<-do.call(rbind,rows)
  stopifnot(nrow(rows)==24L,max(rows$jacobian_error)<1e-7,max(rows$score_gradient_error)<1e-7)
  roce_write_atomic_directory(args[3],function(stage) {
    write.csv(rows,file.path(stage,"unequal_basis_checks.csv"),row.names=FALSE)
    writeLines(c("fixture_is_refitted_C2_C3=FALSE","inference_validated=FALSE",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/check_source_unequal_bases.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="unequal source basis fixture audit")
  print(aggregate(cbind(jacobian_error,score_gradient_error,native_score_error)~arm+reduced_basis,rows,max))
}
if(sys.nframe()==0L) main()
