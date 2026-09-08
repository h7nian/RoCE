#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=4L)stop("usage: decompose_shared_calibration_gap.R LIBRARY SEED PILOT_ROOT OUTPUT")
  .libPaths(c(args[1],.libPaths()));suppressPackageStartupMessages(library(RoCE,lib.loc=args[1]))
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/shared_calibration_basis.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  loss<-new.env(parent=globalenv());sys.source("diagnosis/tate_common_weight/audit_shared_final_calibration.R",loss)
  checked<-function(root,name) {
    lines<-readLines(file.path(root,"sha256.txt"));h<-substr(lines[substring(lines,67L)==name],1L,64L)
    stopifnot(length(h)==1L,roce_sha256_file(file.path(root,name))==h);readRDS(file.path(root,name))
  }
  artifact<-checked(args[2],"artifacts.rds")$group_result$artifacts[["0"]]
  rows<-summaries<-list();clip<-function(x)pmax(-5,pmin(5,x))
  for(config in c("C2","C3")) {
    saved<-checked(file.path(args[3],paste0("shared_final_",config,"_v1")),"shared_final_states.rds")
    data<-.source_working_basis(artifact$data_split,config)
    Ut<-cbind(1,.apply_shared_calibration_basis(saved$basis_spec,data$t$W_outcome,data$t$Z_site))
    Us<-cbind(1,.apply_shared_calibration_basis(saved$basis_spec,data$s1$W_outcome,data$s1$Z_site))
    ti<-unlist(saved$target_ids[2:5],use.names=FALSE);si<-unlist(saved$source_ids[2:5],use.names=FALSE)
    for(arm in 0:1) {
      fitted<-saved$results[[paste0("mu",arm)]];alpha<-fitted$final_fit$alpha;gamma<-fitted$final_fit$gamma
      mt<-plogis(drop(Ut[ti,,drop=FALSE]%*%alpha));ms<-plogis(drop(Us[si,,drop=FALSE]%*%alpha))
      ht<-mt*(1-mt);hs<-ms*(1-ms)
      linear<-drop(Us[si,,drop=FALSE]%*%gamma);wr<-exp(-linear);wf<-exp(-clip(linear))
      A<-as.numeric(data$s1$A[si]==arm);residual<-data$s1$Y[si]-ms
      h0t<-h0s<-w0<-numeric();target_summaries<-list()
      for(k in 2:5) {
        key<-paste0("k2_",k);t<-saved$target_ids[[k]];s<-saved$source_ids[[k]]
        pt<-plogis(drop(cbind(1,data$t$W_outcome[t,,drop=FALSE])%*%fitted$initial_alpha[[key]]))
        ps<-plogis(clip(drop(cbind(1,data$s1$W_outcome[s,,drop=FALSE])%*%fitted$initial_alpha[[key]])))
        h0t<-c(h0t,pt*(1-pt));h0s<-c(h0s,ps*(1-ps))
        w0<-c(w0,exp(-clip(drop(cbind(1,data$s1$Z_site[s,,drop=FALSE])%*%fitted$initial_gamma[[key]]))))
        target_summaries[[key]]<-colMeans(Ut[t,,drop=FALSE]*(pt*(1-pt)))
      }
      gt<-function(x)colMeans(Ut[ti,,drop=FALSE]*x)
      gs<-function(x)colMeans(Us[si,,drop=FALSE]*x)
      calibration<-loss$.shared_final_losses(Us[si,,drop=FALSE],A,data$s1$Y[si],w0,h0s,
        Reduce(`+`,target_summaries)/4,alpha,gamma)
      alpha_gradient<-gt(ht)-gs(A*wf*hs)
      gamma_gradient<- -gs(A*wf*residual*(abs(linear)<5))
      alpha_terms<-list(calibration_residual=calibration$tilting_gradient,
        target_derivative_change=gt(ht-h0t),target_fold_normalization=gt(h0t)-Reduce(`+`,target_summaries)/4,
        source_derivative_change= -gs(A*wf*(hs-h0s)),inference_weight_clipping= -gs(A*(wf-wr)*h0s))
      gamma_terms<-list(calibration_residual=calibration$outcome_gradient,
        initial_final_weight_change= -gs(A*(wf-w0)*residual),
        inference_derivative_clipping=gs(A*wf*residual*(abs(linear)>=5)))
      for(block in c("alpha_final","gamma_final")) {
        terms<-if(block=="alpha_final")alpha_terms else gamma_terms
        gradient<-if(block=="alpha_final")alpha_gradient else gamma_gradient
        error<-max(abs(Reduce(`+`,terms)-gradient));stopifnot(error<1e-12)
        for(term in names(terms)) {
          rows[[length(rows)+1L]]<-data.frame(config,arm,block,term,coordinate=seq_along(gradient),value=terms[[term]])
          summaries[[length(summaries)+1L]]<-data.frame(config,arm,block,term,
            maximum_component=max(abs(terms[[term]])),maximum_gradient=max(abs(gradient)),identity_error=error)
        }
      }
    }
  }
  summaries<-do.call(rbind,summaries)
  roce_write_atomic_directory(args[4],function(stage) {
    write.csv(do.call(rbind,rows),file.path(stage,"coordinate_decomposition.csv"),row.names=FALSE)
    write.csv(summaries,file.path(stage,"component_summary.csv"),row.names=FALSE)
    writeLines(c("inference_validated=FALSE","parameters_changed=FALSE",
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/decompose_shared_calibration_gap.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="shared calibration gradient gap decomposition")
  print(summaries[summaries$block=="gamma_final",],row.names=FALSE)
}
if(sys.nframe()==0L)main()
