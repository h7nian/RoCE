#!/usr/bin/env Rscript
main<-function(args=commandArgs(trailingOnly=TRUE)) {
  if(length(args)!=2L)stop("usage: summarize_fixed_source_population.R INPUT OUTPUT")
  source("scripts/slurm/result_provenance.R");source("scripts/slurm/atomic_output.R")
  manifest<-readLines(file.path(args[1],"sha256.txt"))
  for(name in c("batch_scores.csv","batch_gradients.rds","metadata.txt")) {
    expected<-substr(manifest[substring(manifest,67L)==name],1L,64L)
    stopifnot(length(expected)==1L,roce_sha256_file(file.path(args[1],name))==expected)
  }
  scores<-read.csv(file.path(args[1],"batch_scores.csv"))
  gradients<-readRDS(file.path(args[1],"batch_gradients.rds"))
  keys<-with(scores,paste(batch,site,config,arm,sep=":"))
  expected<-expand.grid(batch=1:25,site=c("target","source"),config=c("C2","C3"),arm=0:1,stringsAsFactors=FALSE)
  expected_keys<-with(expected,paste(batch,site,config,arm,sep=":"))
  stopifnot(nrow(scores)==200L,!anyDuplicated(keys),setequal(keys,expected_keys),
    setequal(names(gradients),expected_keys),all(scores$n==2000L),
    all(is.finite(as.matrix(scores[c("original_mean","corrected_mean","conditional_truth")]))))
  batches<-summaries<-gradient_rows<-list()
  mcse<-function(x)sd(x)/sqrt(length(x))
  for(config in c("C2","C3")) {
    values<-do.call(rbind,lapply(1:25,function(batch) {
      x<-scores[scores$config==config & scores$batch==batch,]
      sign<-ifelse(x$arm==1,1,-1)
      truth<-sum(sign[x$site=="target"]*x$conditional_truth[x$site=="target"])
      original<-sum(sign*x$original_mean);corrected<-sum(sign*x$corrected_mean)
      data.frame(config,batch,truth,original,corrected,original_error=original-truth,
        corrected_error=corrected-truth,change=corrected-original)
    }))
    batches[[config]]<-values
    summaries[[config]]<-data.frame(config,n_integration_batches=25L,
      truth=mean(values$truth),truth_integration_mcse=mcse(values$truth),
      original_expected_tate=mean(values$original),corrected_expected_tate=mean(values$corrected),
      original_fixed_fit_error=mean(values$original_error),original_error_mcse=mcse(values$original_error),
      corrected_fixed_fit_error=mean(values$corrected_error),corrected_error_mcse=mcse(values$corrected_error),
      paired_change=mean(values$change),paired_change_mcse=mcse(values$change))
    for(arm in 0:1) {
      matrix<-do.call(rbind,lapply(1:25,function(batch)
        gradients[[paste(batch,"target",config,arm,sep=":")]]+
        gradients[[paste(batch,"source",config,arm,sep=":")]]))
      stopifnot(all(is.finite(matrix)))
      means<-colMeans(matrix);errors<-apply(matrix,2,mcse)
      gradient_rows[[paste(config,arm)]]<-data.frame(config,arm,coordinate=seq_along(means),
        expected_gradient=means,integration_mcse=errors)
    }
  }
  summaries<-do.call(rbind,summaries)
  roce_write_atomic_directory(args[2],function(stage) {
    write.csv(summaries,file.path(stage,"fixed_fit_summary.csv"),row.names=FALSE)
    write.csv(do.call(rbind,batches),file.path(stage,"paired_batch_errors.csv"),row.names=FALSE)
    write.csv(do.call(rbind,gradient_rows),file.path(stage,"gradient_integration.csv"),row.names=FALSE)
    writeLines(c("inference_validated=FALSE","uncertainty_is_integration_mcse_not_estimator_se=TRUE",
      "new_training_replications=0","parameters_retuned=FALSE",
      paste0("input_manifest_sha256=",roce_sha256_file(file.path(args[1],"sha256.txt"))),
      paste0("script_sha256=",roce_sha256_file("diagnosis/tate_common_weight/summarize_fixed_source_population.R"))),file.path(stage,"metadata.txt"))
    files<-list.files(stage,full.names=TRUE)
    writeLines(paste(vapply(files,roce_sha256_file,""),basename(files),sep="  "),file.path(stage,"sha256.txt"))
  },caller="fixed source population summary")
  print(summaries,row.names=FALSE)
}
if(sys.nframe()==0L)main()
