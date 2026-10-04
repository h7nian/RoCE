#!/usr/bin/env Rscript
# Reconstruct full/ablated estimates and exchange eta on paired saved scores.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: repository completed_calibration_review new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
review_path <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
checks <- jsonlite::fromJSON(file.path(review_path,"pairing_checks.json"),simplifyVector=FALSE)
if(checks$component!="source" || checks$pairs!=600L || checks$crossfit_layers!=2L) stop("Requires complete source-calibration MC200 pairing")
for(path in names(checks$hashes)) if(!identical(digest::digest(file=path,algo="sha256"),checks$hashes[[path]])) stop("Paired input provenance changed: ",path)
library_path <- checks$configurations[[1L]]$library
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/repair/calibration_score_diagnostics.R"))
roots <- unlist(checks$roots,use.names=FALSE)
names(roots) <- c("full","ablation")
manifests <- lapply(roots,function(root) read.csv(file.path(root,"manifest.csv"),stringsAsFactors=FALSE))
select_tasks <- function(manifest) {
  selected <- subset(manifest,p==checks$dimension & K==2 & rho==0 & config %in% c("C1","C2","C3"))
  if(nrow(selected)!=600L || anyDuplicated(paste(selected$config,selected$sim_id)) ||
     any(vapply(split(selected$sim_id,selected$config),function(seeds)
       !identical(sort(as.integer(seeds)),1:200),logical(1L)))) stop("Incomplete MC200 task index")
  selected <- selected[order(selected$config,selected$sim_id),,drop=FALSE]
  rownames(selected) <- NULL
  selected
}
manifests <- lapply(manifests,select_tasks)
stopifnot(identical(manifests$full[c("config","p","K","rho","sim_id","protocol","deviation_mechanism")],
                    manifests$ablation[c("config","p","K","rho","sim_id","protocol","deviation_mechanism")]))
reference <- read.csv(file.path(review_path,"paired_estimates.csv"),stringsAsFactors=FALSE)
reference <- reference[reference$aggregation_mode=="joint_tate",,drop=FALSE]
dir.create(output,recursive=TRUE)
dir.create(file.path(output,"workflow"))
for(name in c("review_calibration_scores.R","calibration_score_diagnostics.R")) {
  if(!file.copy(file.path(repository,"diagnosis/repair",name),file.path(output,"workflow",name))) stop("Could not archive workflow")
}
jsonlite::write_json(list(pairing_review=review_path,library=library_path,
  input_review_sha256=digest::digest(file=file.path(review_path,"paired_estimates.csv"),algo="sha256"),
  workflow_sha256=lapply(list.files(file.path(output,"workflow"),full.names=TRUE),function(path)
    list(file=path,sha256=digest::digest(file=path,algo="sha256")))),file.path(output,"provenance.json"),pretty=TRUE)
rows <- vector("list",600L)
for(index in seq_len(600L)) {
  task <- manifests$full[index,]
  packets <- lapply(names(roots),function(variant) {
    directory <- file.path(roots[[variant]],"tasks",manifests[[variant]]$task_id[index])
    if(!file.exists(file.path(directory,"COMPLETE")) ||
       readLines(file.path(directory,"status.txt"),warn=FALSE)!="COMPLETE") stop("Uncommitted source score packet")
    packet <- readRDS(file.path(directory,"score_derivative_layer_records.rds"))
    packet$data_hash <- readLines(file.path(directory,"data_sha256.txt"),warn=FALSE)
    if(packet$crossfit_layers!=2L || packet$source_validation_method!="outer_fit") stop("Unexpected score training boundary")
    expected <- reference[reference$config==task$config & reference$sim_id==task$sim_id & reference$variant==variant,]
    if(nrow(expected)!=3L || !all(expected$data_sha256==packet$data_hash)) stop("Review/packet data hashes differ")
    packet
  })
  names(packets) <- names(roots)
  validate_paired_score_packets(packets$full,packets$ablation)
  learning <- lapply(packets,fit_packet_joint_weights,aggregation_cutoff=2)
  own <- lapply(names(roots),function(variant) evaluate_packet_weights(
    packets[[variant]],learning[[variant]]$weights,learning[[variant]]))
  names(own) <- names(roots)
  for(variant in names(roots)) {
    expected <- reference[reference$config==task$config & reference$sim_id==task$sim_id & reference$variant==variant,]
    expected <- expected[match(c("mu1","mu0","tate"),expected$quantity),]
    value <- own[[variant]]
    if(max(abs(as.numeric(value$estimates)-expected$estimate))>1e-12 ||
       max(abs(c(diag(value$covariance),value$tate_variance)-expected$variance))>1e-12 ||
       max(abs(c(diag(value$covariance_fixed_weights),value$tate_variance_fixed_weights)-expected$variance_fixed_weights))>1e-12) {
      stop("Saved scores failed to reconstruct the committed estimate/variance")
    }
  }
  evaluated <- c(own,list(
    ablation_full_eta=evaluate_packet_weights(packets$ablation,learning$full$weights,learning$full),
    full_ablation_eta=evaluate_packet_weights(packets$full,learning$ablation$weights,learning$ablation)))
  values <- lapply(names(evaluated),function(variant) {
    value <- evaluated[[variant]]
    data.frame(variant=variant,candidate="aggregate",quantity=c("mu1","mu0","tate"),
      estimate=as.numeric(value$estimates),truth=c(packets$full$arm_truth,
        packets$full$arm_truth[["mu1"]]-packets$full$arm_truth[["mu0"]]),
      variance=c(diag(value$covariance),value$tate_variance),
      variance_fixed_weights=c(diag(value$covariance_fixed_weights),value$tate_variance_fixed_weights))
  })
  for(variant in names(roots)) {
    candidates <- summarize_packet_candidates(packets[[variant]])
    values[[length(values)+1L]] <- data.frame(variant=variant,
      candidate=candidates$candidate,quantity=candidates$quantity,estimate=candidates$estimate,
      truth=candidates$truth,variance=candidates$score_variance,
      variance_fixed_weights=candidates$score_variance)
  }
  rows[[index]] <- cbind(config=task$config,p=task$p,K=task$K,rho=task$rho,sim_id=task$sim_id,
                         data_sha256=packets$full$data_hash,do.call(rbind,values))
  if(index %% 20L==0L) {
    writeLines(sprintf("Validated %d/600 paired score packets",index),file.path(output,"progress.txt"))
    cat(sprintf("Validated %d/600 paired score packets\n",index));flush.console()
  }
}
write.csv(do.call(rbind,rows),file.path(output,"score_estimates.csv"),row.names=FALSE)
writeLines(c("CALIBRATION_SCORE_MC200_CHECKS_PASSED",
  "All600 paired generated-data hashes and actual outer folds agree.",
  "All1200 reconstructed own-weight fits reproduce saved mu1/mu0/TATE estimates and both variance formulas within1e-12.",
  "Eta training scores cover only the corresponding outer-training observations.",
  "Weight exchange is a diagnostic using both fitted programs; it is not a replacement estimator or a full nuisance-refit variance."),file.path(output,"checks.txt"))
file.create(file.path(output,"COMPLETE"))
