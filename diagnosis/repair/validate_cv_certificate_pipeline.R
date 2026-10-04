#!/usr/bin/env Rscript
# One frozen-library pipeline comparison arm; final comparisons use paired data.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=6L) stop("Usage: library output role scenario p protocol")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
role <- match.arg(arguments[3L],c("reference","off","on"))
scenario <- match.arg(arguments[4L],c("C2","C3"))
dimension <- suppressWarnings(as.integer(arguments[5L]))
protocol <- match.arg(arguments[6L],c("one_round","two_round"))
if(is.na(dimension) || !dimension%in%c(10L,200L) ||
   !startsWith(output,"/scratch.global/zhan9381/FACE-HD/")) stop("Invalid pipeline validation case")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
dir.create(output,recursive=TRUE,showWarnings=FALSE)
configuration <- list(library=library_path,role=role,scenario=scenario,p=dimension,
  protocol=protocol,n_per_site=1000L,K=1L,n_folds=10L,nlambda=100L,seed=71301L)
script_path <- sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE))
stopifnot(length(script_path)==1L)
configuration$workflow_hash <- digest::digest(file=script_path,algo="sha256")
package_root <- find.package("RoCE")
code_files <- c(file.path(package_root,"R/RoCE.rdb"),
  list.files(file.path(package_root,"libs"),pattern="[.](so|dll)$",recursive=TRUE,full.names=TRUE))
stopifnot(length(code_files)>=2L)
configuration$package_hashes <- tools::md5sum(code_files)
configuration_path <- file.path(output,"configuration.rds")
if(file.exists(configuration_path)) {
  stopifnot(identical(readRDS(configuration_path),configuration))
} else {
  temporary <- paste0(configuration_path,".tmp.",Sys.getpid())
  saveRDS(configuration,temporary)
  stopifnot(file.rename(temporary,configuration_path))
}
set.seed(configuration$seed)
data <- split_data_by_site(generate_bounded_data(n_target=1000L,n_source_sizes=1000L,
  p=dimension,config=scenario))
data_hash <- digest::digest(data,algo="sha256")
model_path <- file.path(output,"fitted.rds")
if(file.exists(model_path)) {
  saved <- readRDS(model_path)
  stopifnot(identical(saved$configuration,configuration),identical(saved$data_hash,data_hash))
  fitted <- saved$fitted
} else {
  fit_arguments <- list(data_split=data,n_folds=10L,communication_mode=protocol,
    crossfit_layers=3L,nlambda_init=100L,M_tau=12,M_tau_inference=12,
    target_nuisance_method="hou_calibrated",source_validation_method="calibrated",
    calibration_layout="compact",nuisance_tol=1e-10,nuisance_solver="proximal_newton",
    calibration_control=list(recipe="score_derivative",target_propensity_initialization="calibrated",target_radius=12),
    n_cores=2L,parallel_arms=TRUE,verbose=FALSE,checkpoint_dir=file.path(output,"checkpoints"))
  if(role!="reference") fit_arguments$nuisance_cv_certificate <- role=="on"
  elapsed <- system.time(fitted <- do.call(run_tate_crossfit,fit_arguments))[["elapsed"]]
  temporary <- paste0(model_path,".tmp.",Sys.getpid())
  saveRDS(list(fitted=fitted,configuration=configuration,data_hash=data_hash,
    elapsed_this_attempt=elapsed,slurm_restarts=Sys.getenv("SLURM_RESTART_COUNT","0")),temporary)
  stopifnot(file.rename(temporary,model_path))
}
collect_nuisances <- function(object,path="fit") {
  if(is.numeric(object) && !is.null(attr(object,"lambda_used"))) {
    return(setNames(list(list(coefficients=as.numeric(object),lambda=attr(object,"lambda_used"),
      lambda_min=attr(object,"lambda_min"),lambda_1se=attr(object,"lambda_1se"),
      cv_seed=attr(object,"cv_seed"))),path))
  }
  if(!is.list(object)) return(list())
  labels <- names(object)
  if(is.null(labels)) labels <- as.character(seq_along(object))
  missing_label <- is.na(labels) | !nzchar(labels)
  labels[missing_label] <- paste0("[",which(missing_label),"]")
  result <- list()
  for(index in seq_along(object)) {
    result <- c(result,collect_nuisances(object[[index]],paste(path,labels[index],sep="/")))
  }
  result
}
nuisances <- collect_nuisances(fitted)
stopifnot(length(nuisances)>0L,!anyDuplicated(names(nuisances)))
summaries <- lapply(c("common_tate","separate_arms","joint_tate"),function(mode) {
  result <- reaggregate_tate_crossfit(data,fitted,aggregation_mode=mode,M_tau_inference=12)
  list(mode=mode,estimate=result$estimate,se=result$se,weights=result$fold_weights,
    arm_estimates=vapply(fitted$arm_results,`[[`,numeric(1L),"estimate"),
    arm_standard_errors=vapply(fitted$arm_results,`[[`,numeric(1L),"se"))
})
certified <- 0L
for(path in list.files(file.path(output,"checkpoints"),pattern="[.]rds$",recursive=TRUE,full.names=TRUE)) {
  entry <- readRDS(path)
  count <- attr(entry$fit,"cv_certified_nonconvergent_fold_fits")
  if(!is.null(count)) certified <- certified+count
}
if(role!="reference") stopifnot(identical(fitted$nuisance_cv_certificate,role=="on"))
result <- list(configuration=configuration,data_hash=data_hash,nuisances=nuisances,
  summaries=summaries,certified_fits=certified)
temporary <- file.path(output,paste0("comparison.rds.tmp.",Sys.getpid()))
saveRDS(result,temporary)
stopifnot(file.rename(temporary,file.path(output,"comparison.rds")))
writeLines(paste("PIPELINE_ARM_COMPLETE",role,scenario,dimension,protocol),file.path(output,"COMPLETE"))
cat("PIPELINE_ARM_COMPLETE",role,scenario,dimension,protocol,"certified",certified,"\n")
