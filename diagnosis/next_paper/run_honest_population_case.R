#!/usr/bin/env Rscript
# One fitted-data diagnostic repeat; no confidence interval is asserted valid.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=8L) stop("Usage: checked_R_library workflow output scenario K seed rho invalid_sources")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
workflow <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
scenario <- match.arg(arguments[4L],c("C1","C2","C3"))
numbers <- suppressWarnings(as.numeric(arguments[5:8]))
if(any(!is.finite(numbers)) || any(numbers[c(1L,2L,4L)]!=floor(numbers[c(1L,2L,4L)])) ||
   numbers[1L]<1 || numbers[1L]>.Machine$integer.max || numbers[2L]<1 || numbers[2L]>.Machine$integer.max ||
   numbers[3L]<0 || numbers[4L]<0 || numbers[4L]>numbers[1L]) stop("Invalid diagnostic design")
source_count <- as.integer(numbers[1L]); seed <- as.integer(numbers[2L])
rho <- numbers[3L]; invalid_sources <- as.integer(numbers[4L])
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
files <- c("candidate_score_summary.R","honest_candidate_scores.R","honest_population_moments.R",
           "dispersion_bias_intervals.R","joint_dispersion_intervals.R","run_honest_population_case.R")
paths <- file.path(workflow,files)
for(path in paths[-length(paths)]) source(path)
configuration <- list(scenario=scenario,K=source_count,seed=seed,rho=rho,
  invalid_sources=invalid_sources,n_per_site=1000L,p=4L,calibration_folds=3L,
  nlambda=100L,library=library_path,source_hashes=tools::md5sum(paths))
configuration_hash <- digest::digest(configuration,algo="sha256")
dir.create(output,recursive=TRUE,showWarnings=FALSE)
atomic_rds <- function(object,name) {
  path <- file.path(output,name)
  temporary <- paste0(path,".tmp.",Sys.getpid())
  saveRDS(object,temporary)
  if(!file.rename(temporary,path)) stop("Could not commit ",name)
}
configuration_path <- file.path(output,"configuration.rds")
if(file.exists(configuration_path)) {
  if(!identical(readRDS(configuration_path),configuration)) stop("Diagnostic inputs or code changed")
} else atomic_rds(configuration,"configuration.rds")
complete_path <- file.path(output,"COMPLETE.rds")
if(file.exists(complete_path)) {
  completed <- readRDS(complete_path)
  stopifnot(identical(completed$configuration_hash,configuration_hash),
    identical(unname(completed$hashes),unname(tools::md5sum(file.path(output,names(completed$hashes))))))
  quit(status=0L)
}
set.seed(seed)
dgp <- generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,source_count),p=4L,
  config=scenario,ate_deviation=rho,n_deviated_sites=invalid_sources,deviation_mechanism="treated_arm")
data <- split_data_by_site(dgp)
data_hash <- digest::digest(data,algo="sha256")
fitted_path <- file.path(output,"fitted.rds")
if(file.exists(fitted_path)) {
  saved <- readRDS(fitted_path)
  stopifnot(identical(saved$data_hash,data_hash),identical(saved$configuration_hash,configuration_hash))
  fitted <- saved$fitted
} else {
  fitted <- fit_honest_candidates(data,calibration_folds=3L,nlambda=100L,
    checkpoint_dir=file.path(output,"checkpoints"))
  atomic_rds(list(fitted=fitted,data_hash=data_hash,configuration_hash=configuration_hash),"fitted.rds")
}
packet <- evaluate_honest_candidates(fitted,data)
sources <- setdiff(fitted$sites,"t"); arms <- c("mu1","mu0")
shifts <- matrix(0,source_count,2L,dimnames=list(sources,arms))
if(invalid_sources>0L) shifts[seq_len(invalid_sources),"mu1"] <- rho
common_fitted <- fitted
for(arm in arms) for(site in sources) {
  common_fitted$models[[arm]]$source[[site]]$outcome <- fitted$models[[arm]]$common$coefficients
}
common_packet <- evaluate_honest_candidates(common_fitted,data)
stopifnot(max(abs(common_packet$training_shifts))==0)
variants <- list(translated_source=list(fitted=fitted,packet=packet),
  common_outcome=list(fitted=common_fitted,packet=common_packet))
rows <- list(); aggregate_rows <- list(); outputs <- list()
for(construction in names(variants)) {
  variant <- variants[[construction]]
  population <- honest_population_moments(variant$fitted,variant$packet,dgp,24L,shifts)
  coarse <- honest_population_moments(variant$fitted,variant$packet,dgp,16L,shifts)
  quadrature_error <- max(abs(population$means-coarse$means),abs(population$covariance-coarse$covariance))
  stopifnot(quadrature_error<1e-8)
  width <- nrow(population$means)
  valid <- rbind(c(TRUE,TRUE),shifts==0)
  dimnames(valid) <- dimnames(population$means)
  for(site in seq_len(width)) for(quantity in c("mu1","mu0","tate")) {
    contrast <- numeric(2L*width)
    if(quantity!="mu0") contrast[site] <- 1
    if(quantity!="mu1") contrast[width+site] <- if(quantity=="tate") -1 else 1
    truth <- if(quantity=="mu1") population$truth[1L] else if(quantity=="mu0")
      population$truth[2L] else population$truth[1L]-population$truth[2L]
    expected <- sum(contrast*as.numeric(population$means))
    conditional_variance <- drop(crossprod(contrast,population$covariance%*%contrast))
    estimated_variance <- drop(crossprod(contrast,variant$packet$covariance%*%contrast))
    is_valid <- if(quantity=="tate") all(valid[site,]) else valid[site,quantity]
    rows[[length(rows)+1L]] <- data.frame(construction=construction,site=rownames(population$means)[site],
      quantity=quantity,valid=is_valid,conditional_mean=expected,truth=as.numeric(truth),
      conditional_bias=expected-truth,conditional_sd=sqrt(conditional_variance),
      bias_over_conditional_sd=(expected-truth)/sqrt(conditional_variance),
      estimated_to_conditional_variance=estimated_variance/conditional_variance)
  }
  design <- cbind(mu1=rep(c(1,0),each=width),mu0=rep(c(0,1),each=width))
  selected <- which(as.numeric(valid)==1)
  for(covariance_type in c("conditional_population","estimated")) {
    covariance <- if(covariance_type=="estimated") variant$packet$covariance else population$covariance
    oracle <- joint_mean_gls(covariance,design,selected)
    expected <- sum(oracle$coefficients*as.numeric(population$means))
    variance <- drop(crossprod(oracle$coefficients,population$covariance%*%oracle$coefficients))
    truth <- population$truth[1L]-population$truth[2L]
    aggregate_rows[[length(aggregate_rows)+1L]] <- data.frame(construction=construction,
      covariance_type=covariance_type,fixed_weight_bias=expected-truth,
      fixed_weight_sd=sqrt(variance),bias_over_fixed_weight_sd=(expected-truth)/sqrt(variance),
      reported_to_fixed_weight_variance=oracle$variance/variance,
      quadrature_error=quadrature_error)
  }
  outputs[[construction]] <- list(packet=variant$packet,population=population,valid=valid)
}
atomic_rds(list(configuration_hash=configuration_hash,data_hash=data_hash,variants=outputs),"population_packets.rds")
for(item in list(list(name="candidate_moments.csv",rows=rows),list(name="oracle_valid_aggregate.csv",rows=aggregate_rows))) {
  path <- file.path(output,item$name);temporary <- paste0(path,".tmp.",Sys.getpid())
  write.csv(do.call(rbind,item$rows),temporary,row.names=FALSE)
  stopifnot(file.rename(temporary,path))
}
artifact_names <- c("population_packets.rds","candidate_moments.csv","oracle_valid_aggregate.csv")
hashes <- tools::md5sum(file.path(output,artifact_names));names(hashes) <- artifact_names
atomic_rds(list(configuration_hash=configuration_hash,hashes=hashes,
  scope=paste("Conditional fitted-score diagnostic; validity labels and population moments are oracle inputs, not a feasible interval.",
    "Estimated-covariance GLS summaries hold its realized weights fixed; they are not conditional moments of a data-dependent weighted estimator.")),"COMPLETE.rds")
cat("HONEST_POPULATION_CASE_COMPLETE",scenario,source_count,seed,rho,"\n")
