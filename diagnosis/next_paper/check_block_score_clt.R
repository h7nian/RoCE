#!/usr/bin/env Rscript
# Check covariance factorization/whitening on fitted packets and radial density
# bounds numerically. These checks do not replace the cited CLT or its moments.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: repository_root fitted_study_root new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
study <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
source(file.path(repository,"diagnosis/next_paper/dispersion_bias_intervals.R"))
source(file.path(repository,"diagnosis/next_paper/joint_dispersion_intervals.R"))
dir.create(output,recursive=TRUE)
root_matrix <- function(matrix,inverse=FALSE) {
  spectral <- eigen(matrix,symmetric=TRUE)
  scale <- max(abs(spectral$values))
  stopifnot(min(spectral$values)>=-1e-10*scale)
  if(inverse) stopifnot(min(spectral$values)>0)
  values <- if(inverse) 1/sqrt(spectral$values) else sqrt(pmax(0,spectral$values))
  sweep(spectral$vectors,2L,values,"*")%*%t(spectral$vectors)
}
manifest <- read.csv(file.path(study,"manifest.csv"),stringsAsFactors=FALSE)
rows <- list()
for(index in seq_len(nrow(manifest))) {
  task <- manifest[index,]
  directory <- file.path(study,"cases",task$task_id)
  complete <- readRDS(file.path(directory,"COMPLETE.rds"))
  stopifnot(identical(unname(complete$hashes),unname(tools::md5sum(file.path(directory,names(complete$hashes))))))
  packets <- readRDS(file.path(directory,"population_packets.rds"))$variants
  for(construction in names(packets)) {
    population <- packets[[construction]]$population
    covariance <- population$covariance
    width <- task$K+1L
    target_loading <- matrix(0,2L*width,4L)
    target_loading[1L,1L] <- target_loading[width+1L,2L] <- 1
    target_loading[2:width,3L] <- target_loading[width+(2:width),4L] <- 1
    loading <- matrix(0,2L*width,4L+2L*task$K)
    loading[,1:4] <- target_loading%*%root_matrix(population$target_covariance)
    for(site in seq_len(task$K)) {
      positions <- c(site+1L,width+site+1L)
      columns <- 4L+c(2L*site-1L,2L*site)
      loading[positions,columns] <- root_matrix(population$source_moments[[site]]$covariance)
    }
    covariance_error <- max(abs(tcrossprod(loading)-covariance))
    calibration <- make_joint_dispersion_calibration(covariance,task$K%/%2L,task$K,"shared_prediction_bound")
    map <- calibration$bias_map
    contrast_covariance <- map%*%covariance%*%t(map)
    whitening <- root_matrix(contrast_covariance,inverse=TRUE)
    transform <- whitening%*%map%*%loading
    whitening_error <- max(abs(tcrossprod(transform)-diag(nrow(map))))
    scalar <- drop(crossprod(calibration$coefficients,loading))/sqrt(calibration$variance)
    scalar_error <- abs(sum(scalar^2)-1)
    stopifnot(covariance_error<1e-12,whitening_error<1e-9,scalar_error<1e-9)
    # Hypothetical component perturbations test the relative-covariance
    # argument; they are not an estimator or a confidence region.
    factors <- c(rep(1.08,4L),rep(1+.1*sin(seq_len(task$K)),each=2L))
    altered_loading <- sweep(loading,2L,sqrt(factors),"*")
    estimated <- tcrossprod(altered_loading)
    estimated_inverse_root <- root_matrix(estimated,inverse=TRUE)
    relative <- estimated_inverse_root%*%covariance%*%estimated_inverse_root
    epsilon <- max(abs(eigen(relative,symmetric=TRUE,only.values=TRUE)$values-1))
    altered <- make_joint_dispersion_calibration(estimated,task$K%/%2L,task$K,"shared_prediction_bound")
    altered_map <- altered$bias_map
    altered_whitening <- root_matrix(altered_map%*%estimated%*%t(altered_map),inverse=TRUE)
    normalized_loading <- altered_whitening%*%altered_map%*%loading
    omega <- tcrossprod(normalized_loading)
    spectrum <- eigen(omega,symmetric=TRUE,only.values=TRUE)$values
    coupling <- sum((sqrt(spectrum)-1)^2)
    kappa <- max(sqrt(1+epsilon)-1,1-sqrt(1-epsilon))
    upper <- nrow(altered_map)*kappa^2
    stopifnot(epsilon<1,min(spectrum)>=1-epsilon-1e-10,
      max(spectrum)<=1+epsilon+1e-10,coupling<=upper+1e-10)
    rows[[length(rows)+1L]] <- cbind(task,data.frame(construction=construction,
      covariance_error=covariance_error,whitening_error=whitening_error,scalar_error=scalar_error,
      relative_covariance_error=epsilon,gaussian_coupling_error=coupling,coupling_upper=upper))
  }
}
density_rows <- list()
for(degrees in c(1L,2L,4L,16L,64L,256L)) for(noncentrality in c(0,.1,4,16,128,1024)) {
  density <- function(radius) {
    if(radius==0) return(if(degrees==1L) 2*dnorm(sqrt(noncentrality)) else 0)
    2*radius*dchisq(radius^2,df=degrees,ncp=noncentrality)
  }
  optimized <- optimize(density,c(0,sqrt(degrees+noncentrality)+10),maximum=TRUE)
  maximum <- max(density(0),optimized$objective)
  stopifnot(maximum<=1+1e-10)
  density_rows[[length(density_rows)+1L]] <- data.frame(degrees=degrees,
    noncentrality=noncentrality,maximum_density=maximum)
}
write.csv(do.call(rbind,rows),file.path(output,"fitted_whitening_checks.csv"),row.names=FALSE)
write.csv(do.call(rbind,density_rows),file.path(output,"radial_density_checks.csv"),row.names=FALSE)
writeLines(c("BLOCK_SCORE_CLT_GEOMETRY_CHECKS_PASSED",paste("Fitted covariance packets:",length(rows)),
  "36 radial-density numerical checks; covariance tensorization and perturbation identities only.",
  "No empirical moment-rate verification or confidence-interval coverage study."),file.path(output,"checks.txt"))
