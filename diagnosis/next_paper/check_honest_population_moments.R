#!/usr/bin/env Rscript
# Independent outcome-law checks for the conditional population diagnostic.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: checked_R_library repository_root new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
repository <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/next_paper/honest_population_moments.R"))
dir.create(output,recursive=TRUE)
errors <- list()
for(mechanism in c("treated_arm","both_arms")) {
  set.seed(59301L)
  dgp <- generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,4L),p=4L,
    config="C1",ate_deviation=.5,n_deviated_sites=2L,deviation_mechanism=mechanism)
  sources <- paste0("s",1:4); arms <- c("mu1","mu0")
  shifts <- matrix(0,4L,2L,dimnames=list(sources,arms))
  shifts[1:2,"mu1"] <- .5
  if(mechanism=="both_arms") shifts[1:2,"mu0"] <- .5
  models <- lapply(c(mu1=1L,mu0=0L),function(value) {
    alpha <- if(value==1L) dgp$alpha1_true else dgp$alpha0_true
    list(target=list(outcome=alpha,weight=c(0,(2*value-1)*c(.35,-.25,.125,-.0625))),
      common=list(coefficients=alpha),source=setNames(lapply(sources,function(site)
        list(outcome=alpha,weight=dgp$gamma_params[[paste0(site,"_",value)]])),sources))
  })
  fitted <- list(models=models,sites=c("t",sources),configuration=list(radius=12))
  packet <- list(training_shifts=shifts*0,evaluation_sizes=setNames(rep(500,5L),fitted$sites),
    covariance=matrix(0,10L,10L,dimnames=list(paste0(rep(arms,each=5L),":",rep(c("target_anchor",sources),2L)),
      paste0(rep(arms,each=5L),":",rep(c("target_anchor",sources),2L)))))
  population <- honest_population_moments(fitted,packet,dgp,24L,shifts)
  quadrature <- RoCE:::.bounded_quadrature(24L)
  design <- cbind(1,quadrature$X)
  target_means <- sapply(list(dgp$alpha1_true,dgp$alpha0_true),function(alpha) plogis(drop(design%*%alpha)))
  expected <- rbind(colSums(target_means*quadrature$weights),
    t(vapply(sources,function(site) colSums(plogis(sweep(qlogis(target_means),2L,shifts[site,],"+"))*quadrature$weights),numeric(2L))))
  mean_error <- max(abs(population$means-expected))
  cross_error <- max(vapply(sources,function(site) {
    delta <- expected[match(site,sources)+1L,]-population$truth
    abs(population$source_moments[[site]]$covariance[1L,2L]+prod(delta)/500)
  },numeric(1L)))
  coarse <- honest_population_moments(fitted,packet,dgp,16L,shifts)
  quadrature_error <- max(abs(population$means-coarse$means),abs(population$covariance-coarse$covariance))
  stopifnot(mean_error<1e-12,cross_error<1e-12,quadrature_error<1e-8)
  rejected <- tryCatch({honest_population_moments(fitted,packet,dgp,8L);FALSE},error=function(error)
    grepl("Source shift does not match",conditionMessage(error),fixed=TRUE))
  stopifnot(rejected)
  errors[[mechanism]] <- data.frame(mechanism=mechanism,mean_error=mean_error,
    private_cross_covariance_error=cross_error,quadrature_error=quadrature_error,
    incorrect_source_law_rejected=rejected)
}
write.csv(do.call(rbind,errors),file.path(output,"oracle_checks.csv"),row.names=FALSE)
writeLines("SHIFTED_HONEST_POPULATION_CHECKS_PASSED",file.path(output,"checks.txt"))
