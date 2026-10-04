#!/usr/bin/env Rscript
# Prototype a sufficient nonconvergence certificate; production CV is unchanged.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: checked_R_library new_scratch_output")
library_path <- normalizePath(arguments[[1L]],mustWork=TRUE)
output <- arguments[[2L]]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
library(testthat)
dir.create(output,recursive=TRUE)

density_cv_kkt_bound <- function(design,moment,coefficients,kkt_tolerance) {
  design <- as.matrix(design)
  if(nrow(design)<1L || ncol(design)!=length(moment) || length(coefficients)!=length(moment) ||
     any(!is.finite(c(design,moment))) || any(design[,1L]!=1) ||
     length(kkt_tolerance)!=1L || !is.finite(kkt_tolerance) || kkt_tolerance<=0) stop("Invalid certificate inputs")
  empty <- list(penalty_bound=0,direction=rep(0,ncol(design)))
  if(ncol(design)<2L || moment[1L]<=0 || any(!is.finite(coefficients))) return(empty)
  scale <- max(abs(coefficients[-1L]))
  if(scale==0) return(empty)
  slopes <- coefficients[-1L]/scale
  slopes <- slopes/sum(abs(slopes))
  padding <- 64*.Machine$double.eps*(nrow(design)+ncol(design)+1L)*
    max(1,abs(design),abs(moment))
  projection <- drop(design[,-1L,drop=FALSE]%*%slopes)
  direction <- c(-min(projection)+padding,slopes)
  roundoff <- padding*(1+sum(abs(moment)))
  bound <- (-sum(moment*direction)-kkt_tolerance*sum(abs(direction))-roundoff)/sum(abs(slopes))
  if(!is.finite(bound) || any(!is.finite(direction))) return(empty)
  list(penalty_bound=max(0,bound),direction=direction)
}

kkt_residual <- function(coefficients,gradient,lambda) {
  slopes <- ifelse(coefficients[-1L]==0,pmax(abs(gradient[-1L])-lambda,0),
    abs(gradient[-1L]+lambda*sign(coefficients[-1L])))
  max(abs(gradient[1L]),slopes)
}

test_that("a multivariate separating direction certifies the KKT obstruction", {
  x <- seq(0,1,length.out=11)
  design <- cbind(1,x,1-x)
  bound <- density_cv_kkt_bound(design,c(1,0,0),c(0,1,1),1e-4)
  expect_equal(bound$penalty_bound,.5-1.5e-4,tolerance=1e-10)
  expect_gte(min(design%*%bound$direction),0)
  scaled <- density_cv_kkt_bound(design,c(1,0,0),c(99,1e200,1e200),1e-4)
  expect_equal(scaled,bound,tolerance=1e-14)
  # Although the mathematical objective is unbounded for lambda<.5, a point
  # can meet a loose CV KKT tolerance. The certificate must not reject it.
  gamma <- c(-1,1,1)
  gradient <- c(1,0,0)-colMeans(design*exp(-drop(design%*%gamma)))
  lambda <- .5-.5e-4
  expect_lte(kkt_residual(gamma,gradient,lambda),1e-4)
  expect_lte(bound$penalty_bound,lambda)
  set.seed(811)
  for(iteration in 1:20) {
    gamma <- rnorm(3)
    weights <- runif(nrow(design))*exp(-pmax(-12,pmin(12,drop(design%*%gamma))))
    gradient <- c(1,0,0)-colMeans(design*weights)
    expect_gt(kkt_residual(gamma,gradient,.1),1e-4)
  }
})

test_that("interior support and unusable directions do not generate certificates", {
  design <- cbind(1,as.matrix(expand.grid(x=c(-1,1),y=c(-1,1))))
  for(direction in list(c(0,1,1),c(0,-1,2),c(0,0,0),c(0,Inf,0))) {
    expect_equal(density_cv_kkt_bound(design,c(1,0,0),direction,1e-4)$penalty_bound,0)
  }
  expect_error(density_cv_kkt_bound(design,c(1,0,0),rep(0,3),0),"Invalid")
  extreme <- cbind(1,rep(1e308,4),rep(1e308,4))
  expect_equal(density_cv_kkt_bound(extreme,c(1e308,0,0),c(0,1,1),1e-4)$penalty_bound,0)
})

set.seed(2841)
features <- RoCE:::.draw_bounded_covariates(560L,200L)
target <- RoCE:::.draw_bounded_covariates(700L,200L)
arm <- rep(c(1L,0L),each=280L)
moment <- c(1,colMeans(target))
design <- cbind(1,features[arm==1L,,drop=FALSE])
previous_solver <- RoCE:::.set_nuisance_solver("proximal_newton")
elapsed <- system.time(initial <- RoCE:::fit_initial_density_ratio_cpp(features,arm,moment,
  lambda=.001,max_iter=100L,tol=1e-4,A_val=1L,M_tau=12,warm_start=numeric()))[["elapsed"]]
certificate <- density_cv_kkt_bound(design,moment,as.numeric(initial$gamma),1e-4)
record <- data.frame(n_source=560,n_arm=280,p=200,initial_lambda=.001,
  witness_fit_seconds=elapsed,witness_converged=initial$converged,
  penalty_bound=certificate$penalty_bound,probe_lambda=NA_real_,probe_converged=NA,
  probe_kkt=NA_real_,probe_seconds=NA_real_)
if(certificate$penalty_bound>0) {
  lambda <- certificate$penalty_bound/2
  elapsed <- system.time(probe <- RoCE:::fit_initial_density_ratio_cpp(features,arm,moment,
    lambda=lambda,max_iter=300L,tol=1e-4,A_val=1L,M_tau=12,warm_start=numeric()))[["elapsed"]]
  gamma <- as.numeric(probe$gamma)
  weights <- exp(-pmax(-12,pmin(12,drop(design%*%gamma))))
  independent_kkt <- kkt_residual(gamma,moment-colSums(design*weights)/nrow(features),lambda)
  stopifnot(!isTRUE(probe$converged),independent_kkt>1e-4,
    abs(independent_kkt-probe$kkt_residual)<1e-8)
  record$probe_lambda <- lambda;record$probe_converged <- probe$converged
  record$probe_kkt <- independent_kkt;record$probe_seconds <- elapsed
}
RoCE:::.restore_nuisance_solver(previous_solver)
write.csv(record,file.path(output,"certificate_probe.csv"),row.names=FALSE)
saveRDS(list(certificate=certificate,witness=initial),file.path(output,"witness.rds"))
writeLines(c("DENSITY_CV_CERTIFICATE_PROTOTYPE_CHECKS_PASSED",
  "Sufficient KKT obstruction only. Full descending-CV-path parity and production integration remain untested."),
  file.path(output,"checks.txt"))
print(record,row.names=FALSE)
