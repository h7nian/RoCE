#!/usr/bin/env Rscript

# Deterministic population example, not Monte Carlo or a change of estimator.
# A correct quadratic outcome model with a reduced linear tilting basis can
# satisfy the implemented calibration equations without zeroing the full
# derivative of the source-assisted mean with respect to outcome coefficients.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: check_calibration_tangent_span.R V19_LIBRARY")
  lib <- normalizePath(args[[1L]], mustWork=TRUE)
  .libPaths(c(lib,.libPaths()))
  suppressPackageStartupMessages(library(RoCE,lib.loc=lib))
  stopifnot(normalizePath(find.package("RoCE"))==file.path(lib,"RoCE"))
  source("scripts/slurm/result_provenance.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  hash <- .weight_calibration_installed_package_fingerprint(lib,roce_sha256_file)
  stopifnot(hash=="2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")

  x <- c(-1,0,1)
  W <- cbind(intercept=1,linear=x,quadratic=x^2)
  Z <- W[,1:2,drop=FALSE]
  target_mass <- rep(1/3,3)
  source_treated_mass <- c(.08,.28,.12)
  alpha <- c(-.2,.7,1.2)
  true_mean <- plogis(drop(W%*%alpha))
  derivative <- true_mean*(1-true_mean)
  target_moment <- drop(crossprod(Z,target_mass*derivative))
  objective <- function(gamma) sum(target_mass*derivative*drop(Z%*%gamma))+
    sum(source_treated_mass*derivative*exp(-drop(Z%*%gamma)))
  gradient <- function(gamma) drop(crossprod(Z,derivative*(target_mass-
    source_treated_mass*exp(-drop(Z%*%gamma)))))
  fitted <- optim(c(0,0),objective,gradient,method="BFGS",
                  control=list(reltol=1e-14,maxit=10000))
  stopifnot(fitted$convergence==0L,max(abs(gradient(fitted$par)))<1e-8)

  # Exact masses represented as source rows for the installed native kernel.
  source_x <- rep(x,c(90L,120L,90L))
  source_a <- c(rep(c(1,0),c(24L,66L)),rep(c(1,0),c(84L,36L)),
                rep(c(1,0),c(36L,54L)))
  native <- RoCE:::fit_unified_density_ratio_cpp(
    matrix(source_x,ncol=1),source_a,target_moment,alpha,0,10000L,1e-12,
    TRUE,5,cbind(source_x,source_x^2),1L,1L,1L,numeric())
  native_gamma <- as.numeric(native$gamma)
  gamma_error <- max(abs(native_gamma-fitted$par))
  balance <- gradient(native_gamma)
  weights <- exp(-drop(Z%*%native_gamma))
  outcome_gradient <- drop(crossprod(W,derivative*(target_mass-source_treated_mass*weights)))
  expected_score <- function(a) {
    prediction <- plogis(drop(W%*%a))
    sum(target_mass*prediction)+sum(source_treated_mass*weights*(true_mean-prediction))
  }
  epsilon <- 1e-5
  numeric_gradient <- vapply(seq_along(alpha),function(j) {
    direction <- rep(0,3); direction[j]<-epsilon
    (expected_score(alpha+direction)-expected_score(alpha-direction))/(2*epsilon)
  },numeric(1L))
  stopifnot(gamma_error<1e-6,max(abs(balance))<1e-8,
    abs(expected_score(alpha)-sum(target_mass*true_mean))<1e-12,
    max(abs(outcome_gradient-numeric_gradient))<1e-8,
    abs(outcome_gradient[3])>1e-3)
  result <- data.frame(direction=colnames(W),
    covered_by_linear_calibration=c(TRUE,TRUE,FALSE),
    analytic_score_derivative=outcome_gradient,
    finite_difference_derivative=numeric_gradient)
  cat("installed_v19_native_vs_independent_optimizer_max_error=",gamma_error,"\n",sep="")
  cat("implemented_calibration_moment_max_error=",max(abs(balance)),"\n",sep="")
  print(result,row.names=FALSE,digits=12)
  cat("Correct outcome gives zero population estimation bias at its true coefficients,\n",
      "but the quadratic nuisance tangent remains first-order sensitive.\n",
      "This is a structural counterexample, not a quantified explanation of the C1 pilot\n",
      "and not a proposed union-basis repair or a coverage validation.\n",sep="")
  invisible(result)
}

if(sys.nframe()==0L) main()
