#!/usr/bin/env Rscript
# Population check: weight correctness can protect orthogonality even when
# initial and calibrated OR projections have different probability limits.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: checked_R_library new_scratch_output")
library_path <- normalizePath(arguments[[1L]],mustWork=TRUE)
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
output <- arguments[[2L]]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
dir.create(output,recursive=TRUE)
evaluate <- function(order) {
  quadrature <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1,quadrature$X)
  mass <- quadrature$weights
  propensity_slopes <- c(.35,-.25,.125,-.0625)
  target_propensity <- plogis(drop(quadrature$X %*% propensity_slopes))
  rows <- list()
  for(scenario in c("C1","C2")) {
    strength <- RoCE:::.face_misspecification_strengths(scenario,.75)$outcome
    predictor <- drop(RoCE:::.bounded_active_features(quadrature$X,strength) %*% c(.25,.20,.10,.05))
    for(arm in 0:1) {
      true_mean <- plogis(-.5+arm+predictor)
      initial_mass <- mass*if(arm==1L) target_propensity else 1-target_propensity
      initial <- glm.fit(design,true_mean,weights=initial_mass,family=quasibinomial(),
                         control=glm.control(epsilon=1e-12,maxit=100))
      calibrated <- glm.fit(design,true_mean,weights=mass,family=quasibinomial(),
                            control=glm.control(epsilon=1e-12,maxit=100))
      stopifnot(initial$converged,calibrated$converged)
      initial_derivative <- initial$fitted.values*(1-initial$fitted.values)
      final_derivative <- calibrated$fitted.values*(1-calibrated$fitted.values)
      gap <- max(abs(initial$coefficients-calibrated$coefficients))
      if(scenario=="C1") stopifnot(gap<1e-9) else stopifnot(gap>1e-4)
      for(site in 1:2) {
        common <- site/2*c(.25,-.20,.15,-.10)
        slopes <- cbind(common-propensity_slopes/2,common+propensity_slopes/2)
        exponential <- exp(quadrature$X %*% slopes)
        normalizer <- sum(mass*rowSums(exponential))
        joint <- exponential[,arm+1L]/normalizer
        log_weight <- c(log(normalizer),-slopes[,arm+1L])
        weights <- exp(drop(design %*% log_weight))
        predictor_bound <- abs(log_weight[1L])+RoCE:::.bounded_covariate_spec()$bound*sum(abs(log_weight[-1L]))
        stopifnot(predictor_bound<12)
        initial_weight_gradient <- drop(crossprod(design,mass*(joint*weights-1)))
        calibrated_weight_gradient <- drop(crossprod(design,mass*initial_derivative*(joint*weights-1)))
        outcome_gradient <- drop(crossprod(design,mass*final_derivative*(1-joint*weights)))
        weight_gradient <- drop(crossprod(design,mass*joint*weights*(true_mean-calibrated$fitted.values)))
        score <- function(alpha,log_weight_coefficients) {
          predicted <- plogis(drop(design %*% alpha))
          fitted_weight <- exp(drop(design %*% log_weight_coefficients))
          sum(mass*(predicted+joint*fitted_weight*(true_mean-predicted)))
        }
        numerical <- vapply(seq_len(ncol(design)),function(index) {
          step <- numeric(ncol(design)); step[index] <- 1e-5
          (score(calibrated$coefficients,log_weight+step)-score(calibrated$coefficients,log_weight-step))/2e-5
        },numeric(1L))
        stopifnot(max(abs(c(initial_weight_gradient,calibrated_weight_gradient,outcome_gradient,weight_gradient)))<1e-9,
                  max(abs(numerical-weight_gradient))<1e-8)
        rows[[length(rows)+1L]] <- data.frame(config=scenario,arm=arm,site=site,initial_final_or_gap=gap,
          outcome_score_gradient=max(abs(outcome_gradient)),weight_score_gradient=max(abs(weight_gradient)),
          initial_weight_loss_gradient=max(abs(initial_weight_gradient)),calibrated_weight_loss_gradient=max(abs(calibrated_weight_gradient)),
          point_bias=score(calibrated$coefficients,log_weight)-sum(mass*true_mean),
          finite_difference_error=max(abs(numerical-weight_gradient)),weight_predictor_bound=predictor_bound)
      }
    }
  }
  do.call(rbind,rows)
}
coarse <- evaluate(16L); fine <- evaluate(24L)
fields <- c("initial_final_or_gap","outcome_score_gradient","weight_score_gradient","point_bias")
error <- max(abs(as.matrix(coarse[fields])-as.matrix(fine[fields])))
stopifnot(error<1e-8)
write.csv(fine,file.path(output,"population_checks.csv"),row.names=FALSE)
writeLines(c(paste("Maximum 16/24 quadrature discrepancy:",error),
  "Population score check only; no fitted-nuisance convergence rate is asserted."),file.path(output,"checks.txt"))
cat("UNION_ORTHOGONALITY_POPULATION_CHECKS_PASSED\n")
