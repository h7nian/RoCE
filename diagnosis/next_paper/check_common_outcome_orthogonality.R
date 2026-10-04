#!/usr/bin/env Rscript
# Independent population equations for the proposed common-outcome score.
# This unpenalized, five-parameter solver is used only for the population audit.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: checked_R_library new_scratch_output")
library_path <- normalizePath(arguments[[1L]],mustWork=TRUE)
output <- arguments[[2L]]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
dir.create(output,recursive=TRUE)

solve_population_calibration <- function(design,source_mass,target_moment) {
  gamma <- c(log(sum(source_mass)/target_moment[1L]),rep(0,ncol(design)-1L))
  objective <- function(value) sum(source_mass*exp(-drop(design%*%value)))+sum(target_moment*value)
  for(iteration in seq_len(100L)) {
    tilted <- source_mass*exp(-drop(design%*%gamma))
    gradient <- target_moment-drop(crossprod(design,tilted))
    if(max(abs(gradient))<1e-11) return(gamma)
    direction <- solve(crossprod(design,design*tilted),gradient)
    step <- 1
    current <- objective(gamma)
    roundoff <- 16*.Machine$double.eps*max(1,abs(current))
    while(objective(gamma-step*direction)>current-1e-4*step*sum(gradient*direction)+roundoff) {
      step <- step/2
      if(step<2^-30) stop("Population Newton line search failed")
    }
    gamma <- gamma-step*direction
  }
  stop("Population calibration did not converge")
}

audit <- function(order) {
  quadrature <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1,quadrature$X)
  mass <- quadrature$weights
  propensity_slopes <- c(.35,-.25,.125,-.0625)
  outcome_slopes <- c(.25,.20,.10,.05)
  rows <- list()
  for(scenario in c("C1","C2","C3")) {
    strength <- RoCE:::.face_misspecification_strengths(scenario,.75)
    assignment <- RoCE:::.bounded_active_features(quadrature$X,strength$propensity)
    target_propensity <- plogis(drop(assignment%*%propensity_slopes))
    outcome_predictor <- drop(RoCE:::.bounded_active_features(quadrature$X,strength$outcome)%*%outcome_slopes)
    for(arm in 0:1) {
      probability <- if(arm==1L) target_propensity else 1-target_propensity
      true_mean <- plogis(-.5+arm+outcome_predictor)
      initial_ps <- solve_population_calibration(design,mass*probability,
        drop(crossprod(design,mass*(1-probability))))
      working_probability <- plogis(drop(design%*%initial_ps))
      fit_projection <- function(weights) {
        fit <- glm.fit(design,true_mean,weights=weights,family=quasibinomial(),
          control=glm.control(epsilon=1e-12,maxit=100))
        stopifnot(fit$converged)
        fit$coefficients
      }
      common_beta <- fit_projection(mass*probability/working_probability)
      initial_beta <- fit_projection(mass*probability)
      initial_mean <- plogis(drop(design%*%initial_beta))
      derivative_weight <- initial_mean*(1-initial_mean)
      if(scenario!="C2") stopifnot(max(abs(common_beta-c(-.5+arm,outcome_slopes)))<1e-9)
      if(scenario!="C3") stopifnot(max(abs(working_probability-probability))<1e-9)
      for(site in 1:2) {
        tilt <- site/2*c(.25,-.20,.15,-.10)
        slopes <- cbind(tilt-propensity_slopes/2,tilt+propensity_slopes/2)
        exponential <- exp(assignment%*%slopes)
        joint <- exponential[,arm+1L]/sum(mass*rowSums(exponential))
        calibrated_gamma <- solve_population_calibration(design,mass*joint*derivative_weight,
          drop(crossprod(design,mass*derivative_weight)))
        variants <- list(common_outcome=list(beta=common_beta,gamma=calibrated_gamma))
        if(scenario=="C2") variants$wrong_target_projection <- list(beta=initial_beta,gamma=calibrated_gamma)
        if(scenario=="C3") variants$uncalibrated_weight <- list(beta=common_beta,
          gamma=solve_population_calibration(design,mass*joint,drop(crossprod(design,mass))))
        for(variant in names(variants)) {
          beta <- variants[[variant]]$beta; gamma <- variants[[variant]]$gamma
          regression <- plogis(drop(design%*%beta))
          weight <- exp(-drop(design%*%gamma))
          gradient_beta <- drop(crossprod(design,mass*(1-joint*weight)*regression*(1-regression)))
          gradient_gamma <- -drop(crossprod(design,mass*joint*weight*(true_mean-regression)))
          score <- function(b,g) {
            predicted <- plogis(drop(design%*%b))
            sum(mass*(predicted+joint*exp(-drop(design%*%g))*(true_mean-predicted)))
          }
          numerical <- sapply(seq_len(ncol(design)),function(column) {
            step <- numeric(ncol(design));step[column] <- 1e-5
            c(beta=(score(beta+step,gamma)-score(beta-step,gamma))/2e-5,
              gamma=(score(beta,gamma+step)-score(beta,gamma-step))/2e-5)
          })
          error <- max(abs(numerical-rbind(gradient_beta,gradient_gamma)))
          bound <- abs(gamma[1L])+RoCE:::.bounded_covariate_spec()$bound*sum(abs(gamma[-1L]))
          stopifnot(error<1e-8,bound<12,abs(score(beta,gamma)-sum(mass*true_mean))<1e-9)
          if(variant=="common_outcome") stopifnot(max(abs(c(gradient_beta,gradient_gamma)))<1e-8)
          else stopifnot(max(abs(c(gradient_beta,gradient_gamma)))>1e-6)
          rows[[length(rows)+1L]] <- data.frame(config=scenario,arm=arm,site=site,variant=variant,
            beta_score_gradient=max(abs(gradient_beta)),gamma_score_gradient=max(abs(gradient_gamma)),
            finite_difference_error=error,weight_predictor_bound=bound,
            initial_common_beta_gap=max(abs(initial_beta-common_beta)))
        }
      }
    }
  }
  do.call(rbind,rows)
}
coarse <- audit(16L); fine <- audit(24L)
fields <- c("beta_score_gradient","gamma_score_gradient","initial_common_beta_gap")
discrepancy <- max(abs(as.matrix(coarse[fields])-as.matrix(fine[fields])))
stopifnot(discrepancy<1e-8)
write.csv(fine,file.path(output,"population_gradients.csv"),row.names=FALSE)
writeLines(c("COMMON_OUTCOME_POPULATION_ORTHOGONALITY_CHECKS_PASSED",
  paste("16/24 quadrature discrepancy:",discrepancy),
  "Positive and negative population controls; not fitted-nuisance convergence or CI coverage validation."),
  file.path(output,"checks.txt"))
cat("COMMON_OUTCOME_POPULATION_ORTHOGONALITY_CHECKS_PASSED\n")
