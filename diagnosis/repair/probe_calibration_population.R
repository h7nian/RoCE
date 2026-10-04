#!/usr/bin/env Rscript
# Population mechanism check. This does not select a DGP by Monte Carlo coverage.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: checked_R_library new_scratch_output")
library_path <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
settings <- rbind(data.frame(config="C1",strength=0,shift=c(1,2)),
  expand.grid(config=c("C2","C3"),strength=c(.75,1),shift=c(1,2),stringsAsFactors=FALSE))
laws <- lapply(seq_len(nrow(settings)),function(index) {
  setting <- settings[index,]
  set.seed(901L)
  generate_bounded_data(n_target=1000L,n_source_sizes=c(1000L,1000L),p=4L,
    config=setting$config,misspecification_strength=setting$strength,
    dgp_control=list(covariate_shift=setting$shift,source_treatment_scale=1))
})
evaluate <- function(order) {
  grid <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1,grid$X); mass <- grid$weights
  bound <- RoCE:::.bounded_covariate_spec()$bound
  regression <- function(response,weights,family) {
    fit <- glm.fit(design,response,weights=weights,family=family,
      control=glm.control(epsilon=1e-12,maxit=100))
    stopifnot(fit$converged,all(is.finite(fit$coefficients)))
    fit$coefficients
  }
  rows <- list()
  for(index in seq_along(laws)) {
    law <- laws[[index]]; setting <- settings[index,]
    strengths <- RoCE:::.face_misspecification_strengths(law$config,law$misspecification_strength)
    outcome_design <- cbind(1,RoCE:::.bounded_active_features(grid$X,strengths$outcome))
    assignment <- RoCE:::.bounded_active_features(grid$X,strengths$propensity)
    parameters <- law$oracle$source_parameters$s2
    theta <- parameters$slopes[,"mu1"]-parameters$slopes[,"mu0"]
    target_propensity <- plogis(drop(assignment%*%theta))
    joint <- exp(assignment%*%parameters$slopes-parameters$log_normalizer)
    stopifnot(abs(sum(mass*rowSums(joint))-1)<1e-8)
    propensity_range <- plogis(RoCE:::.bounded_predictor_range(theta,strengths$propensity))
    for(arm in 0:1) {
      true_mean <- plogis(drop(outcome_design%*%if(arm==1L) law$alpha1_true else law$alpha0_true))
      truth <- sum(mass*true_mean)
      probability <- if(arm==1L) target_propensity else 1-target_propensity
      joint_arm <- joint[,arm+1L]
      initial_outcome <- regression(true_mean,mass*probability,quasibinomial())
      prediction <- plogis(drop(design%*%initial_outcome))
      derivative <- prediction*(1-prediction)
      # Poisson score equations with artificial response1/joint reproduce
      # the unpenalized population tilt loss. gamma is minus the log-rate.
      initial_weight <- -regression(1/joint_arm,mass*joint_arm,quasipoisson())
      calibrated_weight <- -regression(1/joint_arm,mass*joint_arm*derivative,quasipoisson())
      ordinary_outcome <- regression(true_mean,mass*joint_arm,quasibinomial())
      calibrated_outcome <- regression(true_mean,mass*joint_arm*exp(-drop(design%*%initial_weight)),quasibinomial())
      initial_gradient <- crossprod(design,mass*(1-joint_arm*exp(-drop(design%*%initial_weight))))
      calibrated_gradient <- crossprod(design,mass*derivative*(1-joint_arm*exp(-drop(design%*%calibrated_weight))))
      stopifnot(max(abs(c(initial_gradient,calibrated_gradient)))<1e-9)
      gram_min <- min(eigen(crossprod(design,design*(mass*joint_arm)),symmetric=TRUE,only.values=TRUE)$values)
      for(program in c("standard","calibrated")) {
        alpha <- if(program=="standard") ordinary_outcome else calibrated_outcome
        gamma <- if(program=="standard") initial_weight else calibrated_weight
        predicted <- plogis(drop(design%*%alpha)); q <- exp(-drop(design%*%gamma))
        outcome_gradient <- drop(crossprod(design,mass*(1-joint_arm*q)*predicted*(1-predicted)))
        weight_gradient <- -drop(crossprod(design,mass*joint_arm*q*(true_mean-predicted)))
        score <- function(alpha,gamma) {
          prediction <- plogis(drop(design%*%alpha))
          sum(mass*(prediction+joint_arm*exp(-drop(design%*%gamma))*(true_mean-prediction)))
        }
        numerical <- lapply(c("outcome","weight"),function(type) vapply(seq_len(ncol(design)),function(coordinate) {
          delta <- numeric(ncol(design)); delta[coordinate] <- 1e-5
          if(type=="outcome") (score(alpha+delta,gamma)-score(alpha-delta,gamma))/2e-5 else
            (score(alpha,gamma+delta)-score(alpha,gamma-delta))/2e-5
        },numeric(1L)))
        gradient_error <- max(abs(numerical[[1L]]-outcome_gradient),abs(numerical[[2L]]-weight_gradient))
        predictor_bound <- max(vapply(list(initial_outcome,alpha,gamma),function(beta)
          abs(beta[1L])+bound*sum(abs(beta[-1L])),numeric(1L)))
        score_mean <- score(alpha,gamma)
        if(program=="calibrated" || setting$config=="C1") {
          stopifnot(max(abs(c(outcome_gradient,weight_gradient)))<1e-8)
        }
        stopifnot(gradient_error<1e-8,abs(score_mean-truth)<1e-8,
          predictor_bound<12,gram_min>0,propensity_range[1L]>0,propensity_range[2L]<1)
        rows[[length(rows)+1L]] <- data.frame(config=setting$config,
          misspecification_strength=setting$strength,covariate_shift=setting$shift,
          source="s2",source_program=program,arm=arm,population_mean=score_mean,truth=truth,
          population_bias=score_mean-truth,outcome_gradient_max=max(abs(outcome_gradient)),
          weight_gradient_max=max(abs(weight_gradient)),finite_difference_error=gradient_error,
          fitted_predictor_bound=predictor_bound,source_arm_probability=sum(mass*joint_arm),
          propensity_min=propensity_range[1L],propensity_max=propensity_range[2L],gram_min_eigenvalue=gram_min)
      }
    }
    cat("Population check",order,setting$config,setting$strength,setting$shift,"\n")
  }
  do.call(rbind,rows)
}
coarse <- evaluate(16L); fine <- evaluate(24L)
fields <- c("population_mean","truth","outcome_gradient_max","weight_gradient_max","fitted_predictor_bound")
error <- max(abs(as.matrix(coarse[fields])-as.matrix(fine[fields])))
stopifnot(error<1e-8)
dir.create(output,recursive=TRUE)
write.csv(fine,file.path(output,"population_gradients.csv"),row.names=FALSE)
write.csv(settings,file.path(output,"prespecified_settings.csv"),row.names=FALSE)
saveRDS(list(laws=lapply(laws,function(law) law[c("dgp_version","config","misspecification_strength","dgp_control")]),
  quadrature_error=error,library=library_path),file.path(output,"validation.rds"))
writeLines(c("CALIBRATION_POPULATION_GRADIENT_CHECKS_PASSED",
  "Oracle population projections on the four active coordinates; inactive independent coordinates have zero population gradients.",
  "Transportable sources only. This checks nuisance sensitivity, not finite-sample MSE or confidence-interval coverage.",
  "Current settings and stronger in-family controls were prespecified before fitting; no production DGP changed.",
  paste("Maximum16/24 quadrature difference:",error)),file.path(output,"checks.txt"))
