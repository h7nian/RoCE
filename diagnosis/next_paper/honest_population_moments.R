# Conditional population moments for p4 bounded research candidates.
# Quadrature is a diagnostic oracle, never a feasible-estimator input.
honest_population_moments <- function(fitted,packet,dgp,order,source_logit_shifts=NULL) {
  stopifnot(dgp$p==4L,dgp$dgp_version=="bounded_joint_v3")
  quadrature <- RoCE:::.bounded_quadrature(order)
  design <- cbind(1,quadrature$X)
  mass <- quadrature$weights
  strengths <- RoCE:::.face_misspecification_strengths(dgp$config,dgp$misspecification_strength)
  true_outcome_design <- cbind(1,RoCE:::.bounded_active_features(quadrature$X,strengths$outcome))
  true_assignment <- RoCE:::.bounded_active_features(quadrature$X,strengths$propensity)
  propensity_slopes <- c(.35,-.25,.125,-.0625)
  target_probability <- plogis(drop(true_assignment %*% propensity_slopes))
  true_mean <- cbind(mu1=plogis(drop(true_outcome_design %*% dgp$alpha1_true)),
                     mu0=plogis(drop(true_outcome_design %*% dgp$alpha0_true)))
  truth <- colSums(true_mean*mass)
  arms <- c("mu1","mu0"); sources <- setdiff(fitted$sites,"t")
  if(is.null(source_logit_shifts)) {
    source_logit_shifts <- matrix(0,length(sources),2L,dimnames=list(sources,arms))
  }
  if(!is.matrix(source_logit_shifts) || !is.numeric(source_logit_shifts) ||
     !identical(dim(source_logit_shifts),c(length(sources),2L)) ||
     !identical(rownames(source_logit_shifts),sources) ||
     !identical(colnames(source_logit_shifts),arms) || any(!is.finite(source_logit_shifts))) {
    stop("source_logit_shifts must have matching source rows and mu1/mu0 columns")
  }
  for(site in sources) {
    rows <- which(dgp$R==site)
    observed_design <- cbind(1,dgp$W_outcome_true[rows,,drop=FALSE])
    for(arm in arms) {
      coefficients <- if(arm=="mu1") dgp$alpha1_true else dgp$alpha0_true
      expected <- plogis(drop(observed_design%*%coefficients)+source_logit_shifts[site,arm])
      if(max(abs(expected-dgp$oracle[[arm]][rows]))>1e-12) {
        stop("Source shift does not match the generated outcome law: ",site,"/",arm)
      }
    }
  }
  predict <- function(coefficients) plogis(drop(design %*% coefficients))
  clipped <- function(coefficients) pmax(-fitted$configuration$radius,
    pmin(fitted$configuration$radius,drop(design %*% coefficients)))
  anchor_predictions <- sapply(arms,function(arm) predict(fitted$models[[arm]]$target$outcome))
  common_predictions <- sapply(arms,function(arm) predict(fitted$models[[arm]]$common$coefficients))
  propensities <- sapply(arms,function(arm) plogis(clipped(fitted$models[[arm]]$target$weight)))
  base <- cbind(anchor_predictions,common_predictions)
  target_mean <- numeric(4L); target_second <- matrix(0,4L,4L)
  # Enumerate actual (A,Y) laws; the unobserved potential-outcome dependence
  # is irrelevant because only one arm's residual enters each observation.
  for(index in 1:2) for(outcome in 0:1) {
    probability <- if(index==1L) target_probability else 1-target_probability
    weights <- mass*probability*if(outcome==1L) true_mean[,index] else 1-true_mean[,index]
    scores <- base
    scores[,index] <- scores[,index]+(outcome-anchor_predictions[,index])/propensities[,index]
    target_mean <- target_mean+colSums(scores*weights)
    target_second <- target_second+crossprod(scores,scores*weights)
  }
  target_covariance <- (target_second-tcrossprod(target_mean))/packet$evaluation_sizes[["t"]]
  width <- length(fitted$sites)
  columns <- c(1L,rep(3L,length(sources)),2L,rep(4L,length(sources)))
  covariance <- target_covariance[columns,columns]
  means <- original_means <- matrix(0,width,2L,dimnames=list(c("target_anchor",sources),arms))
  means[1L,] <- original_means[1L,] <- target_mean[1:2]
  source_moments <- list()
  identity_error <- 0
  for(index in seq_along(sources)) {
    site <- sources[index]
    parameters <- dgp$oracle$source_parameters[[site]]
    joint <- exp(true_assignment %*% parameters$slopes-parameters$log_normalizer)
    stopifnot(abs(sum(mass*rowSums(joint))-1)<1e-8)
    residual_mean <- residual_second <- original_target <- numeric(2L)
    for(arm in seq_along(arms)) {
      model <- fitted$models[[arms[arm]]]$source[[site]]
      regression <- predict(model$outcome)
      weight <- exp(-clipped(model$weight))
      joint_mass <- mass*joint[,if(arm==1L) 2L else 1L]
      source_mean <- true_mean[,arm]
      if(source_logit_shifts[site,arm]!=0) {
        source_mean <- plogis(qlogis(source_mean)+source_logit_shifts[site,arm])
      }
      difference <- source_mean-regression
      residual_mean[arm] <- sum(joint_mass*weight*difference)
      residual_second[arm] <- sum(joint_mass*weight^2*(source_mean*(1-source_mean)+difference^2))
      original_target[arm] <- sum(mass*regression)
    }
    private <- (diag(residual_second)-tcrossprod(residual_mean))/packet$evaluation_sizes[[site]]
    positions <- c(index+1L,width+index+1L)
    covariance[positions,positions] <- covariance[positions,positions]+private
    original_means[site,] <- original_target+residual_mean
    means[site,] <- target_mean[3:4]+residual_mean+packet$training_shifts[site,]
    shift <- packet$training_shifts[site,]-(original_target-target_mean[3:4])
    identity_error <- max(identity_error,abs(means[site,]-original_means[site,]-shift))
    source_moments[[site]] <- list(mean=residual_mean,covariance=private)
  }
  dimnames(covariance) <- dimnames(packet$covariance)
  stopifnot(max(abs(truth-c(dgp$mu1_true,dgp$mu0_true)))<1e-8,
    min(eigen(covariance,symmetric=TRUE,only.values=TRUE)$values)>-1e-12,identity_error<1e-12)
  list(means=means,original_means=original_means,truth=truth,covariance=covariance,
       target_covariance=target_covariance,source_moments=source_moments,identity_error=identity_error)
}
