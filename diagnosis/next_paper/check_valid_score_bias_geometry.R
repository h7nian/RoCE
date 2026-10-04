#!/usr/bin/env Rscript
# Independent finite-matrix checks for the draft bias-budget construction.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: repository_root new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
source(file.path(repository,"diagnosis/next_paper/dispersion_bias_intervals.R"))
source(file.path(repository,"diagnosis/next_paper/joint_dispersion_intervals.R"))
dir.create(output,recursive=TRUE)
audit_checks <- 0L
check <- function(condition,label) {
  audit_checks <<- audit_checks+1L
  if(!isTRUE(condition)) stop(label,call.=FALSE)
}
close <- function(actual,expected,label,tolerance=1e-9) {
  check(all(is.finite(actual)) && max(abs(actual-expected))<=tolerance*max(1,abs(expected)),label)
}
quadratic <- function(x,matrix) drop(crossprod(x,matrix%*%x))

pool_coefficients <- function(common,private_pool,counts) {
  contrast <- c(1,-1)
  anchor <- common[1:2,1:2]
  cross <- common[1:2,3:4]
  prediction <- common[3:4,3:4]
  curvature <- anchor+prediction-cross-t(cross)
  right <- drop((anchor-t(cross))%*%contrast)
  active <- which(counts>0L)
  beta <- numeric(2L)
  if(length(active)) {
    matrix <- curvature[active,active,drop=FALSE]+diag(private_pool[active],nrow=length(active))
    beta[active] <- solve(matrix,right[active])
  }
  list(anchor=contrast-beta,source=beta)
}

audit_covariance <- function(covariance,counts,label) {
  width <- nrow(covariance)/2L
  count <- width-1L
  calibration <- make_joint_dispersion_calibration(covariance,counts[1L],counts[2L],"shared_prediction_bound")
  design <- calibration$design
  known <- calibration$known_valid
  selected_design <- design[known,,drop=FALSE]
  precision <- solve(covariance[known,known,drop=FALSE])
  information <- crossprod(selected_design,precision%*%selected_design)
  residual <- matrix(0,2L*width,2L*width)
  residual[known,known] <- precision-precision%*%selected_design%*%solve(information)%*%t(selected_design)%*%precision
  projection <- diag(2L*width)-residual%*%covariance
  close(residual%*%design,0,paste(label,"reference residual annihilates means"))
  close(residual%*%covariance%*%residual,residual,paste(label,"residual precision identity"))
  if(nrow(calibration$bias_map)>0L) {
    close(calibration$bias_map%*%design,0,paste(label,"bias map annihilates means"))
    close(calibration$bias_map%*%covariance%*%residual,0,paste(label,"bias-map/reference orthogonality"))
  }
  inflated <- inflate_shared_private_covariance(covariance)$covariance
  private <- shared_prediction_variances(inflated)$private_variances
  check(min(private)>0,paste(label,"positive private variances"))
  anchors <- c(1L,width+1L)
  source_positions <- list(2:width,width+(2:width))
  common <- matrix(0,4L,4L)
  common[1:2,1:2] <- inflated[anchors,anchors]
  common[1:2,3:4] <- inflated[anchors,c(2L,width+2L)]
  common[3:4,1:2] <- t(common[1:2,3:4])
  common[3L,3L] <- inflated[2L,3L]
  common[4L,4L] <- inflated[width+2L,width+3L]
  common[3L,4L] <- common[4L,3L] <- inflated[2L,width+3L]
  budgets <- seq(.003,.017,length.out=2L*width)
  budgets[anchors] <- c(.009,.012)
  transformed <- drop(crossprod(abs(projection),budgets))
  source_budget <- vapply(source_positions,function(index) max(transformed[index]),numeric(1L))
  envelope <- function(pool) {
    coefficients <- pool_coefficients(common,pool,counts)
    sum(abs(coefficients$anchor)*transformed[anchors])+sum(abs(coefficients$source)*source_budget)
  }
  limits <- lapply(1:2,function(arm) {
    if(counts[arm]==0L) return(0)
    ordered <- sort(private[,arm])
    unique(c(1/sum(1/head(ordered,counts[arm])),1/sum(1/tail(ordered,counts[arm]))))
  })
  corners <- expand.grid(limits)
  corner_bound <- max(apply(corners,1L,envelope))
  worst_sites <- lapply(1:2,function(arm) head(order(private[,arm],decreasing=TRUE),counts[arm]))
  worst_indices <- c(anchors,1L+worst_sites[[1L]],width+1L+worst_sites[[2L]])
  variance_bound <- joint_mean_gls(inflated,design,worst_indices)$variance
  reference_budget <- sum(abs(calibration$reference_coefficients)*budgets)
  reference_factor <- calibration$reference_variance-calibration$variance
  shared_factor <- variance_bound-calibration$variance
  close(min(reference_factor,shared_factor),calibration$bias_factor,paste(label,"zero-budget variance bound"))
  if("valid_score_bias"%in%names(formals(make_joint_dispersion_calibration))) {
    budgeted <- make_joint_dispersion_calibration(covariance,counts[1L],counts[2L],
      "shared_prediction_bound",valid_score_bias=budgets)
    close(budgeted$reference_bias_allowance,reference_budget,paste(label,"implemented reference allowance"))
    if(length(calibration$uncertain) && !calibration$anchor_only) {
      close(budgeted$valid_subset_bias_bounds[1L,"bias"],corner_bound,paste(label,"implemented corner envelope"))
      close(budgeted$valid_subset_bias_bounds[1L,"variance"],shared_factor,paste(label,"uncapped paired variance"))
    }
  }
  sets <- lapply(counts,function(number) if(number==0L) list(integer()) else combn(count,number,simplify=FALSE))
  maximum_envelope_ratio <- 0
  maximum_projection_error <- 0
  unprojected_failures <- 0L
  subset_count <- 0L
  for(first in sets[[1L]]) for(second in sets[[2L]]) {
    subset_count <- subset_count+1L
    indices <- sort(c(anchors,first+1L,second+width+1L))
    check(all(known%in%indices),paste(label,"subset contains reference"))
    comparison <- joint_mean_gls(inflated,design,indices)
    coefficients <- comparison$coefficients
    projected <- drop(projection%*%coefficients)
    difference <- projected-calibration$coefficients
    close(crossprod(design,projected),c(1,-1),paste(label,"projected design constraint"))
    outside <- setdiff(seq_len(2L*width),indices)
    if(length(outside)) close(projected[outside],0,paste(label,"projected support"))
    close(residual%*%covariance%*%projected,0,paste(label,"projected reference residuals"))
    if(nrow(calibration$bias_map)>0L) {
      reconstruction <- drop(t(calibration$bias_map)%*%calibration$bias_precision%*%
        calibration$bias_map%*%covariance%*%difference)
      maximum_projection_error <- max(maximum_projection_error,abs(reconstruction-difference))
      close(reconstruction,difference,paste(label,"reduced-map row space"))
    } else close(difference,0,paste(label,"all-valid projection is full GLS"))
    variance <- quadratic(projected,covariance)
    raw_variance <- quadratic(coefficients,covariance)
    check(variance<=raw_variance+1e-12 && raw_variance<=comparison$variance+1e-12,
      paste(label,"projection and inflation variance order"))
    check(comparison$variance<=variance_bound+1e-12,paste(label,"worst subset variance"))
    close(quadratic(difference,covariance),variance-calibration$variance,paste(label,"BLUE variance difference"))
    pool <- c(if(length(first)) 1/sum(1/private[first,1L]) else 0,
              if(length(second)) 1/sum(1/private[second,2L]) else 0)
    collapsed <- pool_coefficients(common,pool,counts)
    close(coefficients[anchors],collapsed$anchor,paste(label,"collapsed anchor coefficients"))
    close(c(sum(coefficients[first+1L]),sum(coefficients[second+width+1L])),collapsed$source,
      paste(label,"collapsed source coefficients"))
    actual_budget <- sum(abs(projected)*budgets)
    triangle_budget <- sum(abs(coefficients)*transformed)
    simple_budget <- envelope(pool)
    check(actual_budget<=triangle_budget+1e-12 && triangle_budget<=simple_budget+1e-12 &&
      simple_budget<=corner_bound+1e-12,paste(label,"four-corner bias envelope"))
    maximum_envelope_ratio <- max(maximum_envelope_ratio,actual_budget/corner_bound)
    errors <- list(sign(projected)*budgets,
      drop(design%*%c(min(budgets)/2,-min(budgets)/3)))
    if(length(outside)) {
      adverse <- sign(projected)*budgets
      adverse[outside] <- rep(c(-.4,.6),length.out=length(outside))
      errors[[3L]] <- adverse
    }
    for(error in errors) {
      noncentrality <- if(nrow(calibration$bias_map))
        quadratic(drop(calibration$bias_map%*%error),calibration$bias_precision) else 0
      noise <- sqrt(max(0,variance-calibration$variance)*max(0,noncentrality))
      check(abs(sum(difference*error))<=noise+1e-10,paste(label,"projected Cauchy bound"))
      allowance <- min(reference_budget+sqrt(max(0,reference_factor*noncentrality)),
        corner_bound+sqrt(max(0,shared_factor*noncentrality)))
      check(abs(sum(calibration$coefficients*error))<=allowance+1e-10,paste(label,"paired bias/variance bound"))
    }
    raw_difference <- coefficients-calibration$coefficients
    blind_error <- drop(covariance%*%residual%*%covariance%*%raw_difference)
    if(max(abs(blind_error))>1e-15) {
      blind_error <- blind_error*.5*min(budgets)/max(abs(blind_error))
      lambda <- if(nrow(calibration$bias_map))
        quadratic(drop(calibration$bias_map%*%blind_error),calibration$bias_precision) else 0
      if(abs(sum(raw_difference*blind_error))>
         sqrt(max(0,(raw_variance-calibration$variance)*lambda))+1e-10) {
        unprojected_failures <- unprojected_failures+1L
      }
    }
  }
  data.frame(case=label,K=count,valid_mu1=counts[1L],valid_mu0=counts[2L],
    subsets=subset_count,corner_bound=corner_bound,reference_budget=reference_budget,
    max_exact_to_corner_ratio=maximum_envelope_ratio,max_row_space_error=maximum_projection_error,
    unprojected_cauchy_failures=unprojected_failures)
}

results <- list()
for(count in c(2L,4L,6L)) for(seed in 1:3) for(pattern in c("zero","positive","negative","mixed")) {
  set.seed(801L+seed)
  factor <- matrix(rnorm(16L,sd=.25),4L,4L)
  diag(factor) <- c(.8,1,.45,.55)
  common <- tcrossprod(factor)/500
  width <- count+1L
  map <- c(1L,rep(3L,count),2L,rep(4L,count))
  covariance <- common[map,map]
  private <- cbind(seq(.4,1.8,length.out=count),seq(2,.5,length.out=count))/500
  correlation <- switch(pattern,zero=rep(0,count),positive=rep(.65,count),negative=rep(-.65,count),
    mixed=rep(c(-.6,.5),length.out=count))
  for(site in seq_len(count)) {
    positions <- c(site+1L,width+site+1L)
    cross <- correlation[site]*sqrt(prod(private[site,]))
    covariance[positions,positions] <- covariance[positions,positions]+matrix(c(private[site,1L],cross,cross,private[site,2L]),2L)
  }
  for(counts in list(c(0L,0L),c(0L,count/2L),c(count/2L,0L),c(count/2L,count/2L),
                     c(count/2L,count),c(count,count/2L),c(count,count))) {
    label <- paste("synthetic",count,seed,pattern,paste(counts,collapse="_"),sep="_")
    results[[length(results)+1L]] <- audit_covariance(covariance,counts,label)
  }
}
table <- do.call(rbind,results)
check(sum(table$unprojected_cauchy_failures)>0L,"negative controls must expose the missing projection")
write.csv(table,file.path(output,"geometry_checks.csv"),row.names=FALSE)
writeLines(c("VALID_SCORE_BIAS_GEOMETRY_CHECKS_PASSED",paste("Assertions:",audit_checks),
  paste("Covariance/count cases:",nrow(table)),paste("Enumerated subsets:",sum(table$subsets)),
  "Known-covariance finite-matrix identities and externally supplied bias budgets only; no fitted-CI validity claim."),
  file.path(output,"checks.txt"))
cat("VALID_SCORE_BIAS_GEOMETRY_CHECKS_PASSED",audit_checks,"assertions\n")
