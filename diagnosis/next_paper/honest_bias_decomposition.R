# Decompose valid-coordinate population drift; validity and moments are oracle
# inputs in the fitted study. Under a known-covariance Gaussian score model,
# residual_mahalanobis is the residual statistic's noncentrality parameter.

decompose_honest_bias <- function(means,covariance,truth,valid) {
  if(!is.matrix(means) || !is.numeric(means) || ncol(means)!=2L ||
     !is.matrix(valid) || !is.logical(valid) || !identical(dim(means),dim(valid)) ||
     anyNA(valid) || !all(valid[1L,]) || length(truth)!=2L ||
     any(!is.finite(c(means,truth,covariance)))) stop("Invalid candidate means or validity labels")
  width <- nrow(means)
  if(!is.matrix(covariance) || !identical(dim(covariance),rep(2L*width,2L)) ||
     max(abs(covariance-t(covariance)))>1e-12) stop("Covariance dimensions or symmetry differ")
  design <- cbind(mu1=rep(c(1,0),each=width),mu0=rep(c(0,1),each=width))
  selected <- which(as.numeric(valid)==1)
  selected_design <- design[selected,,drop=FALSE]
  precision <- chol2inv(chol(covariance[selected,selected,drop=FALSE]))
  information <- crossprod(selected_design,precision%*%selected_design)
  mean_covariance <- chol2inv(chol(information))
  drift <- as.numeric(sweep(means,2L,truth,"-"))[selected]
  mean_shift <- drop(mean_covariance%*%crossprod(selected_design,precision%*%drift))
  residual <- drift-drop(selected_design%*%mean_shift)
  residual_mahalanobis <- drop(crossprod(residual,precision%*%residual))
  level_mahalanobis <- drop(crossprod(mean_shift,information%*%mean_shift))
  total_mahalanobis <- drop(crossprod(drift,precision%*%drift))
  identity_error <- abs(total_mahalanobis-residual_mahalanobis-level_mahalanobis)
  if(identity_error>1e-10*max(1,total_mahalanobis)) stop("Bias projection identity failed")
  degrees_freedom <- length(selected)-2L
  contrast <- c(1,-1)
  variance <- drop(crossprod(contrast,mean_covariance%*%contrast))
  list(mean_shift=setNames(mean_shift,c("mu1","mu0")),tate_shift=mean_shift[1L]-mean_shift[2L],
    oracle_sd=sqrt(variance),residual_mahalanobis=residual_mahalanobis,
    residual_degrees_freedom=degrees_freedom,
    residual_over_sqrt_df=if(degrees_freedom>0L) residual_mahalanobis/sqrt(degrees_freedom) else NA_real_,
    level_mahalanobis=level_mahalanobis,total_mahalanobis=total_mahalanobis,
    identity_error=identity_error)
}
