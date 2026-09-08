# Candidate adapter: choose a declared existing dictionary that spans both inputs.
# This does not infer new features or alter either initial nuisance fit.
.fit_shared_calibration_basis <- function(outcome_basis, site_basis, final_basis) {
  valid <- function(x) is.matrix(x) && is.numeric(x) && nrow(x)>0L && ncol(x)>0L && all(is.finite(x))
  if(!valid(outcome_basis) || !valid(site_basis) || nrow(outcome_basis)!=nrow(site_basis))
    stop("invalid initial nuisance basis matrices")
  if(!is.character(final_basis) || length(final_basis)!=1L || is.na(final_basis) ||
     !final_basis %in% c("outcome","site")) stop("final_basis must explicitly name outcome or site")
  U<-cbind(1,if(final_basis=="outcome") outcome_basis else site_basis)
  decomposition<-qr(U,tol=1e-12)
  if(decomposition$rank!=ncol(U)) stop("declared final basis is not full column rank on training data")
  map<-function(basis) {
    initial<-cbind(1,basis)
    coefficients<-qr.coef(decomposition,initial)
    residual<-max(abs(U%*%coefficients-initial))
    if(!all(is.finite(coefficients)) || residual>1e-10*max(1,max(abs(initial))))
      stop("declared final basis does not span an initial nuisance basis")
    coefficients
  }
  list(final_basis=final_basis,outcome_map=map(outcome_basis),site_map=map(site_basis),
       outcome_width=ncol(outcome_basis),site_width=ncol(site_basis))
}

.apply_shared_calibration_basis <- function(spec,outcome_basis,site_basis) {
  if(!is.matrix(outcome_basis) || !is.matrix(site_basis) || !is.numeric(outcome_basis) ||
    !is.numeric(site_basis) || ncol(outcome_basis)!=spec$outcome_width ||
    ncol(site_basis)!=spec$site_width || nrow(outcome_basis)!=nrow(site_basis) ||
    nrow(outcome_basis)<1L || any(!is.finite(outcome_basis)) || any(!is.finite(site_basis)))
    stop("new nuisance bases do not match the declared training specification")
  U<-if(spec$final_basis=="outcome") outcome_basis else site_basis
  augmented<-cbind(1,U)
  for(name in c("outcome","site")) {
    original<-cbind(1,if(name=="outcome")outcome_basis else site_basis)
    if(max(abs(augmented%*%spec[[paste0(name,"_map")]]-original))>1e-10*max(1,max(abs(original))))
      stop("training basis mapping does not hold on new observations")
  }
  U
}
