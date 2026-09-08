# Primitive source-assisted arm score and final-parameter derivative rows.
.calibrated_site_score <- function(W,Z,A,Y,alpha,gamma,site,M=5) {
  stopifnot(site %in% c("target","source"),is.numeric(M),length(M)==1L,is.finite(M),M>0,
    nrow(W)==nrow(Z),ncol(W)==length(alpha),ncol(Z)==length(gamma),
    length(A)==nrow(W),length(Y)==nrow(W))
  prediction<-plogis(drop(W%*%alpha));log_weight<-drop(Z%*%gamma)
  if(site=="target") {
    score<-prediction
    gradient<-cbind(W*(prediction*(1-prediction)),matrix(0,nrow(Z),ncol(Z)))
  } else {
    weight<-exp(-pmax(-M,pmin(M,log_weight)))
    score<-A*weight*(Y-prediction)
    gradient<-cbind(-W*(A*weight*prediction*(1-prediction)),
      -Z*(score*(abs(log_weight)<M)))
  }
  stopifnot(all(is.finite(score)),all(is.finite(gradient)))
  list(score=score,gradient=gradient,n_clipped=if(site=="source")sum(abs(log_weight)>M)else 0L)
}
