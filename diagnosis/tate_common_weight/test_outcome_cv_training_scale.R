library(testthat)
test_that("outcome CV keeps source-population training scale for both arms and families", {
  x<-c(as.vector(rbind(-(1:31)/31,(1:31)/31)),as.vector(rbind(-(1:19)/19,(1:19)/19)))
  X<-matrix(x,ncol=1);A<-rep(c(1,0),c(62,38));grid<-c(.01,.02)
  y_binary<-unlist(lapply(c(31L,19L),function(n) {
    positive<-as.numeric(seq_len(n)%%3!=0)
    as.vector(rbind(1-positive,positive))
  }))
  for(arm in 0:1)for(family in 0:1)for(calibrated in c(FALSE,TRUE)) {
    Y<-if(family==0L)2*x+.3*sin(3*x)else y_binary
    ids<-which(A==arm);q<-length(ids)/length(A)
    folds<-rep(rep(1:5,length.out=length(ids)/2),each=2)
    expected<-vapply(grid,function(lambda) {
      sum(vapply(1:5,function(k) {
        train<-ids[folds!=k];val<-ids[folds==k]
        if(family==0L) {
          slope<-max(q*mean(x[train]*Y[train])-lambda,0)/(q*mean(x[train]^2))
          losses<-.5*(Y[val]-slope*x[val])^2
        } else {
          gradient<-function(beta)q*mean(x[train]*(plogis(beta*x[train])-Y[train]))+lambda
          slope<-if(gradient(0)>=0)0 else uniroot(gradient,c(0,20),tol=1e-12)$root
          eta<-slope*x[val]
          losses<-pmax(eta,0)+log1p(exp(-abs(eta)))-Y[val]*eta
        }
        sum(losses)/length(A)/5
      },numeric(1L)))
    },numeric(1L))
    observed<-if(calibrated)RoCE:::select_lambda_cv_calibrated_outcome_cpp(
      X,Y,A,c(0,0),grid,5L,10000L,1e-10,arm,5,X,family,family,folds)else
      RoCE:::select_lambda_cv_general_refined_outcome_cpp(X,Y,A,c(0,0),family,family,
        grid,5L,10000L,1e-10,arm,X,folds)
    expect_equal(as.numeric(observed$cv_scores),expected,tolerance=1e-7,
      info=paste("arm",arm,"family",family,"calibrated",calibrated))
    expect_equal(observed$invalid_fold_fits,0L)
  }
})
