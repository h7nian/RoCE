library(testthat)
source(if(file.exists("source_working_basis.R")) "source_working_basis.R" else
  "diagnosis/tate_common_weight/source_working_basis.R")
basis_fixture <- function() {
  X <- matrix(seq_len(12)/10,4,3)
  site <- list(n=4L,X=X,X_dagger=X,A=c(0,1,0,1),Y=c(1,1,0,0),
    W_outcome=cbind(sweep(X,2,c(.1,.2,.3),"-"),X^2), Z_site=cbind(X,X^2),
    W_outcome_true=cbind(X,X^2),Z_site_true=cbind(X,X^2))
  list(t=site,s1=site)
}
test_that("working-basis changes leave observations and truth untouched", {
  x<-basis_fixture()
  expect_identical(.source_working_basis(x),x)
  for(config in c("C2","C3")) {
    y<-.source_working_basis(x,config)
    changed<-if(config=="C2") "W_outcome" else "Z_site"
    for(site in names(x)) {
      expect_identical(y[[site]][[changed]],x[[site]]$X)
      preserved<-setdiff(names(x[[site]]),changed)
      expect_identical(y[[site]][preserved],x[[site]][preserved])
    }
  }
  expect_identical(x,basis_fixture())
})
test_that("malformed and already reduced inputs are rejected", {
  x<-basis_fixture()
  expect_error(.source_working_basis(x,"C4"),"config")
  expect_error(.source_working_basis(x,NA_character_),"config")
  expect_error(.source_working_basis(.source_working_basis(x,"C3"),"C2"))
  expect_error(.source_working_basis(.source_working_basis(x,"C2"),"C3"))
  x$t$Z_site[1,1]<-99
  expect_error(.source_working_basis(x,"C3"))
})
