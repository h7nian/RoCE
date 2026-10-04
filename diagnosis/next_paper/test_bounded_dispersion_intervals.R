#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: joint_interval_file bounded_interval_file")
source(arguments[[1L]]); source(arguments[[2L]])
library(testthat)

block_matrix <- function(count,n) {
  site <- rep(seq_len(count),each=n)
  same <- outer(site,site,`==`)
  kernel <- ifelse(same,1/(count*n*(n-1)),-1/(count*(count-1)*n^2))
  diag(kernel) <- 0
  kernel
}

test_that("linear-cost summaries match the exact quadratic kernel and its norms", {
  set.seed(9901)
  for(count in c(2L,3L,5L)) for(n in c(2L,3L,5L)) {
    values <- matrix(runif(n*count),n,count)
    kernel <- block_matrix(count,n)
    computed <- bounded_site_dispersion(t(colSums(values)),t(colSums(values^2)),n,1)
    expect_equal(computed$statistic,drop(crossprod(as.vector(values),kernel %*% as.vector(values))),tolerance=1e-12)
    expect_equal(computed$frobenius_squared,sum(kernel^2),tolerance=1e-12)
    expect_equal(computed$operator_norm,max(abs(eigen(kernel,symmetric=TRUE,only.values=TRUE)$values)),tolerance=1e-12)
    expect_equal(rowSums(kernel),rep(0,n*count),tolerance=1e-12)
  }
})

test_that("population variance and one-sided inversion hold for exact Bernoulli laws", {
  count <- 3L; n <- 2L; means <- c(.2,.5,.9)
  kernel <- block_matrix(count,n)
  probabilities <- rep(means,each=n)
  observations <- as.matrix(expand.grid(rep(list(0:1),n*count)))
  mass <- apply(observations,1,function(row) prod(ifelse(row==1,probabilities,1-probabilities)))
  values <- rowSums((observations %*% kernel)*observations)
  truth <- var(means)
  variance <- probabilities*(1-probabilities)
  exact_variance <- 4*sum((kernel %*% probabilities)^2*variance)+2*sum(kernel^2*outer(variance,variance))
  expect_equal(sum(mass*values),truth,tolerance=1e-12)
  expect_equal(sum(mass*(values-truth)^2),exact_variance,tolerance=1e-12)
  summaries <- bounded_site_dispersion(matrix(0,1,count),matrix(0,1,count),n,1)
  expect_gte(summaries$linear_variance_factor*truth+summaries$null_variance_bound,exact_variance-1e-12)
  summaries$statistic <- values
  for(delta in c(.1,.25,.5)) {
    upper <- cantelli_dispersion_upper(summaries,delta)
    expect_lte(sum(mass[upper<truth-1e-12]),delta)
  }
})

test_that("interval special cases and finite range requirements are explicit", {
  target <- c(10,12); sums <- rbind(c(9,10,12),c(8,11,13)); squares <- sums
  none <- bounded_source_interval(target,20,sums,squares,20,1,0)
  all <- bounded_source_interval(target,20,sums,squares,20,1,3)
  expect_equal(none$estimate,target/20)
  expect_equal(all$estimate,(target+rowSums(sums))/80)
  expect_equal(all$upper-all$lower,rep(2*sqrt(log(40)/160),2))
  partial <- bounded_source_interval(target,20,sums,squares,20,1,2)
  expect_true(all(partial$lower<=partial$upper))
  expect_true(all(partial$upper-partial$lower<=2*sqrt(log(80)/40)+1e-12))
  expect_error(bounded_site_dispersion(sums,0*sums,20,1),"Inconsistent")
})

test_that("arm contrasts preserve their error allocation without assuming arm independence", {
  target <- cbind(mu1=c(10,12),mu0=c(6,8))
  sums <- list(mu1=rbind(c(9,10,12),c(8,11,13)),mu0=rbind(c(5,7,6),c(7,8,4)))
  result <- bounded_tate_interval(target,20,sums,sums,20,c(mu1=1,mu0=1),c(mu1=3,mu0=3))
  expect_equal(result$estimate,(target[,1]+rowSums(sums$mu1)-target[,2]-rowSums(sums$mu0))/80)
  expect_equal(result$lower,result$arm_intervals$mu1$lower-result$arm_intervals$mu0$upper)
  expect_equal(result$upper,result$arm_intervals$mu1$upper-result$arm_intervals$mu0$lower)
  expect_equal(result$upper-result$lower,rep(4*sqrt(log(80)/160),2))
  expect_error(bounded_tate_interval(target[,2:1],20,sums,sums,20,c(mu1=1,mu0=1),c(mu1=3,mu0=3)),"aligned")
})
cat("BOUNDED_DISPERSION_TESTS_PASSED\n")
