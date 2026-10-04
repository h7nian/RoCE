#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: bound_source new_scratch_output")
source(arguments[[1L]])
output <- arguments[[2L]]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
dir.create(output,recursive=TRUE)
library(testthat)
test_that("finite binary likelihood mixtures match exhaustive observed-data probabilities", {
  outcome <- .6; treatment <- .5; control <- .4; probability <- .8; ratio <- probability/(1-probability)
  shift <- .02
  base <- c((1-treatment)*(1-control),(1-treatment)*control,treatment*(1-outcome),treatment*outcome)
  changed <- function(delta) c(base[1:2],treatment*(1-outcome-delta),treatment*(outcome+delta))
  for(n in 1:4) {
    observations <- as.matrix(expand.grid(rep(list(1:4),n)))
    mass <- apply(observations,1,function(row) prod(base[row]))
    first <- apply(observations,1,function(row) prod(changed(shift)[row]/base[row]))
    second <- apply(observations,1,function(row) prod(changed(-ratio*shift)[row]/base[row]))
    mixture <- probability*first+(1-probability)*second
    standardized <- shift/sqrt(outcome*(1-outcome)/(n*treatment))
    exact <- bernoulli_mixture_chisquare(n,standardized,probability)
    expect_equal(sum(mass*(mixture-1)^2),exact,tolerance=1e-12)
    expect_equal(sum(mass*first*second),(1-treatment*ratio*shift^2/(outcome*(1-outcome)))^n,tolerance=1e-12)
  }
  u <- .01
  expect_equal(bernoulli_mixture_chisquare(2,u,probability),ratio^2*u^4/4,tolerance=1e-14)
  expect_equal(bernoulli_mixture_chisquare(1000,u,probability)/u^4,(1-1/1000)*ratio^2/2,tolerance=.002)
})
settings <- expand.grid(num_sources=c(4L,16L,64L,256L,1024L),valid_fraction=c(.6,.75))
results <- do.call(rbind,lapply(seq_len(nrow(settings)),function(index)
  bernoulli_adaptation_lower_bound(1000L,settings$num_sources[index],settings$valid_fraction[index])))
stopifnot(all(results$expected_length_lower_bound>0))
write.csv(results,file.path(output,"finite_bounds.csv"),row.names=FALSE)
writeLines("Exact observed-data likelihood identities passed; finite positive lower bounds retained.",file.path(output,"COMPLETE"))
cat("BERNOULLI_ADAPTATION_BOUND_CHECKED\n")
