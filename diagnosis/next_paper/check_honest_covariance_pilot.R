#!/usr/bin/env Rscript
# Validate sample boundaries and site-specific covariance scaling using a saved fit.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=3L) stop("Usage: repository_root saved_case new_scratch_output")
repository <- normalizePath(arguments[1L],mustWork=TRUE)
directory <- normalizePath(arguments[2L],mustWork=TRUE)
output <- arguments[3L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
configuration <- readRDS(file.path(directory,"configuration.rds"))
.libPaths(c(configuration$library,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=configuration$library))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(configuration$library,"RoCE"))))
files <- c("candidate_score_summary.R","honest_candidate_scores.R","honest_covariance_pilot.R")
for(file in files) source(file.path(repository,"diagnosis/next_paper",file))
library(testthat)
set.seed(configuration$seed)
dgp <- generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,configuration$K),p=4L,
  config=configuration$scenario,ate_deviation=configuration$rho,
  n_deviated_sites=configuration$invalid_sources,deviation_mechanism="treated_arm")
data <- split_data_by_site(dgp)
saved <- readRDS(file.path(directory,"fitted.rds"))
stopifnot(identical(digest::digest(data,algo="sha256"),saved$data_hash),
  identical(digest::digest(configuration,algo="sha256"),saved$configuration_hash))
fitted <- saved$fitted
previous <- readRDS(file.path(directory,"population_packets.rds"))$variants$translated_source$packet
packet <- evaluate_honest_candidates(fitted,data)
split <- split_honest_evaluation(fitted)
pilot <- evaluate_honest_candidates(fitted,data,keep_scores=TRUE,evaluation_indices=split$covariance)
evaluation <- evaluate_honest_candidates(fitted,data,keep_scores=TRUE,evaluation_indices=split$evaluation)

test_that("default evaluation remains exactly compatible with the frozen packet", {
  expect_identical(packet,previous)
  expect_identical(rescale_honest_covariance(packet,packet$evaluation_sizes)$covariance,packet$covariance)
})

test_that("sample roles are exhaustive, disjoint and independent of random state", {
  rng <- .Random.seed
  expect_identical(split,split_honest_evaluation(fitted))
  expect_identical(.Random.seed,rng)
  for(site in fitted$sites) {
    expect_length(intersect(split$covariance[[site]],split$evaluation[[site]]),0L)
    expect_equal(sort(c(split$covariance[[site]],split$evaluation[[site]])),sort(fitted$partition[[site]]$evaluation))
    expect_length(intersect(c(split$covariance[[site]],split$evaluation[[site]]),fitted$partition[[site]]$training),0L)
    expect_length(split$covariance[[site]],250L)
    expect_length(split$evaluation[[site]],250L)
  }
  expect_identical(pilot$training_shifts,evaluation$training_shifts)
  reordered <- split$evaluation[rev(fitted$sites)]
  expect_identical(evaluate_honest_candidates(fitted,data,TRUE,reordered),evaluation)
  for(fraction in list(0,1,NA_real_,"half",c(.3,.5),.999)) {
    expect_error(split_honest_evaluation(fitted,fraction))
  }
})

test_that("malformed subsets and training leakage are rejected", {
  expect_error(evaluate_honest_candidates(fitted,data,evaluation_indices=unname(split$evaluation)),"named list")
  expect_error(evaluate_honest_candidates(fitted,data,evaluation_indices=split$evaluation[-1L]),"named list")
  duplicated_site <- split$evaluation
  names(duplicated_site)[2L] <- "t"
  expect_error(evaluate_honest_candidates(fitted,data,evaluation_indices=duplicated_site),"named list")
  for(rows in list(c(1L,1L),NA_real_,c(1.1,2),0L,1001L,integer(0),"1")) {
    invalid <- split$evaluation; invalid$t <- rows
    expect_error(evaluate_honest_candidates(fitted,data,evaluation_indices=invalid),"indices")
  }
  for(rows in list(fitted$partition$t$training[1:2],split$evaluation$t[1L])) {
    invalid <- split$evaluation; invalid$t <- rows
    expect_error(evaluate_honest_candidates(fitted,data,evaluation_indices=invalid),"held-out")
  }
})

mutate_rows <- function(rows) {
  changed <- data
  for(site in fitted$sites) {
    changed[[site]]$Y[rows[[site]]] <- 1L-changed[[site]]$Y[rows[[site]]]
    changed[[site]]$W_outcome[rows[[site]],] <- -changed[[site]]$W_outcome[rows[[site]],]
    changed[[site]]$Z_site[rows[[site]],] <- -changed[[site]]$Z_site[rows[[site]],]
  }
  changed
}
test_that("pilot and final scores use only their own observations", {
  changed_pilot <- mutate_rows(split$covariance)
  expect_identical(evaluate_honest_candidates(fitted,changed_pilot,TRUE,split$evaluation),evaluation)
  changed <- evaluate_honest_candidates(fitted,changed_pilot,TRUE,split$covariance)
  expect_gt(max(abs(changed$covariance-pilot$covariance)),1e-10)
  changed_evaluation <- mutate_rows(split$evaluation)
  expect_identical(evaluate_honest_candidates(fitted,changed_evaluation,TRUE,split$covariance),pilot)
  changed <- evaluate_honest_candidates(fitted,changed_evaluation,TRUE,split$evaluation)
  expect_gt(max(abs(changed$means-evaluation$means)),1e-10)
  changed_training <- data
  row <- fitted$partition$t$training[1L]
  changed_training$t$Y[row] <- 1L-changed_training$t$Y[row]
  expect_error(evaluate_honest_candidates(fitted,changed_training,TRUE,split$evaluation),"Training data changed")
})

test_that("unequal site sizes require separate target and private covariance scaling", {
  pilot_rows <- final_rows <- split$evaluation
  for(index in seq_along(fitted$sites)) {
    site <- fitted$sites[index]
    pilot_rows[[site]] <- head(split$covariance[[site]],100L+index)
    final_rows[[site]] <- head(split$evaluation[[site]],190L-index)
  }
  unequal <- evaluate_honest_candidates(fitted,data,TRUE,pilot_rows)
  final_sizes <- as.numeric(lengths(final_rows)); names(final_sizes) <- fitted$sites
  scaled <- rescale_honest_covariance(unequal,final_sizes)
  direct <- function(scores,size) {
    centered <- sweep(scores,2L,colMeans(scores),"-")
    crossprod(centered)/(nrow(scores)*size)
  }
  expected <- direct(unequal$scores$target,final_sizes[["t"]])
  sources <- setdiff(fitted$sites,"t"); width <- length(fitted$sites)
  for(index in seq_along(sources)) {
    site <- sources[index]; columns <- c(index+1L,width+index+1L)
    private <- direct(unequal$scores$source[[site]],final_sizes[[site]])
    expect_equal(scaled$source_private_covariances[[site]],private,tolerance=1e-14)
    expected[columns,columns] <- expected[columns,columns]+private
  }
  expect_equal(scaled$covariance,expected,tolerance=1e-14)
  expect_equal(scaled$common_covariance,
    direct(cbind(unequal$scores$anchor,unequal$scores$common),final_sizes[["t"]]),tolerance=1e-14)
  expect_identical(rescale_honest_covariance(unequal,rev(final_sizes)),scaled)
  expect_gt(max(abs(scaled$covariance-unequal$covariance*
    unequal$evaluation_sizes[["t"]]/final_sizes[["t"]])),1e-10)
  expect_error(rescale_honest_covariance(unequal,unname(final_sizes)),"sizes")
  expect_error(rescale_honest_covariance(unequal,final_sizes[-1L]),"sizes")
  invalid_sizes <- final_sizes; invalid_sizes[1L] <- 1.5
  expect_error(rescale_honest_covariance(unequal,invalid_sizes),"sizes")
  invalid <- unequal; invalid$source_private_covariances[[1L]][1L,2L] <- NA_real_
  expect_error(rescale_honest_covariance(invalid,final_sizes),"covariance block")
})

dir.create(output,recursive=TRUE)
saveRDS(list(split=split,configuration=configuration,
  source_hashes=tools::md5sum(file.path(repository,"diagnosis/next_paper",c(files,"check_honest_covariance_pilot.R")))),
  file.path(output,"validation.rds"))
writeLines(c("HONEST_COVARIANCE_PILOT_CHECKS_PASSED",
  "Frozen default packet exactly reproduced; disjoint roles, data-mutation boundaries and unequal-size scaling checked.",
  "No nuisance refits or production algorithm changes. These checks do not establish interval coverage."),
  file.path(output,"checks.txt"))
