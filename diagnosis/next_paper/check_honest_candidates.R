#!/usr/bin/env Rscript
# Implementation checks for a research interface, not a coverage study.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)<3L || length(arguments)>6L) stop("Usage: checked_R_library repository_root new_scratch_output [C1|C2|C3] [source_count] [one_round|two_round]")
library_path <- normalizePath(arguments[[1L]],mustWork=TRUE)
repository <- normalizePath(arguments[[2L]],mustWork=TRUE)
output <- arguments[[3L]]
scenario <- if(length(arguments)>=4L) match.arg(arguments[[4L]],c("C1","C2","C3")) else "C2"
source_count <- if(length(arguments)>=5L) suppressWarnings(as.numeric(arguments[[5L]])) else 1L
protocol <- if(length(arguments)>=6L) match.arg(arguments[[6L]],c("one_round","two_round")) else "one_round"
if(!is.finite(source_count) || source_count<1L || source_count!=floor(source_count) ||
   source_count>.Machine$integer.max) stop("source_count must be a positive integer")
source_count <- as.integer(source_count)
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
dir.create(output,recursive=TRUE)
.libPaths(c(library_path,.libPaths()))
suppressPackageStartupMessages(library(RoCE,lib.loc=library_path))
stopifnot(identical(normalizePath(find.package("RoCE")),normalizePath(file.path(library_path,"RoCE"))))
source(file.path(repository,"diagnosis/next_paper/candidate_score_summary.R"))
source(file.path(repository,"diagnosis/next_paper/honest_candidate_scores.R"))
library(testthat)
set.seed(49231)
data <- split_data_by_site(generate_bounded_data(n_target=1000L,n_source_sizes=rep(1000L,source_count),
  K=source_count,p=4L,config=scenario))
stopifnot(length(data)==source_count+1L)
fitted <- fit_honest_candidates(data,calibration_folds=3L,nlambda=100L,
  communication_mode=protocol,
  checkpoint_dir=file.path(output,"checkpoints"))
packet <- evaluate_honest_candidates(fitted,data,keep_scores=TRUE)

test_that("held-out data do not change the split or any training input", {
  changed <- data
  for(site in fitted$sites) {
    rows <- fitted$partition[[site]]$evaluation
    changed[[site]]$A[rows] <- 1L-changed[[site]]$A[rows]
    changed[[site]]$Y[rows] <- 1L-changed[[site]]$Y[rows]
    changed[[site]]$W_outcome[rows,] <- -changed[[site]]$W_outcome[rows,]
    changed[[site]]$Z_site[rows,] <- -changed[[site]]$Z_site[rows,]
  }
  refitted <- fit_honest_candidates(changed,calibration_folds=3L,nlambda=100L,
    communication_mode=protocol,checkpoint_dir=file.path(output,"mutated_evaluation_checkpoints"))
  expect_identical(fitted$partition,refitted$partition)
  expect_identical(fitted$training_hashes,refitted$training_hashes)
  for(arm in c("mu1","mu0")) {
    expect_equal(as.numeric(fitted$models[[arm]]$common$coefficients),
                 as.numeric(refitted$models[[arm]]$common$coefficients),tolerance=0)
    for(type in c("outcome","weight")) {
      expect_equal(as.numeric(fitted$models[[arm]]$target[[type]]),
                   as.numeric(refitted$models[[arm]]$target[[type]]),tolerance=0)
      for(site in names(fitted$models[[arm]]$source)) {
        expect_equal(as.numeric(fitted$models[[arm]]$source[[site]][[type]]),
                     as.numeric(refitted$models[[arm]]$source[[site]][[type]]),tolerance=0)
      }
    }
  }
  changed_packet <- evaluate_honest_candidates(fitted,changed)
  expect_equal(packet$training_shifts,changed_packet$training_shifts,tolerance=0)
  changed$t$Y[fitted$partition$t$training[1L]] <- 1L-changed$t$Y[fitted$partition$t$training[1L]]
  expect_error(evaluate_honest_candidates(fitted,changed),"Training data changed")
})

gradient_checks <- list()
test_that("common projection solves full-target IPW loss, with both arms checked", {
  train <- honest_subset(data$t,fitted$partition$t$training)
  design <- cbind(1,train$W_outcome)
  for(arm in c("mu1","mu0")) {
    model <- fitted$models[[arm]]$common
    coefficients <- as.numeric(model$coefficients)
    weights <- model$inverse_propensity
    arm_value <- as.integer(arm=="mu1")
    expect_identical(as.integer(model$training_indices),as.integer(train$original_idx))
    gradient <- drop(crossprod(design,as.numeric(train$A==arm_value)*weights*
      (plogis(drop(design %*% coefficients))-train$Y)))/train$n
    lambda <- as.numeric(attr(model$coefficients,"lambda_used"))
    expect_length(lambda,1L)
    residual <- ifelse(abs(coefficients[-1L])>1e-7,
      abs(gradient[-1L]+lambda*sign(coefficients[-1L])),pmax(abs(gradient[-1L])-lambda,0))
    kkt <- max(abs(gradient[1L]),residual)
    expect_lt(kkt,1e-7)
    gradient_checks[[arm]] <<- data.frame(arm=arm,lambda=lambda,kkt_residual=kkt,
      minimum_weight=min(weights),maximum_weight=max(weights))
  }
})

test_that("source messages and covariance reconstruct individual score statistics", {
  expect_lt(packet$shift_identity_error,1e-12)
  expected_sizes <- setNames(rep(500,source_count+1L),fitted$sites)
  expect_equal(packet$evaluation_sizes,expected_sizes)
  expect_equal(packet$training_sizes,expected_sizes)
  target <- packet$scores$target
  centered <- function(x) sweep(x,2L,colMeans(x),"-")
  independent <- crossprod(centered(target))/nrow(target)^2
  for(index in seq_len(source_count)) {
    site <- names(packet$scores$source)[index]
    source <- packet$scores$source[[site]]
    columns <- c(index+1L,source_count+index+2L)
    independent[columns,columns] <- independent[columns,columns]+
      crossprod(centered(source))/nrow(source)^2
    message <- packet$source_messages[[site]]
    expect_equal(as.numeric(message$adjusted_sum),as.numeric(colSums(source)),tolerance=1e-12)
    expect_equal(as.numeric(message$adjusted_square_sum),as.numeric(colSums(source^2)),tolerance=1e-12)
    expect_equal(as.numeric(message$adjusted_cross_sum),sum(source[,1L]*source[,2L]),tolerance=1e-12)
    expect_equal(as.numeric(packet$means[site,]),as.numeric(colMeans(packet$scores$common)+colMeans(source)),tolerance=1e-12)
  }
  expect_equal(packet$covariance,independent,tolerance=1e-14)
  expect_gte(min(eigen(packet$covariance,symmetric=TRUE,only.values=TRUE)$values),-1e-14)
})

test_that("auxiliary proximal Newton and coordinate descent agree at the same penalty", {
  train <- honest_subset(data$t,fitted$partition$t$training)
  design <- cbind(1,train$W_outcome)
  fit_coordinate_descent <- function(model,arm_value,lambda) {
    previous <- RoCE:::.set_nuisance_solver("coordinate_descent")
    on.exit(RoCE:::.restore_nuisance_solver(previous))
    RoCE::fit_unified_outcome(train$W_outcome,train$Y,train$A,A_val=arm_value,
      gamma_s=c(0,1),Z_site=matrix(-log(model$inverse_propensity),ncol=1L),
      lambda=lambda,tol=1e-10,calibrated=FALSE,M_tau=12)
  }
  for(arm in c("mu1","mu0")) {
    model <- fitted$models[[arm]]$common
    arm_value <- as.integer(arm=="mu1")
    lambda <- as.numeric(attr(model$coefficients,"lambda_used"))
    coefficients <- fit_coordinate_descent(model,arm_value,lambda)
    objective <- function(beta) {
      predictor <- drop(design %*% beta)
      loss <- pmax(predictor,0)+log1p(exp(-abs(predictor)))-train$Y*predictor
      mean(as.numeric(train$A==arm_value)*model$inverse_propensity*loss)+lambda*sum(abs(beta[-1L]))
    }
    expect_equal(as.numeric(coefficients),as.numeric(model$coefficients),tolerance=1e-6)
    expect_lt(abs(objective(coefficients)-objective(model$coefficients)),1e-9)
  }
})

saveRDS(list(fitted=fitted,packet=packet),file.path(output,"checked_candidates.rds"))
write.csv(do.call(rbind,gradient_checks),file.path(output,"projection_kkt.csv"),row.names=FALSE)
writeLines(c("HONEST_CANDIDATE_IMPLEMENTATION_CHECKS_PASSED",
  paste("Single",scenario,"dataset, n=1000/site, p=4,",source_count,"sources,",protocol,"3 calibration blocks, 100 nuisance lambdas; evaluation mutation uses an independent refit cache."),
  "These checks do not establish nuisance rates, conditional bias bounds, or confidence-interval coverage."),
  file.path(output,"checks.txt"))
