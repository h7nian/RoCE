#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(2L,3L)) stop("Usage: scalar_dispersion_file joint_dispersion_file [fitted_fixture_root]")
fitted_fixture_root <- if(length(arguments)==3L) normalizePath(arguments[[3L]],mustWork=TRUE) else NULL
source(arguments[[1L]])
source(arguments[[2L]])
library(testthat)

make_structured_test_covariance <- function(count, seed) {
  set.seed(seed)
  basis <- matrix(rnorm(24), 4, 6)
  common <- tcrossprod(basis) / 6 + diag(.2, 4)
  loading <- matrix(0, 2L * (count + 1L), 4L)
  loading[1L, 1L] <- loading[count + 2L, 2L] <- 1
  loading[2:(count + 1L), 3L] <- loading[(count + 3L):nrow(loading), 4L] <- 1
  covariance <- loading %*% common %*% t(loading)
  sources <- setdiff(seq_len(nrow(loading)), c(1L, count + 2L))
  diag(covariance)[sources] <- diag(covariance)[sources] + runif(2L * count, .4, 1.8)
  covariance
}

test_that("zero valid-score budgets retain the exact existing interval program", {
  covariance <- make_structured_test_covariance(4L,9631L)
  estimates <- matrix(seq(-.4,.5,length.out=30L),3L)
  for(method in c("subset_exact","shared_prediction_exact","shared_prediction_bound","known_valid_upper")) {
    original <- make_joint_dispersion_calibration(covariance,2L,4L,method)
    zero <- make_joint_dispersion_calibration(covariance,2L,4L,method,valid_score_bias=0)
    zero_matrix <- make_joint_dispersion_calibration(covariance,2L,4L,method,
      valid_score_bias=matrix(0,5L,2L,dimnames=list(NULL,c("mu1","mu0"))))
    expect_identical(original,zero)
    expect_identical(original,zero_matrix)
    expect_identical(joint_dispersion_intervals(estimates,original),joint_dispersion_intervals(estimates,zero))
  }
  for(invalid in list(-.01,NA,Inf,"0.01",c(.01,.02),matrix(.01,2L,5L))) {
    expect_error(make_joint_dispersion_calibration(covariance,2L,4L,valid_score_bias=invalid),"valid_score_bias")
  }
})

test_that("mean-error allowances bound arbitrary invalid-coordinate errors in every small subset", {
  for(count in c(2L,4L,6L)) for(sign in c(-1,1)) {
    covariance <- make_structured_test_covariance(count,9640L+count)/500
    private <- shared_prediction_variances(covariance)$private_variances
    width <- count+1L
    for(site in seq_len(count)) {
      positions <- c(site+1L,width+site+1L)
      covariance[positions[1L],positions[2L]] <- covariance[positions[1L],positions[2L]]+
        sign*.5*sqrt(prod(private[site,]))
      covariance[positions[2L],positions[1L]] <- covariance[positions[1L],positions[2L]]
    }
    budgets <- seq(.002,.013,length.out=2L*width)
    for(counts in list(c(count%/%2L,count),c(count%/%2L,count%/%2L),c(0L,count%/%2L))) {
      calibrations <- lapply(c("subset_exact","shared_prediction_bound","known_valid_upper"),function(method)
        make_joint_dispersion_calibration(covariance,counts[1L],counts[2L],method,valid_score_bias=budgets))
      reference <- make_joint_dispersion_calibration(covariance,counts[1L],counts[2L],"shared_prediction_bound")
      for(calibration in calibrations) {
        expect_identical(calibration$coefficients,reference$coefficients)
        expect_identical(calibration$bias_map,reference$bias_map)
        expect_identical(calibration$bias_precision,reference$bias_precision)
        expect_equal(calibration$reference_bias_allowance,
          sum(abs(calibration$reference_coefficients)*budgets),tolerance=1e-12)
      }
      sets <- lapply(counts,function(number) if(number==0L) list(integer()) else combn(count,number,simplify=FALSE))
      maximum_violation <- rep(-Inf,length(calibrations))
      for(first in sets[[1L]]) for(second in sets[[2L]]) {
        valid <- c(1L,width+1L,first+1L,second+width+1L)
        error <- rep(c(-.4,.3),length.out=2L*width)
        error[valid] <- budgets[valid]*sign(reference$coefficients[valid])
        contrast <- drop(reference$bias_map%*%error)
        noncentrality <- drop(crossprod(contrast,reference$bias_precision%*%contrast))
        bias <- abs(sum(reference$coefficients*error))
        for(index in seq_along(calibrations)) {
          allowance <- joint_score_bias_allowance(calibrations[[index]],noncentrality)
          maximum_violation[index] <- max(maximum_violation[index],bias-allowance)
        }
      }
      expect_true(all(maximum_violation<1e-10))
      common_shift <- drop(reference$design%*%c(min(budgets)/2,-min(budgets)/3))
      expect_lt(max(abs(reference$bias_map%*%common_shift)),1e-12)
      for(calibration in calibrations) {
        expect_gte(joint_score_bias_allowance(calibration,0),abs(sum(calibration$coefficients*common_shift))-1e-12)
      }
    }
  }
})

test_that("all-valid and anchor-only intervals include their reference mean-error allowance", {
  covariance <- make_structured_test_covariance(3L,9651L)/500
  estimates <- matrix(seq(.1,.8,length.out=8L),nrow=1L)
  for(counts in list(c(3L,3L),c(0L,0L))) {
    plain <- make_joint_dispersion_calibration(covariance,counts[1L],counts[2L])
    budgeted <- make_joint_dispersion_calibration(covariance,counts[1L],counts[2L],valid_score_bias=.01)
    old <- joint_dispersion_intervals(estimates,plain)
    new <- joint_dispersion_intervals(estimates,budgeted)
    expect_equal(old$lower-new$lower,budgeted$reference_bias_allowance,tolerance=1e-12)
    expect_equal(new$upper-old$upper,budgeted$reference_bias_allowance,tolerance=1e-12)
    expect_equal(new$bias_allowance,budgeted$reference_bias_allowance,tolerance=0)
    expect_identical(old$estimate,new$estimate)
  }
})

test_that("bias budgets enlarge the un-intersected radius without changing Q or GLS", {
  covariance <- make_structured_test_covariance(4L,9654L)/500
  plain <- make_joint_dispersion_calibration(covariance,2L,4L,"shared_prediction_bound")
  budgeted <- make_joint_dispersion_calibration(covariance,2L,4L,"shared_prediction_bound",valid_score_bias=.01)
  estimates <- matrix(seq(-.2,.3,length.out=30L),nrow=3L)
  old <- joint_dispersion_intervals(estimates,plain,reference_fraction=NULL)
  new <- joint_dispersion_intervals(estimates,budgeted,reference_fraction=NULL)
  expect_identical(old$estimate,new$estimate)
  expect_identical(old$statistic,new$statistic)
  expect_identical(old$noncentrality_upper,new$noncentrality_upper)
  expect_true(all(new$lower<=old$lower & new$upper>=old$upper))
  expect_true(all(new$bias_allowance>=old$bias_allowance))
})

test_that("private-correlation inflation bounds every small-K subset without changing the statistic", {
  for(count in c(2L,4L,6L)) for(sign in c(-1,1)) {
    covariance <- make_structured_test_covariance(count,9540L+count)
    private <- shared_prediction_variances(covariance)$private_variances
    correlations <- sign*seq(.15,.7,length.out=count)
    for(site in seq_len(count)) {
      positions <- c(site+1L,count+site+2L)
      value <- correlations[site]*sqrt(prod(private[site,]))
      covariance[positions[1L],positions[2L]] <- covariance[positions[1L],positions[2L]]+value
      covariance[positions[2L],positions[1L]] <- covariance[positions[1L],positions[2L]]
    }
    inflated <- inflate_shared_private_covariance(covariance)
    expect_gte(min(eigen(inflated$covariance-covariance,symmetric=TRUE,only.values=TRUE)$values),-1e-12)
    expect_equal(inflated$inflation,1+abs(correlations),tolerance=1e-12)
    for(validity in list(c(1,1),c(count%/%2,count),c(count,count%/%2),c(0,1),c(count,count),c(0,0))) {
      exact <- make_joint_dispersion_calibration(covariance,validity[1L],validity[2L],"subset_exact")
      bound <- make_joint_dispersion_calibration(covariance,validity[1L],validity[2L],"shared_prediction_bound")
      expect_gte(bound$bias_factor,exact$bias_factor-1e-10)
      expect_lte(bound$bias_factor,bound$reference_variance-bound$variance+1e-12)
      for(field in c("covariance","coefficients","variance","reference_coefficients","reference_variance","bias_map","bias_precision")) {
        expect_identical(bound[[field]],exact[[field]])
      }
      if(length(bound$uncertain)) {
        upper <- joint_mean_gls(inflated$covariance,bound$design,bound$bound_subset)$variance
        expect_equal(bound$bias_factor,min(upper,exact$reference_variance)-exact$variance,tolerance=1e-12)
      } else expect_equal(bound$bias_factor,0)
    }
  }
})

test_that("the new bound reduces to the exact zero-correlation calculation", {
  covariance <- make_structured_test_covariance(5L,9550L)
  exact <- make_joint_dispersion_calibration(covariance,2L,5L,"shared_prediction_exact")
  bound <- make_joint_dispersion_calibration(covariance,2L,5L,"shared_prediction_bound")
  expect_equal(bound$bias_factor,exact$bias_factor,tolerance=0)
  expect_identical(bound$coefficients,exact$coefficients)
  broken <- covariance
  broken[1L,2L] <- broken[2L,1L] <- broken[1L,2L]+.01
  expect_error(make_joint_dispersion_calibration(broken,2L,5L,"shared_prediction_bound"),"shared predictions")
})

test_that("actual fitted-score covariance packets satisfy the spectral subset bound", {
  skip_if(is.null(fitted_fixture_root),"Provide fitted_fixture_root for the saved-data integration check")
  root <- fitted_fixture_root
  for(version in c(3L,4L)) {
    checked <- readRDS(file.path(root,paste0("honest_candidate_checks_v",version),"checked_candidates.rds"))
    covariance <- checked$packet$covariance
    exact <- make_joint_dispersion_calibration(covariance,1L,2L,"subset_exact")
    bound <- make_joint_dispersion_calibration(covariance,1L,2L,"shared_prediction_bound")
    expect_gte(bound$bias_factor,exact$bias_factor-1e-14)
    expect_identical(bound$bias_map,exact$bias_map)
    expect_identical(bound$bias_precision,exact$bias_precision)
    expect_gte(min(bound$private_variance_inflation),1)
  }
})

test_that("joint bias estimates remove only known-valid residual noise", {
  covariance <- make_structured_test_covariance(4L, 9401)
  for (validity in list(c(2, 2), c(2, 4), c(4, 4), c(0, 0))) {
    calibration <- make_joint_dispersion_calibration(covariance, validity[1L], validity[2L])
    design <- calibration$design
    precision <- solve(covariance)
    null_projection <- precision - precision %*% design %*%
      solve(crossprod(design, precision %*% design)) %*% t(design) %*% precision
    expect_equal(unname(drop(crossprod(design, calibration$coefficients))), c(1, -1), tolerance = 1e-12)
    expect_equal(calibration$degrees_freedom, 8L - 4L*sum(validity == 4L))
    if (length(calibration$uncertain)) {
      expect_equal(unname(calibration$bias_map %*% design), matrix(0, length(calibration$uncertain), 2), tolerance = 1e-12)
      expect_equal(unname(calibration$bias_precision),
        unname(null_projection[calibration$uncertain, calibration$uncertain]), tolerance = 1e-12)
      expect_equal(drop(crossprod(calibration$coefficients, covariance %*% t(calibration$bias_map))),
                   rep(0, length(calibration$uncertain)), tolerance = 1e-12)
      set.seed(9402)
      values <- matrix(rnorm(30*nrow(covariance)), 30)
      known <- calibration$known_valid
      known_design <- design[known, , drop = FALSE]
      known_precision <- solve(covariance[known, known])
      known_projection <- known_precision - known_precision %*% known_design %*%
        solve(crossprod(known_design, known_precision %*% known_design)) %*% t(known_design) %*% known_precision
      direct <- rowSums((values %*% null_projection) * values) -
        rowSums((values[, known, drop = FALSE] %*% known_projection) * values[, known, drop = FALSE])
      contrasts <- values %*% t(calibration$bias_map)
      statistic <- rowSums((contrasts %*% calibration$bias_precision) * contrasts)
      expect_equal(statistic, direct, tolerance = 1e-11)
    }
  }
})

test_that("general covariance bias bound is sharp on a maximizing support", {
  set.seed(9403)
  basis <- matrix(rnorm(100), 10, 10)
  covariance <- tcrossprod(basis) + diag(1, 10)
  calibration <- make_joint_dispersion_calibration(covariance, 2L, 3L)
  precision <- solve(covariance)
  design <- calibration$design
  projection <- precision - precision %*% design %*% solve(crossprod(design, precision %*% design)) %*% t(design) %*% precision
  support <- setdiff(1:10, calibration$bound_subset)
  bias <- numeric(10)
  bias[support] <- solve(projection[support, support], calibration$coefficients[support])
  noncentrality <- drop(crossprod(bias, projection %*% bias))
  mean_bias <- sum(calibration$coefficients * bias)
  expect_equal(mean_bias^2, calibration$bias_factor*noncentrality, tolerance = 1e-10)
  upper <- make_joint_dispersion_calibration(covariance, 2L, 3L, "known_valid_upper")
  expect_gte(upper$bias_factor + 1e-12, calibration$bias_factor)
  expect_equal(calibration$subsets_examined, choose(4, 2)*choose(4, 3))
})

test_that("structured source selection matches exhaustive subset search", {
  for (seed in 9411:9414) for (validity in list(c(1, 2), c(2, 2), c(2, 4), c(0, 3))) {
    covariance <- make_structured_test_covariance(4L, seed)
    exact <- make_joint_dispersion_calibration(covariance, validity[1L], validity[2L], "subset_exact")
    fast <- make_joint_dispersion_calibration(covariance, validity[1L], validity[2L], "shared_prediction_exact")
    expect_equal(fast$bias_factor, exact$bias_factor, tolerance = 1e-12)
    expect_equal(fast$coefficients, exact$coefficients, tolerance = 1e-12)
    expect_equal(fast$subsets_examined, 1L)
    expect_lt(fast$maximum_structure_error, 1e-12)
  }
  covariance <- make_structured_test_covariance(4L, 9411)
  covariance[2L, 7L] <- covariance[7L, 2L] <- covariance[2L, 7L] + .02
  expect_error(make_joint_dispersion_calibration(covariance, 2, 2, "shared_prediction_exact"), "shared predictions")
  expect_error(make_joint_dispersion_calibration(diag(34), 8, 8, max_subsets = 100), "enumeration")
})

test_that("known-valid reference uses exactly the declared arm information", {
  covariance <- make_structured_test_covariance(4L, 9417)
  partial <- make_joint_dispersion_calibration(covariance, 2, 4)
  expect_setequal(partial$known_valid, c(1, 6:10))
  expect_equal(partial$reference_coefficients[2:5], rep(0, 4), tolerance = 1e-14)
  target_variance <- covariance[1, 1]+covariance[6, 6]-2*covariance[1, 6]
  expect_lte(partial$reference_variance, target_variance + 1e-12)
  values <- matrix(seq(-.2, .7, length.out = 30), 3)
  interval <- joint_dispersion_intervals(values, partial)
  reference <- drop(values %*% partial$reference_coefficients)
  radius <- qnorm(.9875)*sqrt(partial$reference_variance)
  expect_true(all(interval$upper-interval$lower <= 2*radius+1e-12))
  overlap <- !interval$empty_intersection
  expect_true(all(interval$lower[overlap] >= reference[overlap]-radius & interval$upper[overlap] <= reference[overlap]+radius))
  for (validity in list(c(0, 0), c(4, 4))) {
    calibration <- make_joint_dispersion_calibration(covariance, validity[1L], validity[2L])
    interval <- joint_dispersion_intervals(values, calibration)
    reference <- drop(values %*% calibration$reference_coefficients)
    radius <- qnorm(.975)*sqrt(calibration$reference_variance)
    expect_equal(interval$lower, reference-radius, tolerance = 1e-12)
    expect_equal(interval$upper, reference+radius, tolerance = 1e-12)
    expect_equal(interval$bias_allowance, rep(0, 3))
  }
})

test_that("shorter empty-intersection fallback retains both component length bounds", {
  shorter <- intersect_reference_interval(c(0,0,0),c(2,8,3),c(5,10,2),c(10,12,4))
  expect_equal(shorter$lower,c(0,10,2))
  expect_equal(shorter$upper,c(2,12,3))
  expect_identical(shorter$reference_fallback,c(FALSE,TRUE,FALSE))
  expect_identical(shorter$empty_intersection,c(TRUE,TRUE,FALSE))
  reference <- intersect_reference_interval(c(0,0,0),c(2,8,3),c(5,10,2),c(10,12,4),"reference")
  expect_equal(reference$lower,c(5,10,2))
  expect_equal(reference$upper,c(10,12,3))
})

test_that("the fast marginal comparator preserves exhaustive armwise intervals", {
  covariance <- make_structured_test_covariance(4L, 9421)
  values <- matrix(seq(-.3, .8, length.out = 30), 3)
  for (validity in list(c(2, 2), c(2, 4), c(3, 3))) {
    exact <- make_arm_dispersion_calibration(covariance, validity[1L], validity[2L], "subset_exact")
    fast <- make_shared_prediction_arm_calibration(covariance, validity[1L], validity[2L])
    expect_equal(fast$mu1$bias_factor, exact$mu1$bias_factor, tolerance = 1e-12)
    expect_equal(fast$mu0$bias_factor, exact$mu0$bias_factor, tolerance = 1e-12)
    expect_equal(arm_dispersion_intervals(values, fast), arm_dispersion_intervals(values, exact), tolerance = 1e-12)
  }
})

cat("JOINT_DISPERSION_TESTS_PASSED\n")
