#!/usr/bin/env Rscript
# Known-Gaussian shared-prediction reference; no fitted-RoCE inference claim.
arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(3L,4L)) stop("Usage: frozen_source_directory scratch_output repeats [shorter|reference]")
source_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
source_files <- file.path(source_directory, c("dispersion_bias_intervals.R", "joint_dispersion_intervals.R",
                                              "joint_dispersion_validation_settings.R"))
for (file in source_files) source(file)
output <- arguments[[2L]]
repeats <- as.integer(arguments[[3L]])
empty_intersection <- match.arg(if(length(arguments)==4L) arguments[[4L]] else "shorter",c("shorter","reference"))
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/") || !is.finite(repeats) || repeats < 100L || repeats > 10000L) {
  stop("Use FACE-HD scratch and 100–10000 repetitions")
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "cells"), showWarnings = FALSE)
settings <- make_joint_dispersion_validation_settings()
configuration <- list(settings = settings, repeats = repeats, n_per_site = 1000L, mu1 = .6, mu0 = .4,
  alpha = .05, dispersion_fraction = .5, reference_fraction = .5, empty_intersection = empty_intersection,
  code_md5 = tools::md5sum(c(source_files, sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))))
configuration_path <- file.path(output, "configuration.rds")
if (file.exists(configuration_path)) stopifnot(identical(readRDS(configuration_path), configuration)) else saveRDS(configuration, configuration_path)
calibrations <- list()
all_metrics <- list()
truth <- configuration$mu1 - configuration$mu0
for (cell in seq_len(nrow(settings))) {
  path <- file.path(output, "cells", paste0(cell, ".rds"))
  if (file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(saved$cell == cell, identical(saved$configuration, configuration))
    if (is.null(calibrations[[saved$calibration_key]])) {
      calibrations[[saved$calibration_key]] <- list(
        joint = make_joint_dispersion_calibration(saved$covariance,saved$valid_mu1,saved$valid_mu0,"shared_prediction_exact"),
        armwise = make_shared_prediction_arm_calibration(saved$covariance,saved$valid_mu1,saved$valid_mu0))
    }
    all_metrics[[cell]] <- saved$metrics
    next
  }
  setting <- settings[cell, ]
  count <- setting$num_sources
  width <- count + 1L
  # An OR-correct limiting score has a common target-prediction component and
  # source-local residual variance, with zero cross-arm local residual covariance.
  prediction <- matrix(c(.02, (.04-setting$shared_tate_variance)/2,
                        (.04-setting$shared_tate_variance)/2, .02), 2L)
  covariance <- kronecker(prediction, matrix(1, width, width))
  private <- c(.5, .45*(1+.2*sin(seq_len(count))), .5, .45*(1+.2*cos(seq_len(count))))
  diag(covariance) <- diag(covariance) + private
  covariance <- covariance / configuration$n_per_site
  q1 <- setting$valid_mu1_minimum
  q0 <- setting$valid_mu0_minimum
  valid1 <- seq_len(setting$unshifted_mu1_count)
  valid0 <- seq.int(count-setting$unshifted_mu0_count+1L,count)
  bias <- numeric(2L*width)
  bias[1L+setdiff(seq_len(count),valid1)] <- setting$local_bias/sqrt(configuration$n_per_site)
  bias[width+1L+setdiff(seq_len(count),valid0)] <- -setting$local_bias/sqrt(configuration$n_per_site)
  population <- c(rep(configuration$mu1,width),rep(configuration$mu0,width)) + bias
  key <- paste(count, setting$shared_tate_variance, q1, q0)
  if (is.null(calibrations[[key]])) {
    calibrations[[key]] <- list(joint = make_joint_dispersion_calibration(covariance,q1,q0,"shared_prediction_exact"),
      armwise = make_shared_prediction_arm_calibration(covariance,q1,q0))
  }
  joint <- calibrations[[key]]$joint
  armwise <- calibrations[[key]]$armwise
  # Match the underlying normal draws across all bias and validity settings.
  set.seed(88400L+count)
  draws <- sweep(matrix(rnorm(repeats*2L*width),repeats) %*% chol(covariance),2L,population,"+")
  known <- drop(draws %*% joint$reference_coefficients)
  known_radius <- qnorm(configuration$alpha/2,lower.tail=FALSE)*sqrt(joint$reference_variance)
  target <- draws[,1L]-draws[,width+1L]
  target_variance <- covariance[1L,1L]+covariance[width+1L,width+1L]-2*covariance[1L,width+1L]
  target_radius <- qnorm(.975)*sqrt(target_variance)
  valid <- which(bias == 0)
  oracle <- joint_mean_gls(covariance,joint$design,valid)
  oracle_mean <- drop(draws %*% oracle$coefficients)
  pooled <- drop(draws %*% joint$coefficients)
  intervals <- list(target=data.frame(lower=target-target_radius,upper=target+target_radius),
    guaranteed_gls=data.frame(lower=known-known_radius,upper=known+known_radius),
    oracle_gls=data.frame(lower=oracle_mean-qnorm(.975)*sqrt(oracle$variance),upper=oracle_mean+qnorm(.975)*sqrt(oracle$variance)),
    pooled_gls=data.frame(lower=pooled-qnorm(.975)*sqrt(joint$variance),upper=pooled+qnorm(.975)*sqrt(joint$variance)))
  intervals$joint_dispersion <- joint_dispersion_intervals(draws,joint,alpha=configuration$alpha,
    dispersion_fraction=configuration$dispersion_fraction,reference_fraction=configuration$reference_fraction,
    empty_intersection=configuration$empty_intersection)
  intervals$joint_unanchored <- joint_dispersion_intervals(draws,joint,alpha=configuration$alpha,
    dispersion_fraction=configuration$dispersion_fraction,reference_fraction=NULL,
    empty_intersection=configuration$empty_intersection)
  intervals$armwise_target_reference <- arm_dispersion_intervals(draws,armwise,
    alpha=configuration$alpha,dispersion_fraction=configuration$dispersion_fraction,anchor_fraction=configuration$reference_fraction)
  arm_interval <- arm_dispersion_intervals(draws,armwise,
    alpha=configuration$alpha*(1-configuration$reference_fraction),dispersion_fraction=configuration$dispersion_fraction)
  radius <- qnorm(configuration$alpha*configuration$reference_fraction/2,lower.tail=FALSE)*sqrt(joint$reference_variance)
  combined <- intersect_reference_interval(arm_interval$lower,arm_interval$upper,known-radius,known+radius,
                                           configuration$empty_intersection)
  intervals$armwise_known_reference <- data.frame(lower=combined$lower,upper=combined$upper,
                                                 reference_fallback=combined$reference_fallback)
  records <- do.call(rbind,lapply(names(intervals),function(method) {
    interval <- intervals[[method]]
    data.frame(iteration=seq_len(repeats),method=method,lower=interval$lower,upper=interval$upper)
  }))
  stopifnot(all(is.finite(records$lower)),all(is.finite(records$upper)),all(records$lower<=records$upper))
  metrics <- do.call(rbind,lapply(split(records,records$method),function(rows) {
    covered <- rows$lower<=truth & rows$upper>=truth
    coverage_interval <- binom.test(sum(covered),repeats)$conf.int
    lengths <- rows$upper-rows$lower
    data.frame(cell=cell,setting,method=rows$method[1L],repeats=repeats,
      coverage=mean(covered),coverage_lower=coverage_interval[1L],coverage_upper=coverage_interval[2L],
      mean_length=mean(lengths),length_to_known_reference=mean(lengths)/(2*known_radius),
      length_to_target=mean(lengths)/(2*target_radius),length_to_oracle=mean(lengths)/(2*qnorm(.975)*sqrt(oracle$variance)),
      length_ratio_mcse=sd(lengths)/sqrt(repeats)/(2*known_radius),row.names=NULL)
  }))
  saved <- list(cell=cell,configuration=configuration,records=records,metrics=metrics,covariance=covariance,
    population_means=population,valid_mu1=q1,valid_mu0=q0,data_sha256=digest::digest(draws,algo="sha256"),
    joint_degrees_freedom=joint$degrees_freedom,joint_bias_factor=joint$bias_factor,
    joint_details=intervals$joint_dispersion,calibration_key=key)
  temporary <- paste0(path,".tmp")
  saveRDS(saved,temporary)
  stopifnot(file.rename(temporary,path))
  all_metrics[[cell]] <- metrics
  write.csv(do.call(rbind,all_metrics),file.path(output,"metrics.csv"),row.names=FALSE)
  cat("Completed joint-dispersion cell",cell,"of",nrow(settings),"at",format(Sys.time()),"\n")
}
saveRDS(calibrations,file.path(output,"calibrations.rds"))
write.csv(do.call(rbind,all_metrics),file.path(output,"metrics.csv"),row.names=FALSE)
writeLines("Prespecified shared-prediction Gaussian study complete. Estimated covariance and fitted nuisance inference remain unproved.",file.path(output,"COMPLETE"))
