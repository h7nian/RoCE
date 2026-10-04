#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=4L) stop("Usage: scalar_source joint_source completed_v8_run new_report")
source(arguments[[1L]]); source(arguments[[2L]])
root <- arguments[[3L]]; output <- arguments[[4L]]
if(any(!startsWith(c(root,output),"/scratch.global/zhan9381/FACE-HD/")) ||
   !file.exists(file.path(root,"COMPLETE")) || dir.exists(output)) stop("Use a completed reference and new scratch output")
dir.create(output,recursive=TRUE); dir.create(file.path(output,"cells"))
configuration <- readRDS(file.path(root,"configuration.rds"))
if(!is.null(configuration$empty_intersection) && configuration$empty_intersection!="reference") stop("This comparison requires a reference-fallback baseline")
calibrations <- readRDS(file.path(root,"calibrations.rds"))
rows <- list()
for(cell in seq_len(nrow(configuration$settings))) {
  saved <- readRDS(file.path(root,"cells",paste0(cell,".rds")))
  calibration <- calibrations[[saved$calibration_key]]$joint
  detail <- saved$joint_details
  noise_alpha <- configuration$alpha*(1-configuration$reference_fraction)*(1-configuration$dispersion_fraction)
  radius <- qnorm(noise_alpha/2,lower.tail=FALSE)*sqrt(calibration$variance)+detail$bias_allowance
  lower <- detail$estimate-radius; upper <- detail$estimate+radius
  reference <- saved$records[saved$records$method=="guaranteed_gls",]
  center <- (reference$lower+reference$upper)/2
  half_width <- (reference$upper-reference$lower)/2*
    qnorm(configuration$alpha*configuration$reference_fraction/2,lower.tail=FALSE)/qnorm(configuration$alpha/2,lower.tail=FALSE)
  old <- intersect_reference_interval(lower,upper,center-half_width,center+half_width,"reference")
  error <- max(abs(old$lower-detail$lower),abs(old$upper-detail$upper))
  stopifnot(error<1e-12)
  new <- intersect_reference_interval(lower,upper,center-half_width,center+half_width,"shorter")
  lengths <- new$upper-new$lower
  old_lengths <- old$upper-old$lower
  stopifnot(all(lengths<=old_lengths+1e-12),all(lengths<=upper-lower+1e-12),all(lengths<=2*half_width+1e-12))
  truth <- configuration$mu1-configuration$mu0
  covered <- new$lower<=truth & new$upper>=truth
  ci <- binom.test(sum(covered),length(covered))$conf.int
  rows[[cell]] <- data.frame(cell=cell,configuration$settings[cell,],coverage=mean(covered),
    coverage_lower=ci[1L],coverage_upper=ci[2L],
    old_coverage=mean(old$lower<=truth & old$upper>=truth),
    length_to_known_reference=mean(lengths)/mean(reference$upper-reference$lower),
    old_length_to_known_reference=mean(old_lengths)/mean(reference$upper-reference$lower),
    paired_length_change=mean(lengths-old_lengths),paired_change_mcse=sd(lengths-old_lengths)/sqrt(length(lengths)),
    empty_rate=mean(new$empty_intersection),reference_fallback_rate=mean(new$reference_fallback),
    reconstruction_error=error,row.names=NULL)
  saveRDS(list(cell=cell,data_sha256=saved$data_sha256,intervals=new),file.path(output,"cells",paste0(cell,".rds")))
}
write.csv(do.call(rbind,rows),file.path(output,"metrics.csv"),row.names=FALSE)
writeLines("Same Gaussian draws; all old intervals reconstructed and both component length bounds verified.",file.path(output,"COMPLETE"))
cat("SHORTER_FALLBACK_CHECKS_PASSED\n")
