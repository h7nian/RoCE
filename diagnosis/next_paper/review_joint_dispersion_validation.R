#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if (length(arguments)!=2L) stop("Usage: completed_gaussian_run new_scratch_report")
root <- arguments[[1L]]; output <- arguments[[2L]]
if (any(!startsWith(c(root,output),"/scratch.global/zhan9381/FACE-HD/")) ||
    !file.exists(file.path(root,"COMPLETE")) || dir.exists(output)) stop("Use a completed run and a new FACE-HD report directory")
dir.create(output,recursive=TRUE)
configuration <- readRDS(file.path(root,"configuration.rds"))
metrics <- read.csv(file.path(root,"metrics.csv"),stringsAsFactors=FALSE)
comparisons <- vector("list",nrow(configuration$settings))
for (cell in seq_len(nrow(configuration$settings))) {
  saved <- readRDS(file.path(root,"cells",paste0(cell,".rds")))
  stopifnot(saved$cell==cell,identical(saved$configuration,configuration))
  records <- saved$records
  for (method in unique(records$method)) {
    observed <- records[records$method==method,]
    reported <- metrics[metrics$cell==cell & metrics$method==method,]
    truth <- configuration$mu1-configuration$mu0
    stopifnot(nrow(observed)==configuration$repeats,nrow(reported)==1L,
      abs(mean(observed$lower<=truth & observed$upper>=truth)-reported$coverage)<1e-12,
      abs(mean(observed$upper-observed$lower)-reported$mean_length)<1e-12)
  }
  first <- records[records$method=="joint_dispersion",]
  second <- records[records$method=="armwise_known_reference",]
  stopifnot(identical(first$iteration,second$iteration))
  length_difference <- (first$upper-first$lower)-(second$upper-second$lower)
  coverage_difference <- as.numeric(first$lower<=truth & first$upper>=truth)-as.numeric(second$lower<=truth & second$upper>=truth)
  reference_length <- metrics$mean_length[metrics$cell==cell & metrics$method=="guaranteed_gls"]
  comparisons[[cell]] <- data.frame(cell=cell,configuration$settings[cell,],
    paired_length_difference=mean(length_difference),paired_length_mcse=sd(length_difference)/sqrt(nrow(first)),
    paired_length_ratio_difference=mean(length_difference)/reference_length,
    paired_coverage_difference=mean(coverage_difference),paired_coverage_mcse=sd(coverage_difference)/sqrt(nrow(first)),
    reference_fallback_rate=mean(saved$joint_details$reference_fallback),row.names=NULL)
}
comparisons <- do.call(rbind,comparisons)
write.csv(comparisons,file.path(output,"paired_comparisons.csv"),row.names=FALSE)
write.csv(metrics,file.path(output,"metrics.csv"),row.names=FALSE)

colors <- c(joint_dispersion="#0072B2",armwise_known_reference="#D55E00",pooled_gls="#888888",guaranteed_gls="#222222",oracle_gls="#009E73")
labels <- c(joint_dispersion="Joint bias bound",armwise_known_reference="Armwise bias bounds",pooled_gls="Naive pooling",guaranteed_gls="Guaranteed-valid reference",oracle_gls="Oracle valid subset")
select_rows <- function(profile, count=NULL, bias=NULL) {
  rows <- metrics[metrics$shared_tate_variance==.004 & metrics$profile==profile,]
  if (profile=="majority_guarantee") {
    shared <- metrics[metrics$shared_tate_variance==.004 & metrics$profile=="three_quarters" & metrics$num_sources==4L,]
    rows <- rbind(rows,shared)
  }
  if (profile %in% c("treated_only","overstated_validity")) {
    shared_profile <- if (profile=="treated_only") "opposite_halves" else "three_quarters"
    rows <- rbind(rows,metrics[metrics$shared_tate_variance==.004 & metrics$profile==shared_profile & metrics$local_bias==0,])
  }
  if (!is.null(count)) rows <- rows[rows$num_sources==count,]
  if (!is.null(bias)) rows <- rows[rows$local_bias==bias,]
  rows
}

draw_figure <- function(filename, profiles, titles, boundary=FALSE) {
  pdf(file.path(output,filename),width=10.5,height=7.1)
  layout(matrix(c(1,2,3,4,5,6,7,7,7),3,3,byrow=TRUE),heights=c(1,1,.18))
  par(mar=c(3.8,4.2,2.3,.8),oma=c(0,0,2.8,0),mgp=c(2.3,.7,0),cex=.86)
  for (measure in c("coverage","length_to_known_reference")) for (index in seq_along(profiles)) {
    data <- select_rows(profiles[index],count=if (boundary) 64L else NULL,bias=if (boundary) NULL else .5)
    axis_name <- if (boundary) "local_bias" else "num_sources"
    limits <- if (boundary) c(0,4) else c(3.5,73)
    plot(NA,xlim=limits,ylim=if(measure=="coverage") c(0,1.01) else c(0,1.25),
      log=if(boundary) "" else "x",xaxt="n",xlab=if(boundary) "h (mean shift = h / sqrt(1000))" else "Number of sources K",
      ylab=if(measure=="coverage") "Coverage" else "Length / guaranteed-valid reference",main=titles[index])
    axis(1,at=if(boundary) c(0,.5,2,4) else c(4,16,64),labels=if(boundary) c("0","0.5","2","4") else c("4","16","64"))
    abline(h=if(measure=="coverage") .95 else 1,lty=3,col="#999999")
    methods <- if(measure=="coverage") names(colors)[1:4] else names(colors)
    for (method in methods) {
      rows <- data[data$method==method,]
      rows <- rows[order(rows[[axis_name]]),]
      x <- rows[[axis_name]]
      lines(x,rows[[measure]],type="b",pch=if(method=="joint_dispersion") 16 else 1,
        col=colors[[method]],lty=if(method %in% c("guaranteed_gls","oracle_gls")) 2 else 1,lwd=1.5)
      if(measure=="coverage" && method=="joint_dispersion") segments(x,rows$coverage_lower,x,rows$coverage_upper,col=colors[[method]])
    }
  }
  mtext(if(boundary) "Known-covariance Gaussian boundaries: K=64" else "Known-covariance Gaussian reference: weak shifts h=0.5",
        side=3,outer=TRUE,line=1.6,font=2,cex=1.1)
  mtext("1000/site; shared TATE variance = 0.004/1000; 1000 paired draws per setting",side=3,outer=TRUE,line=.1,cex=.85)
  par(mar=c(0,0,0,0))
  plot.new(); legend("center",legend=labels,col=colors,lty=c(1,1,1,2,2),pch=c(16,1,1,1,1),horiz=TRUE,bty="n",cex=.82)
  dev.off()
}
draw_figure("weak_bias_growth.pdf",c("treated_half","opposite_halves","majority_guarantee"),
            c("Half treated; all control valid","Half valid in each arm","Only a majority guaranteed"))
draw_figure("bias_boundaries.pdf",c("three_quarters","treated_only","overstated_validity"),
            c("Three-quarters valid per arm","Only treated means shift","Wrong validity bound (h > 0)"),TRUE)

selected <- metrics[metrics$method=="joint_dispersion" & metrics$num_sources==64 &
                    metrics$shared_tate_variance==.004 & metrics$local_bias==.5,]
valid <- metrics[metrics$method=="joint_dispersion" & metrics$assumption_valid,]
invalid <- metrics[metrics$method=="joint_dispersion" & !metrics$assumption_valid,]
lines <- c("# Joint Gaussian dispersion reference", "",
  sprintf("Completed %d settings with %d paired draws per setting; matched noise is reused across settings.",nrow(configuration$settings),configuration$repeats),
  "This is a known-covariance reference, not fitted high-dimensional RoCE validation. No scalar/source count is inferred from observed agreement.","",
  sprintf("Under the declared count assumptions, coverage ranges %.1f–%.1f%% and mean length is %.3f–%.3f times the equally informed guaranteed-valid reference.",
    100*min(valid$coverage),100*max(valid$coverage),min(valid$length_to_known_reference),max(valid$length_to_known_reference)),
  sprintf("With overstated valid counts, coverage ranges %.1f–%.1f%%. These failures are retained as assumption-boundary evidence.",100*min(invalid$coverage),100*max(invalid$coverage)),"",
  "| K=64, h=0.5, small shared variance | Coverage | Length/reference | Count assumption holds |",
  "| --- | ---: | ---: | --- |")
for (i in seq_len(nrow(selected))) lines <- c(lines,sprintf("| %s | %.1f%% | %.3f | %s |",selected$profile[i],100*selected$coverage[i],selected$length_to_known_reference[i],selected$assumption_valid[i]))
lines <- c(lines,"", "The K=4 majority-guarantee curve reuses the identical three-quarter setting; null curves reuse only configurations with identical declared counts and data.",
  "Coverage ranges are descriptive, not independent replications across settings. Paired length and coverage differences with MCSE are in paired_comparisons.csv.",
  "The joint and armwise methods use the same known-valid reference and error budget. The older ordinary-target-reference armwise variant remains in metrics.csv.",
  "Some small-K cells remain longer than the reference; gains are not uniform. Oracle length is reported separately and has not been uniformly attained.",
  "Estimated covariance, fitted nuisance remainders and a growing-K Gaussian reduction remain to be established before using this as the final RoCE method.")
writeLines(lines,file.path(output,"report.md"))
cat("JOINT_DISPERSION_REVIEW_COMPLETE\n")
