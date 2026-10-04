#!/usr/bin/env Rscript
# Completed joint-TATE results with one treated-arm-deviated source.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)<2L || length(arguments)>3L) stop("Usage: completed_source_count_report new_scratch_figure_directory [p=100]")
report <- normalizePath(arguments[[1L]],mustWork=TRUE)
output <- arguments[[2L]]
dimension <- if(length(arguments)==3L) suppressWarnings(as.numeric(arguments[[3L]])) else 100L
if(!is.finite(dimension) || !dimension %in% c(100L,200L)) stop("p must be 100 or 200")
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use a new FACE-HD scratch directory")
metrics <- read.csv(file.path(report,"metrics.csv"),stringsAsFactors=FALSE)
metrics <- subset(metrics,method=="one_round_crossfit_ate_joint_tate" & p==dimension)
weights <- read.csv(file.path(report,"source_weights/weight_summary.csv"),stringsAsFactors=FALSE)
weights <- subset(weights,p==dimension)
stopifnot(nrow(metrics)==24L,nrow(weights)==24L,all(metrics$repeats==50L),
  all(metrics$p==dimension),all(metrics$n_deviated_sites==1L),
  !anyDuplicated(metrics[c("config","K","rho")]),!anyDuplicated(weights[c("config","K","rho")]))
dir.create(output,recursive=TRUE)
colors <- c("#0072B2","#D55E00"); symbols <- c(16,17)
source_counts <- c(4L,6L); scenarios <- c("C1","C2","C3")
shifts <- c(0,.5,1,2)
pdf(file.path(output,"source_count_boundary.pdf"),width=10,height=9,useDingbats=FALSE,
  title="Fixed source counts and outcome-departure boundary")
par(mfrow=c(3,3),mar=c(3.5,4,2,1),oma=c(3,0,3,0),mgp=c(2.5,.65,0))
for(quantity in c("coverage","bias","weight")) for(scenario in scenarios) {
  limits <- switch(quantity,coverage=c(.65,1.01),bias=c(-.01,.027),weight=c(-.005,.245))
  # A common scale within each row makes the scenario panels comparable.
  panel <- if(quantity=="weight") weights else metrics
  if(quantity=="coverage") {
    limits <- range(limits,panel$coverage_lower,panel$coverage_upper)
  } else {
    values <- if(quantity=="bias") panel$bias else panel$mu1_mean_weight
    errors <- if(quantity=="bias") panel$bias_mcse else panel$mu1_weight_mcse
    limits <- range(limits,values-1.96*errors,values+1.96*errors)
  }
  label <- switch(quantity,coverage="95% CI coverage",bias="TATE bias",weight="Treated-arm weight on source 1")
  plot(NA,xlim=c(-.06,2.06),ylim=limits,xlab=expression(rho),ylab=label,xaxt="n",main=scenario)
  axis(1,at=shifts,labels=shifts)
  abline(h=if(quantity=="coverage") .95 else 0,lty=2,col="gray45")
  for(index in 1:2) {
    rows <- subset(if(quantity=="weight") weights else metrics,config==scenario & K==source_counts[index])
    rows <- rows[match(shifts,rows$rho),]
    stopifnot(nrow(rows)==4L,all(is.finite(rows$rho)))
    if(quantity=="coverage") {
      y <- rows$coverage; lower <- rows$coverage_lower; upper <- rows$coverage_upper
    } else {
      y <- if(quantity=="bias") rows$bias else rows$mu1_mean_weight
      error <- if(quantity=="bias") rows$bias_mcse else rows$mu1_weight_mcse
      lower <- y-1.96*error; upper <- y+1.96*error
    }
    x <- rows$rho+c(-.012,.012)[index]
    nonzero <- upper>lower
    arrows(x[nonzero],lower[nonzero],x[nonzero],upper[nonzero],angle=90,code=3,length=.035,col=colors[index])
    lines(x,y,type="b",col=colors[index],pch=symbols[index],lwd=1.2)
  }
  if(quantity=="coverage" && scenario=="C1") legend("bottomright",c("K = 4","K = 6"),col=colors,pch=symbols,lty=1,bty="n")
}
mtext(paste0("p = ",dimension,", n = 1000/site, joint TATE; one source shifted when rho > 0"),side=3,outer=TRUE,line=1,cex=1.05)
mtext("50 seeds/cell. Coverage: Wilson intervals; other bars: +/-1.96 MCSE. K also changes total N and source tilts.",side=1,outer=TRUE,line=1,cex=.85)
dev.off()
writeLines(c(report,"Only completed, matched method/baseline runs; fixed K, not a growing-K experiment.",
  "Source 1's control arm remains valid; its mean weight is reported in the accompanying CSV."),
  file.path(output,"provenance.txt"))
