#!/usr/bin/env Rscript
# Export a completed paired variance panel without refitting.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)<2L || length(arguments)>3L) stop("Usage: layer_report_directory new_scratch_figure_directory [p]")
report <- normalizePath(arguments[[1L]],mustWork=TRUE)
output <- arguments[[2L]]
dimension <- if(length(arguments)==3L) suppressWarnings(as.numeric(arguments[[3L]])) else 10L
if(length(dimension)!=1L || !is.finite(dimension) || dimension<1 || dimension!=floor(dimension)) stop("p must be a positive integer")
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use a new FACE-HD scratch directory")
metrics <- read.csv(file.path(report,"paired_metrics.csv"),stringsAsFactors=FALSE)
components <- read.csv(file.path(report,"paired_variance_components.csv"),stringsAsFactors=FALSE)
metrics <- subset(metrics,p==dimension & aggregation_mode=="joint_tate" & quantity=="tate")
components <- subset(components,p==dimension & aggregation_mode=="joint_tate" & component=="tate" & inference=="eta_sensitivity")
stopifnot(nrow(metrics)==18L,nrow(components)==18L,
  all(metrics$n_success==50L),all(components$n_pairs==50L),
  !anyDuplicated(metrics[c("config","K","crossfit_layers")]),
  !anyDuplicated(components[c("config","K","crossfit_layers")]))
dir.create(output,recursive=TRUE)
colors <- c("#0072B2","#D55E00")
symbols <- c(16,17)
positions <- c(-.035,.035)
scenarios <- c("C1","C2","C3")
counts <- c(2,4,6)
ratio_limits <- extendrange(c(1,metrics$variance_ratio),f=.12)
pdf(file.path(output,"layer_variance_comparison.pdf"),width=10,height=4.8,
  title="Paired two/three-layer TATE variance diagnosis",useDingbats=FALSE)
par(mfrow=c(1,3),mar=c(4,4,2.4,1),oma=c(3,0,3,0),mgp=c(2.5,.7,0))
for(scenario in scenarios) {
  plot(NA,xlim=c(.8,3.2),ylim=ratio_limits,xaxt="n",xlab="Number of sources K",
    ylab="Mean reported variance / MC variance",main=scenario)
  axis(1,at=1:3,labels=counts)
  abline(h=1,lty=2,col="gray45")
  for(index in 1:2) {
    rows <- subset(metrics,config==scenario & crossfit_layers==index+1L)
    rows <- rows[match(counts,rows$K),]
    lines((1:3)+positions[index],rows$variance_ratio,type="b",pch=symbols[index],col=colors[index],lwd=1.4)
  }
  if(scenario=="C1") legend("topright",c("Two layers","Three layers"),col=colors,pch=symbols,lty=1,bty="n",cex=.9)
}
mtext(sprintf("TATE variance: n = 1000/site, p = %d, joint aggregation",dimension),side=3,outer=TRUE,line=1,cex=1.15)
mtext("50 matched seeds/cell; identical data and outer fits. Descriptive ratios; uncertainty on page 2.",side=1,outer=TRUE,line=1,cex=.9)

scale <- 1e4
limit <- max(abs(components$reported_minus_empirical)+1.96*components$discrepancy_jackknife_se)*scale*1.08
for(scenario in scenarios) {
  plot(NA,xlim=c(.8,3.2),ylim=c(-limit,limit),xaxt="n",xlab="Number of sources K",
    ylab=expression(10^4~"x (reported variance - MC variance)"),main=scenario)
  axis(1,at=1:3,labels=counts)
  abline(h=0,lty=2,col="gray45")
  for(index in 1:2) {
    rows <- subset(components,config==scenario & crossfit_layers==index+1L)
    rows <- rows[match(counts,rows$K),]
    x <- (1:3)+positions[index]
    y <- scale*rows$reported_minus_empirical
    half_width <- scale*1.96*rows$discrepancy_jackknife_se
    arrows(x,y-half_width,x,y+half_width,angle=90,code=3,length=.045,col=colors[index])
    points(x,y,pch=symbols[index],col=colors[index],cex=1.1)
  }
  if(scenario=="C1") legend("topright",c("Two layers","Three layers"),col=colors,pch=symbols,bty="n",cex=.9)
}
mtext("Monte Carlo uncertainty in variance accuracy",side=3,outer=TRUE,line=1,cex=1.15)
mtext("Bars: discrepancy +/- 1.96 delete-one jackknife SE. Exploratory, not multiplicity adjusted.",side=1,outer=TRUE,line=1,cex=.9)
dev.off()
writeLines(c("Source: paired_metrics.csv and paired_variance_components.csv",report,
  paste("Complete p =",dimension,"panel only; no unpaired or incomplete interim results are pooled.")),file.path(output,"provenance.txt"))
