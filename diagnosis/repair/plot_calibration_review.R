#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
if(!length(arguments) %in% 1:2) stop("Usage: completed_calibration_overview [new_figure_directory_name]")
root <- normalizePath(arguments[1L],mustWork=TRUE)
if(!startsWith(root,"/scratch.global/zhan9381/FACE-HD/")) stop("Figures belong in FACE-HD scratch")
values <- read.csv(file.path(root,"joint_tate_summary.csv"),stringsAsFactors=FALSE)
stopifnot(nrow(values)==9L,all(is.finite(values$mse_difference)),
          all(values$mse_difference_mcse>=0))
figure_name <- if(length(arguments)==2L) arguments[2L] else "figures"
if(!grepl("^[A-Za-z0-9_]+$",figure_name)) stop("Use a simple figure directory name")
figures <- file.path(root,figure_name)
if(dir.exists(figures)) stop("Do not overwrite previous figures")
dir.create(figures)
labels <- c(main_source="Source calibration, main DGP",
            main_target="Target calibration, main DGP",
            shift2_source="Source calibration, shift x2")
colors <- c(main_source="#2463A5",main_target="#A65514",shift2_source="#5A8B3E")
positions <- 9:1
row_labels <- paste(unname(labels[values$study]),values$config,sep=" / ")
plot_mse <- function() {
  par(mar=c(4.3,18,3.5,1),mgp=c(2.5,.7,0),las=1)
  differences <- 1e5*values$mse_difference
  half_width <- 1.96e5*values$mse_difference_mcse
  plot(differences,positions,type="n",ylim=c(.5,9.5),
       xlim=range(c(differences-half_width,differences+half_width,0)),yaxt="n",
       xlab="Paired MSE difference x 100,000 (full minus ablation)",ylab="",
       main="Calibration: small TATE MSE differences at p=100")
  abline(v=0,lty=2,col="gray50")
  abline(h=c(3.5,6.5),col="gray90")
  segments(differences-half_width,positions,differences+half_width,positions,
           col=colors[values$study],lwd=2)
  points(differences,positions,pch=19,col=colors[values$study])
  axis(2,at=positions,labels=row_labels,tick=FALSE,cex.axis=.8)
  mtext("Bars: +/- 1.96 paired Monte Carlo SE; 200 seeds per scenario. Negative favors full calibration.",
        side=1,line=3.1,cex=.72)
}
plot_variance <- function() {
  par(mar=c(4.3,18,3.5,1),mgp=c(2.5,.7,0),las=1)
  plot(values$full_variance_ratio,positions,type="n",ylim=c(.5,9.5),
       xlim=range(c(values$full_variance_ratio,values$ablated_variance_ratio,1)),yaxt="n",
       xlab="Mean reported variance / empirical variance",ylab="",
       main="Calibration: observed TATE variance accuracy")
  abline(v=1,lty=2,col="gray50")
  abline(h=c(3.5,6.5),col="gray90")
  segments(values$ablated_variance_ratio,positions,values$full_variance_ratio,positions,
           col=colors[values$study],lwd=2)
  points(values$ablated_variance_ratio,positions,pch=1,col=colors[values$study])
  points(values$full_variance_ratio,positions,pch=19,col=colors[values$study])
  axis(2,at=positions,labels=row_labels,tick=FALSE,cex.axis=.8)
  legend("bottomright",legend=c("Full calibration","Ablation"),pch=c(19,1),bty="n",cex=.75)
  mtext("Descriptive ratios from 200 paired seeds; paired jackknife uncertainty is in the report.",
        side=1,line=3.1,cex=.72)
}
pdf(file.path(figures,"calibration_p100_mc200.pdf"),width=10,height=5.6)
plot_mse();plot_variance();dev.off()
cat("Saved",file.path(figures,"calibration_p100_mc200.pdf"),"\n")
