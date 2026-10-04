#!/usr/bin/env Rscript
# Standalone diagnostic figure; plotted points are training draws, not CI limits.
arguments <- commandArgs(trailingOnly=TRUE)
if(length(arguments)!=2L) stop("Usage: complete_review_directory new_figure_directory")
review <- normalizePath(arguments[1L],mustWork=TRUE)
output <- arguments[2L]
if(!startsWith(output,"/scratch.global/zhan9381/FACE-HD/") || dir.exists(output)) stop("Use new FACE-HD scratch output")
bias <- read.csv(file.path(review,"bias_decomposition.csv"),stringsAsFactors=FALSE)
variance <- read.csv(file.path(review,"oracle_valid_aggregate.csv"),stringsAsFactors=FALSE)
variance <- variance[variance$covariance_type=="estimated",]
colors <- c(translated_source="#3465A4",common_outcome="#C46B35")
source_counts <- c(2L,8L,16L)
for(data in list(bias,variance)) {
  counts <- table(data$scenario,data$rho,data$K,data$construction)
  if(!all(counts==5L) || !setequal(data$K,source_counts)) stop("Figure requires all five prescribed training draws per cell")
}
dir.create(output,recursive=TRUE)
pdf(file.path(output,"fitted_population_diagnostics.pdf"),width=10.5,height=7.5,onefile=TRUE)
for(quantity in c("bias","variance")) {
  layout(matrix(c(1:6,7,7,7),nrow=3L,byrow=TRUE),heights=c(1,1,.16))
  par(mar=c(3.2,4.2,2.4,.7),oma=c(2,0,3,0),mgp=c(2.2,.65,0),
    cex.lab=.85,cex.axis=.8,cex.main=.95)
  data <- if(quantity=="bias") bias else variance
  values <- if(quantity=="bias") abs(data$standardized_tate_shift) else data$reported_to_fixed_weight_variance
  limits <- if(quantity=="bias") c(0,max(.5,values)*1.12) else range(c(.9,1.1,values))*c(.95,1.05)
  for(rho in c(0,.5)) for(scenario in c("C1","C2","C3")) {
    plot(NA,xlim=c(.6,3.4),ylim=limits,xaxt="n",xlab="Sources (K)",
      ylab=if(quantity=="bias") "|TATE drift| / oracle SD" else "Reported / actual variance",
      main=paste(scenario,"  rho =",rho))
    axis(1,at=1:3,labels=source_counts)
    abline(h=if(quantity=="bias") 0 else 1,col="grey70",lty=2)
    for(index in seq_along(colors)) {
      construction <- names(colors)[index]
      medians <- numeric(3L)
      for(k in seq_along(source_counts)) {
        selected <- which(data$scenario==scenario & data$rho==rho & data$K==source_counts[k] & data$construction==construction)
        selected <- selected[order(data$seed[selected])]
        positions <- k+c(-.12,.12)[index]+seq(-.035,.035,length.out=5L)
        points(positions,values[selected],pch=16,col=adjustcolor(colors[index],alpha.f=.7),cex=.7)
        medians[k] <- median(values[selected])
      }
      lines(1:3+c(-.12,.12)[index],medians,col=colors[index],lwd=1.5)
    }
  }
  mtext(if(quantity=="bias") "Fitted valid-source bias: oracle covariance and validity labels" else
    "Covariance estimation: realized GLS weights held fixed",outer=TRUE,line=1,font=2)
  mtext("n = 1000/site; p = 4; 500 training / 500 evaluation; points: five training draws; lines: medians",outer=TRUE,side=1,line=.5,cex=.8)
  par(mar=rep(0,4))
  plot.new()
  legend("center",legend=c("Translated source OR","Common target IPW OR"),col=colors,lty=1,pch=16,horiz=TRUE,bty="n",cex=.9)
}
dev.off()
