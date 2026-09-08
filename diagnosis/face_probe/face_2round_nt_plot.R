#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
files <- list.files("diagnosis/face_probe/validation/chunks2r_nt", pattern="^nt_K4_C1_.*\\.csv$", full.names=TRUE)
if(!length(files)) stop("no NT csvs")
raw <- do.call(rbind, lapply(files, function(f) tryCatch(read.csv(f), error=function(e) NULL)))
keep <- c("one_round_crossfit","two_round_crossfit","target_only")
raw <- raw[raw$method %in% keep, c("method","rho","bias","coverage","ci_width")]
agg <- aggregate(cbind(bias,coverage,ci_width)~method+rho, raw, mean)
agg$rmse <- aggregate(bias~method+rho, raw, function(x) sqrt(mean(x^2)))$bias
agg$n    <- aggregate(bias~method+rho, raw, length)$bias
lab <- c(one_round_crossfit="1-round", two_round_crossfit="2-round", target_only="Target-only")
agg$Method <- factor(lab[agg$method], levels=c("Target-only","1-round","2-round"))
mk <- function(m,v) data.frame(rho=agg$rho, Method=agg$Method, metric=m, value=v)
long <- rbind(mk("Bias",agg$bias),mk("RMSE",agg$rmse),mk("Coverage",agg$coverage),mk("95% CI length",agg$ci_width))
long$metric <- factor(long$metric, levels=c("Bias","RMSE","Coverage","95% CI length"))
href <- data.frame(metric=factor(c("Bias","Coverage"),levels=levels(long$metric)), yint=c(0,0.95))
cm <- c("Target-only"="#1B9E77","1-round"="#E7298A","2-round"="#D95F02")
lt <- c("Target-only"="longdash","1-round"="dashed","2-round"="solid")
p <- ggplot(long, aes(rho,value,colour=Method,linetype=Method,group=Method))+
  geom_hline(data=href,aes(yintercept=yint),linetype="dashed",colour="grey50",linewidth=0.3)+
  geom_line(linewidth=0.9)+facet_wrap(~metric,scales="free_y",ncol=1)+
  scale_colour_manual(values=cm)+scale_linetype_manual(values=lt)+
  guides(colour=guide_legend(nrow=1),linetype=guide_legend(nrow=1))+
  labs(x=expression("source deviation "*rho*" (one non-transportable source, K=4)"),y=NULL)+
  theme_bw(base_size=14)+theme(legend.position="top",panel.grid.minor=element_blank(),
    strip.background=element_rect(fill="grey92"))
ggsave("diagnosis/face_probe/validation/fig_2round_nt.pdf",p,width=6,height=9)
ggsave("docs/figures/sim_2round_nt.pdf",p,width=6,height=9)
cat("wrote fig_2round_nt.pdf\n")
print(agg[order(agg$method,agg$rho),c("method","rho","n","bias","rmse","coverage","ci_width")],row.names=FALSE)
