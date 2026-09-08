#!/usr/bin/env Rscript
# Aggregate the 1-round vs 2-round EM sweep and plot bias/RMSE/coverage/CI-length
# vs effect-modification strength, for mu1 (the treated-arm mean the EM acts on).
suppressPackageStartupMessages(library(ggplot2))
dir <- "diagnosis/face_probe/validation/chunks2r"
files <- list.files(dir, pattern="^em_K4_C1_.*\\.csv$", full.names=TRUE)
if (!length(files)) stop("no EM chunk CSVs yet")
raw <- do.call(rbind, lapply(files, function(f) tryCatch(read.csv(f), error=function(e) NULL)))
keep <- c("one_round_crossfit","two_round_crossfit","target_only")
raw <- raw[raw$method %in% keep, c("method","em","p","bias","coverage","ci_width")]
agg <- aggregate(cbind(bias,coverage,ci_width) ~ method+em+p, raw, mean)
agg$rmse <- aggregate(bias ~ method+em+p, raw, function(x) sqrt(mean(x^2)))$bias
agg$n    <- aggregate(bias ~ method+em+p, raw, length)$bias
lab <- c(one_round_crossfit="1-round", two_round_crossfit="2-round", target_only="Target-only")
agg$Method <- factor(lab[agg$method], levels=c("Target-only","1-round","2-round"))
agg$Pf <- factor(paste0("p = ", agg$p), levels=c("p = 10","p = 50"))
mk <- function(metric, value) data.frame(em=agg$em, Pf=agg$Pf, Method=agg$Method, metric=metric, value=value)
long <- rbind(mk("Bias",agg$bias), mk("RMSE",agg$rmse),
              mk("Coverage",agg$coverage), mk("95% CI length",agg$ci_width))
long$metric <- factor(long$metric, levels=c("Bias","RMSE","Coverage","95% CI length"))
href <- data.frame(metric=factor(c("Bias","Coverage"),levels=levels(long$metric)), yint=c(0,0.95))
col_map <- c("Target-only"="#1B9E77","1-round"="#E7298A","2-round"="#D95F02")
lt_map  <- c("Target-only"="longdash","1-round"="dashed","2-round"="solid")
p <- ggplot(long, aes(em, value, colour=Method, linetype=Method, group=Method)) +
  geom_hline(data=href, aes(yintercept=yint), linetype="dashed", colour="grey50", linewidth=0.3) +
  geom_line(linewidth=0.9) +
  facet_grid(metric ~ Pf, scales="free_y") +
  scale_colour_manual(values=col_map) + scale_linetype_manual(values=lt_map) +
  guides(colour=guide_legend(nrow=1), linetype=guide_legend(nrow=1)) +
  labs(x="source effect-modification strength  s", y=NULL) +
  theme_bw(base_size=15) +
  theme(legend.position="top", legend.text=element_text(size=14),
        axis.title=element_text(size=15), strip.text=element_text(size=14),
        panel.grid.minor=element_blank(), strip.background=element_rect(fill="grey92"))
ggsave("diagnosis/face_probe/validation/fig_2round_em.pdf", p, width=8, height=8.8)
ggsave("docs/figures/sim_2round_em.pdf", p, width=8, height=8.8)
cat(sprintf("wrote fig_2round_em.pdf (%d chunk files)\n", length(files)))
print(agg[order(agg$p, agg$method, agg$em), c("p","method","em","n","bias","rmse","coverage","ci_width")], row.names=FALSE)
