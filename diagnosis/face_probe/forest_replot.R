#!/usr/bin/env Rscript
# Re-plot the RHC forest from the cached method estimates (numbers unchanged):
# unify to the 6 simulation methods, publication labels, no title, standardized
# x-axis. Package-free (ggplot2 only) so it needs no RoCE reinstall.
suppressPackageStartupMessages(library(ggplot2))
csv <- "results/real_data/rhc_K5_death30_kf10_A1_ninsclas_Mt5_methods.csv"
df <- read.csv(csv, stringsAsFactors = FALSE)
# The re-run's methods.csv already carries the publication labels.
keep <- c("Target-only","SS","IVW","Federated-DR","Pooled-DR","RoCE")
df <- df[df$method %in% keep, ]
df$method <- factor(df$method, levels = rev(keep))
col_map <- c("RoCE"="#D95F02","Target-only"="#1B9E77","SS"="#7570B3",
             "IVW"="#E7298A","Federated-DR"="#66A61E","Pooled-DR"="#E6AB02")
p <- ggplot(df, aes(y = method, x = estimate, colour = method)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25, linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = col_map, guide = "none") +
  labs(x = "Potential-outcome mean", y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor = element_blank())
out <- c("results/real_data/rhc_K5_death30_kf10_A1_ninsclas_Mt5_forest.pdf",
         "docs/figures/rhc_forest.pdf")
for (o in out) ggsave(o, p, width = 6.5, height = 3.6)
cat("wrote forest:", paste(out, collapse=", "), "| methods:", paste(rev(levels(df$method)), collapse=", "), "\n")
