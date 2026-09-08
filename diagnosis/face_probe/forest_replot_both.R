#!/usr/bin/env Rscript
# Two-panel RHC forest: treated-arm mean mu^1 (A=1) and control-arm mean mu^0
# (A=0), the six unified methods, publication labels, simulation palette. Shared
# x-axis so the treated cluster sits visibly to the right of the control cluster
# (the positive ATE). Package-free (ggplot2 only) so it needs no RoCE reinstall.
suppressPackageStartupMessages(library(ggplot2))
read_arm <- function(csv, arm) {
  df <- read.csv(csv, stringsAsFactors = FALSE)
  df$arm <- arm
  df
}
a1 <- "results/real_data/rhc_K5_death30_kf10_A1_ninsclas_Mt5_methods.csv"
a0 <- "results/real_data/rhc_K5_death30_kf10_A0_ninsclas_Mt5_methods.csv"
df <- rbind(read_arm(a1, "Treated arm (RHC)"),
            read_arm(a0, "Control arm (no RHC)"))
keep <- c("Target-only", "SS", "IVW", "Federated-DR", "Pooled-DR", "RoCE")
df <- df[df$method %in% keep, ]
df$method <- factor(df$method, levels = rev(keep))
df$arm <- factor(df$arm, levels = c("Treated arm (RHC)", "Control arm (no RHC)"))
col_map <- c("RoCE"="#D95F02","Target-only"="#1B9E77","SS"="#7570B3",
             "IVW"="#E7298A","Federated-DR"="#66A61E","Pooled-DR"="#E6AB02")
p <- ggplot(df, aes(y = method, x = estimate, colour = method)) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25, linewidth = 1.0) +
  geom_point(size = 3.2) +
  facet_wrap(~ arm, nrow = 1) +
  scale_colour_manual(values = col_map, guide = "none") +
  labs(x = "Potential-outcome mean", y = NULL) +
  theme_bw(base_size = 18) +
  theme(panel.grid.major.y = element_line(colour = "grey92"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA),
        axis.title = element_text(size = 19),
        axis.text = element_text(size = 17),
        strip.text = element_text(face = "bold", size = 18))
out <- c("results/real_data/rhc_forest_both.pdf", "docs/figures/rhc_forest_both.pdf")
for (o in out) ggsave(o, p, width = 8.2, height = 3.4)
cat("wrote two-panel forest:", paste(out, collapse = ", "), "\n")
