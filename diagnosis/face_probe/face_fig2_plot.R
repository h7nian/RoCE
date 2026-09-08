#!/usr/bin/env Rscript
# FACE.pdf Figure-2 analog: RoCE adaptive ensemble weights eta-hat per source
# vs the deviation gradient rho. The deviated source (s1) shrinks toward 0 as it
# becomes non-informative; informative sources retain weight. Reads
# fig2_weights_K*.csv from face_fig2_weights.R.
suppressPackageStartupMessages(library(ggplot2))

val_dir <- "diagnosis/face_probe/validation"
files <- list.files(val_dir, pattern = "^fig2_weights_K[0-9]+\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("no fig2_weights_K*.csv files found yet")
dat <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
dat$Kf <- factor(sprintf("K = %d", dat$K), levels = sprintf("K = %d", sort(unique(dat$K))))

p <- ggplot(dat, aes(rho, eta_mean, colour = source, group = source)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.3) +
  geom_errorbar(aes(ymin = eta_mean - eta_se, ymax = eta_mean + eta_se),
                width = 0.05, linewidth = 0.3) +
  facet_wrap(~ Kf, nrow = 1) +
  labs(x = expression("source deviation " * rho * " (log-odds)"),
       y = expression("adaptive ensemble weight " * hat(eta)[j]),
       colour = "source",
       title = "RoCE adaptive weights shrink the non-informative source (s1)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92"))

out_pdf <- file.path(val_dir, "fig2_weights.pdf")
ggsave(out_pdf, p, width = 3 + 2.6 * length(unique(dat$K)), height = 4)
cat(sprintf("wrote %s (K = {%s})\n", out_pdf, paste(sort(unique(dat$K)), collapse = ",")))
