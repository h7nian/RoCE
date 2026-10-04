#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: plot_production_metrics.R CELL_METRICS_CSV OUTPUT_PDF")
metrics <- read.csv(args[1L], stringsAsFactors = FALSE)
methods <- c("one_round_crossfit_ate", "one_round_crossfit_ate_armwise", "target_only_ate")
labels <- c("Direct TATE", "Armwise TATE", "Target only")
colors <- c("#2468B4", "#218559", "#555555")
panels <- data.frame(family = c(rep("negative_transfer", 3L), "shared_shift"),
                     configuration = c("C1", "C2", "C3", "C1"))
grDevices::pdf(args[2L], width = 11, height = 8, onefile = TRUE)
for (metric in c("coverage", "rmse_ratio_target")) {
  par(mfrow = c(3, 4), mar = c(3.1, 3.4, 2.2, .6), oma = c(1, 0, 3.8, 0),
      mgp = c(1.8, .6, 0), cex = .75)
  for (source_count in c(2L, 4L, 8L)) for (panel in seq_len(nrow(panels))) {
    panel_metrics <- metrics[metrics$family == panels$family[panel] &
                        metrics$configuration == panels$configuration[panel] &
                        metrics$K == source_count & metrics$method %in% methods, ]
    limits <- if (metric == "coverage") c(.6, 1) else c(.3, 1.5)
    plot(NA, xlim = c(0, 2.5), ylim = limits, xlab = expression(rho),
         ylab = if (metric == "coverage") "Coverage" else "RMSE / target-only RMSE",
         main = sprintf("%s %s, K=%d; MC n=%s",
                        if (panels$family[panel] == "shared_shift") "Shared" else "NT",
                        panels$configuration[panel], source_count,
                        paste(unique(panel_metrics$n), collapse = "/")))
    abline(h = if (metric == "coverage") .95 else 1, col = "#999999", lty = 2)
    for (index in seq_along(methods)) {
      values <- panel_metrics[panel_metrics$method == methods[index], ]
      values <- values[order(values$rho), ]
      if (metric == "coverage" && index == 1L) {
        arrows(values$rho, pmax(values$coverage_lower, limits[1L]),
               values$rho, values$coverage_upper, angle = 90, code = 3,
               length = .03, col = "#8CAED4")
      }
      lines(values$rho, values[[metric]], col = colors[index], pch = 14 + index,
            type = "b", lwd = 1.2, cex = .65)
    }
    if (source_count == 2L && panel == 1L) {
      legend("bottomleft", labels, col = colors, pch = 15:17, lty = 1,
             bty = "n", cex = .7)
    }
  }
  mtext("TATE at 1000 observations/site, p=100, ten outer folds", outer = TRUE,
        side = 3, line = 2, cex = 1.05)
  mtext("Current v4 raw results; incomplete cells retain their actual counts. Direct-TATE bars: pointwise 95% Wilson intervals.",
        outer = TRUE, side = 3, line = .6, cex = .75)
}
invisible(dev.off())
