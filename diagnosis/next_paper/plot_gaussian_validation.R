#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: metrics_csv output_directory")
metrics <- read.csv(arguments[[1L]], stringsAsFactors = FALSE)
output <- arguments[[2L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Plots must be on FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(device = function(...) stop("An explicit graphics device on FACE-HD scratch is required"))
methods <- c("target_only", "oracle_gls", "naive_screen_gls", "search_marginal_bound",
             "search_known_correlation", "search_wrong_independence")
labels <- c("Target only", "Oracle GLS", "Naive screen + GLS", "Marginal reference",
            "Known-correlation search", "Wrong independence")
colors <- c("#4D4D4D", "#117733", "#CC6677", "#DDCC77", "#0072B2", "#AA4499")
symbols <- c(1, 2, 4, 5, 16, 17)
selected <- metrics[metrics$assumption_valid & metrics$bias_scale == 1 &
                    metrics$bias_pattern == "positive", ]
for (metric in c("coverage", "length_ratio")) {
  pdf_path <- file.path(output, paste0(metric, ".pdf"))
  pdf(pdf_path, width = 11, height = 4.4, useDingbats = FALSE)
  par(mfrow = c(1, 3), mar = c(4.2, 4, 2.6, 1), oma = c(3.8, 0, 1.5, 0))
  for (correlation in c(0, .5, .9)) {
    cell <- selected[selected$shared_correlation == correlation, ]
    plot(NA, xlim = c(4, 64), ylim = if (metric == "coverage") c(.3, 1) else c(0, 1.25),
         log = "x", xaxt = "n", xlab = "Number of sources K",
         ylab = if (metric == "coverage") "Coverage" else "Mean length / target-only length",
         main = paste("Shared correlation =", correlation))
    axis(1, at = c(4, 8, 16, 32, 64), labels = c(4, 8, 16, 32, 64))
    abline(h = if (metric == "coverage") .95 else 1, col = "gray60", lty = 2)
    for (index in seq_along(methods)) {
      rows <- cell[cell$method == methods[index], ]
      rows <- rows[order(rows$num_sources), ]
      if (metric == "coverage") {
        values <- rows$coverage
      } else {
        reference <- cell[cell$method == "target_only", ]
        values <- rows$mean_length / reference$mean_length[match(rows$num_sources, reference$num_sources)]
      }
      lines(rows$num_sources, values, type = "b", col = colors[index], pch = symbols[index], lwd = 1.5)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("bottom", inset = .01, legend = labels, col = colors, pch = symbols, lty = 1,
         ncol = 3, bty = "n", cex = .85)
  mtext("Gaussian reference experiment; local bias scale = 1; 1,000 repeats per setting", side = 3, line = -1.2, cex = .85)
  dev.off()
  converter <- Sys.which("gs")
  if (nzchar(converter)) {
    preview <- file.path(output, paste0(metric, ".png"))
    # The R module's shared-library path can conflict with system Ghostscript.
    status <- system2(converter, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150",
                                 paste0("-sOutputFile=", shQuote(preview)), shQuote(pdf_path)),
                      env = "LD_LIBRARY_PATH=")
    if (status != 0L || !file.exists(preview)) stop("PDF preview rendering failed")
  }
}
