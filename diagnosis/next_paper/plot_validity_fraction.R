#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L || !startsWith(arguments[[2L]], "/scratch.global/zhan9381/FACE-HD/")) {
  stop("Usage: metrics_csv scratch_output")
}
metrics <- read.csv(arguments[[1L]], stringsAsFactors = FALSE)
output <- arguments[[2L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(device = function(...) stop("An explicit scratch graphics device is required"))
selected <- subset(metrics, shared_correlation == 0 & bias_scale == 2)
series <- data.frame(method = c("target_only", "oracle_gls", "naive_screen_gls",
  "search_known_correlation", "search_known_correlation"),
  vote_rule = c("partial", "partial", "partial", "partial", "all_valid"),
  label = c("Target only", "Oracle GLS", "Naive screen + GLS", "Search: partial votes", "Search: all q votes"),
  color = c("#555555", "#228833", "#CC6677", "#0077BB", "#AA3377"),
  symbol = c(1L, 2L, 4L, 16L, 17L), stringsAsFactors = FALSE)
pdf_path <- file.path(output, "validity_fraction.pdf")
pdf(pdf_path, width = 11, height = 7, useDingbats = FALSE)
par(mfrow = c(2, 3), mar = c(4, 4.2, 2.6, 1), oma = c(3.5, 0, 1.5, 0))
for (metric in c("coverage", "length_ratio_to_oracle")) {
  for (fraction in c(.25, .5, .75)) {
    cell <- selected[selected$unshifted_fraction == fraction, ]
    ylim <- if (metric == "coverage") c(0, 1) else c(0, max(selected$length_ratio_to_oracle[selected$method %in% series$method]))
    plot(NA, xlim = c(4, 64), ylim = ylim, log = "x", xaxt = "n",
      xlab = "Number of sources K", ylab = if (metric == "coverage") "Coverage" else "Mean length / oracle length",
      main = paste("Guaranteed valid fraction:", fraction))
    axis(1, at = c(4, 8, 16, 32, 64), labels = c(4, 8, 16, 32, 64))
    abline(h = if (metric == "coverage") .95 else 1, lty = 2, col = "gray60")
    for (index in seq_len(nrow(series))) {
      rows <- cell[cell$method == series$method[index] & cell$vote_rule == series$vote_rule[index], ]
      rows <- rows[order(rows$num_sources), ]
      lines(rows$num_sources, rows[[metric]], type = "b", col = series$color[index],
            pch = series$symbol[index], lwd = 1.5)
    }
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", inset = .01, legend = series$label, col = series$color,
       pch = series$symbol, lty = 1, ncol = 3, bty = "n", cex = .9)
mtext("Known-variance Gaussian reference; local bias scale = 2; 1,000 repeats per setting", side = 3, line = -1.2, cex = .9)
dev.off()
converter <- Sys.which("gs")
if (nzchar(converter)) {
  preview <- file.path(output, "validity_fraction.png")
  status <- system2(converter, c("-q", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150",
    paste0("-sOutputFile=", shQuote(preview)), shQuote(pdf_path)), env = "LD_LIBRARY_PATH=")
  if (status != 0L || !file.exists(preview)) stop("PDF preview rendering failed")
}
