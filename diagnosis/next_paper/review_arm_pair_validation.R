#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: completed_gaussian_run new_scratch_output")
input <- arguments[[1L]]
output <- arguments[[2L]]
if (any(!startsWith(arguments, "/scratch.global/zhan9381/FACE-HD/"))) stop("Use FACE-HD scratch")
if (!file.exists(file.path(input, "COMPLETE"))) stop("The prescribed study must be complete")
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
configuration <- readRDS(file.path(input, "configuration.rds"))
truth <- configuration$mu1 - configuration$mu0
metrics <- read.csv(file.path(input, "metrics.csv"), stringsAsFactors = FALSE)
records <- list()
maximum_identity_error <- 0
target_by_design <- list()
for (cell in seq_len(nrow(configuration$settings))) {
  saved <- readRDS(file.path(input, "cells", paste0(cell, ".rds")))
  stopifnot(saved$cell == cell, identical(saved$configuration, configuration))
  rows <- saved$records
  target <- rows[rows$method == "target", ]
  target <- target[order(target$iteration), ]
  stopifnot(nrow(target) == configuration$repeats)
  setting <- configuration$settings[cell, ]
  key <- paste(setting$num_sources, setting$profile)
  endpoints <- as.matrix(target[, c("lower", "upper")])
  if (is.null(target_by_design[[key]])) target_by_design[[key]] <- endpoints else
    stopifnot(max(abs(target_by_design[[key]] - endpoints)) < 1e-12)
  target_length <- target$upper - target$lower
  if (setting$local_bias == 0) {
    oracle <- rows[rows$method == "oracle_arm_gls", c("lower", "upper")]
    pooled <- rows[rows$method == "pooled_arm_gls", c("lower", "upper")]
    maximum_identity_error <- max(maximum_identity_error, abs(as.matrix(oracle) - as.matrix(pooled)))
  }
  for (method in unique(rows$method)) {
    selected <- rows[rows$method == method, ]
    selected <- selected[order(selected$iteration), ]
    stopifnot(identical(selected$iteration, target$iteration))
    covered <- selected$lower <= truth & selected$upper >= truth
    reference <- metrics[metrics$cell == cell & metrics$method == method, ]
    stopifnot(nrow(reference) == 1L, abs(mean(covered) - reference$coverage) < 1e-12)
    length_ratio <- (selected$upper - selected$lower) / target_length
    record <- cbind(reference, data.frame(mean_length_ratio = mean(length_ratio),
      length_ratio_mcse = sd(length_ratio) / sqrt(length(length_ratio)),
      undercoverage_p = pbinom(sum(covered), length(covered), .95)))
    records[[length(records) + 1L]] <- record
    if (grepl("disjoint", setting$profile) && method == "scalar_tate") {
      stopifnot(max(abs(as.matrix(selected[, c("lower", "upper")]) - endpoints)) < 1e-12)
    }
  }
}
stopifnot(maximum_identity_error < 1e-12)
result <- do.call(rbind, records)
reference_methods <- result$method != "pooled_arm_gls"
result$undercoverage_p_holm <- NA_real_
result$undercoverage_p_holm[reference_methods] <- p.adjust(result$undercoverage_p[reference_methods], "holm")
write.csv(result, file.path(output, "metrics.csv"), row.names = FALSE)

methods <- c("pair_half", "pair_all", "scalar_tate", "target", "oracle_arm_gls", "pooled_arm_gls")
colors <- c("#0072B2", "#56B4E9", "#009E73", "#333333", "#CC79A7", "#D55E00")
draw <- function(field, label, limits) {
  layout(matrix(c(1:6, 7, 7, 7), 3, 3, byrow = TRUE), heights = c(1, 1, .18))
  par(mar = c(3.5, 4.2, 2, .9), oma = c(0, 0, 1, 0), mgp = c(1.9, .55, 0),
      cex.axis = .8, cex.lab = .8)
  for (count in c(4L, 16L)) for (profile in c("treated_half", "disjoint_half", "disjoint_quarter")) {
    selected <- result[result$num_sources == count & result$profile == profile, ]
    plot(NA, xlim = c(0, 4), ylim = limits, xlab = "Local bias multiplier", ylab = label,
         main = paste("K =", count, ":", gsub("_", " ", profile)), cex.main = .85)
    abline(h = if (field == "coverage") .95 else 1, lty = 3, col = "gray50")
    for (index in seq_along(methods)) {
      values <- selected[selected$method == methods[index], ]
      values <- values[order(values$local_bias), ]
      lines(values$local_bias, values[[field]], col = colors[index], lty = index, lwd = 1.5, type = "b", pch = index)
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = gsub("_", " ", methods), col = colors, lty = seq_along(methods),
         pch = seq_along(methods), ncol = 3, bty = "n", cex = .75)
}
pdf(file.path(output, "arm_pair_reference.pdf"), width = 10, height = 7, onefile = TRUE)
draw("coverage", "Coverage", c(0, 1))
draw("mean_length_ratio", "Mean length / target length", range(c(1, result$mean_length_ratio)))
dev.off()
ghostscript <- Sys.which("gs")
if (nzchar(ghostscript)) {
  # R's MSI module supplies an OpenSSL library incompatible with system gs.
  # Remove that search path only in the renderer's child environment.
  status <- system2("env", c("-u", "LD_LIBRARY_PATH", shQuote(ghostscript),
    "-dSAFER", "-dBATCH", "-dNOPAUSE", "-sDEVICE=png16m", "-r150",
    "-dFirstPage=2", "-dLastPage=2", paste0("-sOutputFile=", shQuote(file.path(output, "interval_length.png"))),
    shQuote(file.path(output, "arm_pair_reference.pdf"))))
  if (status != 0L) stop("PDF rendering failed; the original PDF and review tables are retained")
}
writeLines(c("# Completed arm-pair Gaussian reference", "",
  "24 settings, 1000 draws each, known full joint arm covariance. The six methods share each draw.",
  "Mean interval lengths are compared with the ordinary 95% target interval.",
  "The pair references cover conservatively and are usually longer than target-only inference.",
  "They are useful validity references, not a solution to the efficiency problem.",
  "No fitted high-dimensional coverage claim follows from this experiment.", "",
  paste("Maximum all-valid oracle/pooled endpoint identity error:", maximum_identity_error),
  paste("Holm-adjusted undercoverage flags among reference methods:",
    sum(result$undercoverage_p_holm[reference_methods] < .05))), file.path(output, "README.md"))
cat("ARM_PAIR_REFERENCE_REVIEW_PASSED\n")
