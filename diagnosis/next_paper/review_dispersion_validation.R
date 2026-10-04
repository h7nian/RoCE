#!/usr/bin/env Rscript
# Review new intervals against references exploiting the same validity assumptions.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: previous_gaussian_run dispersion_run new_scratch_output")
if (any(!startsWith(arguments, "/scratch.global/zhan9381/FACE-HD/"))) stop("Use FACE-HD scratch")
previous <- arguments[[1L]]
input <- arguments[[2L]]
output <- arguments[[3L]]
if (!file.exists(file.path(input, "COMPLETE"))) stop("The dispersion study must be complete")
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
configuration <- readRDS(file.path(previous, "configuration.rds"))
metrics <- read.csv(file.path(input, "metrics.csv"), stringsAsFactors = FALSE)
truth <- configuration$mu1 - configuration$mu0
records <- list()
baseline_records <- list()
for (cell in seq_len(nrow(configuration$settings))) {
  reference <- readRDS(file.path(previous, "cells", paste0(cell, ".rds")))
  saved <- readRDS(file.path(input, "cells", paste0(cell, ".rds")))
  stopifnot(identical(reference$configuration, configuration), identical(saved$configuration, configuration))
  setting <- configuration$settings[cell, ]
  count <- setting$num_sources
  width <- count + 1L
  q1 <- if (setting$profile == "disjoint_quarter") count / 4 else count / 2
  q0 <- if (setting$profile == "treated_half") count else count - q1
  covariance <- reference$covariance
  set.seed(68000L + match(count, c(4L, 16L)) * 100L +
    match(setting$profile, c("treated_half", "disjoint_half", "disjoint_quarter")))
  draws <- sweep(matrix(rnorm(configuration$repeats * 2L * width), configuration$repeats, 2L * width) %*%
    chol(covariance), 2L, reference$population_means, "+")
  stopifnot(identical(digest::digest(draws, algo = "sha256"), saved$data_sha256))
  # Only target coordinates are always known valid. An entire source arm is
  # included iff its DECLARED validity lower bound equals the source count.
  included <- c(if (q1 == count) seq_len(width) else 1L,
                if (q0 == count) width + seq_len(width) else width + 1L)
  design <- cbind(c(rep(1, width), rep(0, width)), c(rep(0, width), rep(1, width)))
  precision_design <- solve(covariance[included, included], design[included, ])
  information <- crossprod(design[included, ], precision_design)
  contrast <- c(1, -1)
  weights <- numeric(2L * width)
  weights[included] <- drop(precision_design %*% solve(information, contrast))
  variance <- drop(crossprod(contrast, solve(information, contrast)))
  target_variance <- covariance[1, 1] + covariance[width + 1, width + 1] - 2 * covariance[1, width + 1]
  stopifnot(variance > 0, variance <= target_variance + 1e-12)
  if (q1 < count && q0 < count) stopifnot(abs(variance - target_variance) < 1e-12)
  point <- drop(draws %*% weights)
  lower <- point - qnorm(.975) * sqrt(variance)
  upper <- point + qnorm(.975) * sqrt(variance)
  covered <- lower <= truth & upper >= truth
  baseline_length <- mean(upper - lower)
  baseline_records[[cell]] <- data.frame(cell = cell, setting, method = "guaranteed_arm_gls",
    repeats = length(point), coverage = mean(covered), mean_length = baseline_length,
    mean_length_ratio = sqrt(variance / target_variance))
  for (method in unique(saved$records$method)) {
    rows <- saved$records[saved$records$method == method, ]
    rows <- rows[order(rows$iteration), ]
    stopifnot(identical(rows$iteration, seq_len(configuration$repeats)))
    successes <- sum(rows$lower <= truth & rows$upper >= truth)
    metric <- metrics[metrics$cell == cell & metrics$method == method, ]
    stopifnot(nrow(metric) == 1L, abs(successes / nrow(rows) - metric$coverage) < 1e-12)
    ratios <- (rows$upper - rows$lower) / baseline_length
    records[[length(records) + 1L]] <- cbind(metric,
      data.frame(guaranteed_arm_coverage = mean(covered), length_ratio_vs_guaranteed_arm = mean(ratios),
        length_ratio_vs_guaranteed_arm_mcse = sd(ratios) / sqrt(length(ratios)),
        undercoverage_p = pbinom(successes, nrow(rows), .95)))
  }
}
records <- do.call(rbind, records)
records$undercoverage_p_holm <- p.adjust(records$undercoverage_p, "holm")
write.csv(records, file.path(output, "paired_metrics.csv"), row.names = FALSE)
write.csv(do.call(rbind, baseline_records), file.path(output, "guaranteed_arm_metrics.csv"), row.names = FALSE)
writeLines(c("# Dispersion confidence-interval reference", "",
  "24 settings x1000 paired draws, with known full Gaussian covariance.",
  "The guaranteed-arm baseline exploits every arm explicitly assumed fully valid; other source coordinates are not used.",
  "It does not inspect the true validity labels inside a partially valid arm.",
  "Do not attribute gains from a known-valid control arm to the new bias bound.",
  "Compare exact and conservative factors, anchored and unanchored intervals, and retain all results.",
  paste("Holm undercoverage flags among the dispersion references:", sum(records$undercoverage_p_holm < .05)),
  "Fitted high-dimensional validity, general efficiency and covariance estimation remain open."),
  file.path(output, "README.md"))
cat("DISPERSION_PAIRED_REVIEW_PASSED\n")
