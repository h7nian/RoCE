#!/usr/bin/env Rscript
# Compare completed Gaussian variance-validation cells on their saved draws.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: source_file input_directory output_directory")
source_file <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(source_file)
input <- arguments[[2L]]
output <- arguments[[3L]]
if (any(!startsWith(c(input, output), "/scratch.global/zhan9381/FACE-HD/"))) stop("Use FACE-HD scratch")
configuration <- readRDS(file.path(input, "configuration.rds"))
expected_hash <- configuration$code_md5[basename(names(configuration$code_md5)) == "source_confidence_sets.R"]
stopifnot(length(expected_hash) == 1L, unname(tools::md5sum(source_file)) == unname(expected_hash))
checkpoints <- list.files(file.path(input, "cells"), pattern = "^[0-9]+[.]rds$", full.names = TRUE)
if (!length(checkpoints)) stop("No completed validation cells")
truth <- configuration$truth
alpha <- .05
records <- list()

for (path in checkpoints) {
  saved <- readRDS(path)
  stopifnot(identical(saved$configuration, configuration),
            length(saved$draws) == configuration$repeats)
  setting <- configuration$settings[saved$cell, ]
  count <- setting$num_sources
  sample_size <- setting$n_per_site
  shared <- .5 / sample_size
  anchor_private <- 1 / sample_size
  private <- exp(seq(log(.25), log(4), length.out = count)) / sample_size
  valid_minimum <- ceiling(.75 * count)
  valid <- seq.int(count - valid_minimum + 1L, count)
  valid_indices <- c(1L, valid + 1L)
  precision <- 1 / c(anchor_private, private)
  oracle_weights <- precision[valid_indices] / sum(precision[valid_indices])
  oracle_se <- sqrt(shared + 1 / sum(precision[valid_indices]))
  calibration <- make_source_calibration(count, valid_minimum, floor(count / 2) + 1L,
    method = "shared_target_gaussian", shared_variance = shared, source_private_variances = private)
  critical <- qnorm(1 - alpha / 2)
  for (iteration in seq_along(saved$draws)) {
    draw <- saved$draws[[iteration]]
    estimates <- draw$estimates
    original <- saved$records[[iteration]]
    stopifnot(abs(original$lower - draw$interval$lower) < 1e-14,
              abs(original$upper - draw$interval$upper) < 1e-14)
    known <- source_search_interval(estimates[1L], sqrt(shared + anchor_private),
      estimates[-1L], sqrt(shared + private), calibration)
    oracle <- sum(oracle_weights * estimates[valid_indices])
    included <- c(TRUE, abs(estimates[-1L] - estimates[1L]) <= 2 * sqrt(anchor_private + private))
    naive_weights <- precision[included] / sum(precision[included])
    naive <- sum(naive_weights * estimates[included])
    naive_se <- sqrt(shared + 1 / sum(precision[included]))
    target_bound_radius <- qnorm(1 - (alpha - configuration$variance_alpha) / 2) *
      sqrt(draw$calibration$anchor_variance_upper)
    limits <- rbind(
      source_variance_region = c(draw$interval$lower, draw$interval$upper),
      source_known_variance = c(known$lower, known$upper),
      target_known_variance = estimates[1L] + c(-1, 1) * critical * sqrt(shared + anchor_private),
      target_variance_bound = estimates[1L] + c(-1, 1) * target_bound_radius,
      oracle_gls = oracle + c(-1, 1) * critical * oracle_se,
      naive_screen_gls = naive + c(-1, 1) * critical * naive_se)
    stopifnot(all(is.finite(limits)), all(limits[, 1] <= limits[, 2]))
    records[[length(records) + 1L]] <- data.frame(cell = saved$cell, num_sources = count,
      n_per_site = sample_size, iteration = iteration, method = rownames(limits),
      lower = limits[, 1L], upper = limits[, 2L], truth = truth, row.names = NULL)
  }
}
records <- do.call(rbind, records)
metrics <- do.call(rbind, lapply(split(records, paste(records$cell, records$method)), function(rows) {
  count <- nrow(rows)
  coverage <- mean(rows$lower <= truth & truth <= rows$upper)
  z <- qnorm(.975)
  denominator <- 1 + z^2 / count
  wilson_center <- (coverage + z^2 / (2 * count)) / denominator
  wilson_radius <- z * sqrt((coverage * (1 - coverage) + z^2 / (4 * count)) / count) / denominator
  center_error <- (rows$lower + rows$upper) / 2 - truth
  data.frame(cell = rows$cell[1L], num_sources = rows$num_sources[1L], n_per_site = rows$n_per_site[1L],
    method = rows$method[1L], repeats = count, coverage = coverage,
    coverage_lower = wilson_center - wilson_radius, coverage_upper = wilson_center + wilson_radius,
    mean_length = mean(rows$upper - rows$lower), center_bias = mean(center_error),
    center_rmse = sqrt(mean(center_error^2)))
}))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(records, file.path(output, "paired_intervals.csv"), row.names = FALSE)
write.csv(metrics, file.path(output, "metrics.csv"), row.names = FALSE)
saveRDS(list(input_configuration = configuration, completed_cells = sort(unique(records$cell)),
  reporting_code_md5 = tools::md5sum(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
  file.path(output, "provenance.rds"))
cat("Compared", length(unique(records$cell)), "completed cells with six methods on identical draws\n")
