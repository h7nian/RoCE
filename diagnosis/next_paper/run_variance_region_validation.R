#!/usr/bin/env Rscript
# Gaussian sample-variance uncertainty; no fitted high-dimensional nuisances.
arguments <- commandArgs(trailingOnly = TRUE)
if (!(length(arguments) %in% c(2L, 3L))) stop("Usage: source_directory output_directory [repeats]")
source_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(file.path(source_directory, "source_confidence_sets.R"))
source(file.path(source_directory, "variance_region_calibration.R"))
output <- arguments[[2L]]
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Use FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
settings <- expand.grid(num_sources = c(4L, 16L), n_per_site = c(100L, 1000L))
repeats <- if (length(arguments) == 3L) as.integer(arguments[[3L]]) else 20L
check_integer(repeats, "repeats", 20L, 100000L)
truth <- .25
variance_alpha <- .005
runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
configuration <- list(settings = settings, repeats = repeats, truth = truth, variance_alpha = variance_alpha,
  code_md5 = tools::md5sum(c(runner, file.path(source_directory,
    c("source_confidence_sets.R", "variance_region_calibration.R")))))
configuration_path <- file.path(output, "configuration.rds")
if (file.exists(configuration_path)) {
  stopifnot(identical(configuration, readRDS(configuration_path)))
} else saveRDS(configuration, configuration_path)
dir.create(file.path(output, "cells"), showWarnings = FALSE)
records <- list()
draws <- list()

for (cell in seq_len(nrow(settings))) {
  checkpoint <- file.path(output, "cells", paste0(cell, ".rds"))
  if (file.exists(checkpoint)) {
    saved <- readRDS(checkpoint)
    stopifnot(identical(saved$configuration, configuration), saved$cell == cell)
    records <- c(records, saved$records)
    draws <- c(draws, saved$draws)
    next
  }
  first_record <- length(records) + 1L
  count <- settings$num_sources[cell]
  sample_size <- settings$n_per_site[cell]
  shared <- .5 / sample_size
  anchor_private <- 1 / sample_size
  private <- exp(seq(log(.25), log(4), length.out = count)) / sample_size
  component_variances <- c(shared, anchor_private, private)
  valid_minimum <- ceiling(.75 * count)
  votes_required <- floor(count / 2) + 1L
  invalid <- seq_len(count - valid_minimum)
  valid <- setdiff(seq_len(count), invalid)
  # Here all q inputs are genuinely valid, so the factor CDF is their exact
  # Gaussian vote probability. This provides an independent domination check.
  reference <- make_source_calibration(valid_minimum, valid_minimum, votes_required,
    alpha = .05 - variance_alpha, method = "shared_target_gaussian",
    shared_variance = shared, source_private_variances = private[valid])
  set.seed(91000L + cell)
  for (iteration in seq_len(repeats)) {
    shared_error <- rnorm(1L, sd = sqrt(shared))
    estimates <- truth + shared_error + rnorm(count + 1L, sd = sqrt(c(anchor_private, private)))
    estimates[invalid + 1L] <- estimates[invalid + 1L] + 2 * sqrt(anchor_private + private[invalid])
    variance_estimates <- component_variances * rchisq(count + 2L, sample_size - 1L) / (sample_size - 1L)
    bounds <- gaussian_variance_bounds(variance_estimates, sample_size - 1L, variance_alpha)
    calibration <- make_source_region_calibration(valid_minimum, votes_required,
      bounds[1L, ], bounds[-c(1L, 2L), , drop = FALSE], sum(bounds[1:2, "upper"]),
      variance_alpha = variance_alpha)
    fitted <- source_region_interval(estimates[1L], estimates[-1L], calibration)
    region_valid <- all(bounds[, "lower"] <= component_variances &
                        component_variances <= bounds[, "upper"])
    if (region_valid) {
      stopifnot(calibration$source_critical + 3e-7 >= reference$source_critical,
                all(calibration$source_total_variances[valid] >= shared + private[valid]),
                calibration$anchor_variance_upper >= shared + anchor_private)
    }
    stopifnot(all(is.finite(c(fitted$lower, fitted$upper))), fitted$lower <= fitted$upper)
    index <- length(records) + 1L
    records[[index]] <- data.frame(cell = cell, num_sources = count, n_per_site = sample_size,
      iteration = iteration, region_valid = region_valid, source_critical = calibration$source_critical,
      valid_subset_critical = reference$source_critical, calibration_bound = calibration$calibration_bound,
      lower = fitted$lower, upper = fitted$upper, truth = truth,
      used_anchor_fallback = fitted$used_anchor_fallback,
      numerical_fallback = !is.null(calibration$integration_fallback_reason))
    draws[[index]] <- list(estimates = estimates, variance_estimates = variance_estimates,
                          bounds = bounds, calibration = calibration, interval = fitted)
  }
  indices <- seq.int(first_record, length(records))
  temporary <- paste0(checkpoint, ".tmp")
  saveRDS(list(configuration = configuration, cell = cell,
    records = records[indices], draws = draws[indices]), temporary)
  stopifnot(file.rename(temporary, checkpoint))
  cat("Completed", cell, "of", nrow(settings), "variance-validation cells\n")
}
records <- do.call(rbind, records)
write.csv(records, file.path(output, "pipeline_checks.csv"), row.names = FALSE)
saveRDS(list(settings = settings, repeats = repeats, records = records, draws = draws),
        file.path(output, "pipeline_checks.rds"))
metrics <- do.call(rbind, lapply(split(records, records$cell), function(rows) {
  covered <- rows$lower <= truth & truth <= rows$upper
  coverage <- mean(covered)
  data.frame(cell = rows$cell[1L], num_sources = rows$num_sources[1L],
    n_per_site = rows$n_per_site[1L], repeats = nrow(rows), coverage = coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / nrow(rows)),
    mean_length = mean(rows$upper - rows$lower), region_coverage = mean(rows$region_valid),
    numerical_fallback_fraction = mean(rows$numerical_fallback),
    marginal_calibration_fraction = mean(rows$calibration_bound == "marginal_bound"))
}))
write.csv(metrics, file.path(output, "metrics.csv"), row.names = FALSE)
writeLines("Completed Gaussian variance-region validation; not fitted-score RoCE inference.",
           file.path(output, "COMPLETE"))
cat("VARIANCE_REGION_VALIDATION_PASSED:", nrow(records), "draws;",
    sum(records$region_valid), "simultaneous regions contain the true variances\n")
