#!/usr/bin/env Rscript
# Known shared-target factor with unequal source residual variances.
# This remains a summary-model experiment, not fitted high-dimensional RoCE.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: source_file output_directory repeats")
source_file <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(source_file)
output <- arguments[[2L]]
repeats <- as.integer(arguments[[3L]])
check_integer(repeats, "repeats", 100, 100000)
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Output must be on FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output, "cells"), showWarnings = FALSE)

settings <- expand.grid(num_sources = c(4L, 8L, 16L, 32L, 64L),
  shared_variance_unit = c(0, .5, 2), variance_pattern = c("equal", "spread"),
  bias_scale = c(0, 1, 2, 4), stringsAsFactors = FALSE)
settings$invalid_location <- "low_variance"
high_variance <- settings[settings$variance_pattern == "spread" & settings$bias_scale > 0, ]
high_variance$invalid_location <- "high_variance"
settings <- rbind(settings, high_variance)
settings$unshifted_fraction <- .75
settings$assumption_valid <- TRUE
violations <- settings[settings$num_sources %in% c(8L, 32L) & settings$shared_variance_unit == .5 &
                        settings$variance_pattern == "spread" & settings$bias_scale == 2 &
                        settings$invalid_location == "low_variance", ]
violations$unshifted_fraction <- .25
violations$assumption_valid <- FALSE
settings <- rbind(settings, violations)
settings$cell_id <- seq_len(nrow(settings))
settings$seed <- 80300L + settings$cell_id
settings$n_per_site <- 1000L
truth <- .25
alpha <- .05
runner <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
configuration <- list(scope = "Known heteroskedastic shared-target Gaussian reference",
  repeats = repeats, truth = truth, alpha = alpha, assumed_valid_fraction = .75,
  source_md5 = unname(tools::md5sum(source_file)), runner_md5 = unname(tools::md5sum(runner)),
  settings = settings)
configuration_path <- file.path(output, "configuration.rds")
if (file.exists(configuration_path)) {
  stopifnot(identical(configuration, readRDS(configuration_path)))
} else saveRDS(configuration, configuration_path)
write.csv(settings, file.path(output, "manifest.csv"), row.names = FALSE)

summarize_interval <- function(method, limits, covered_set = NULL, fallback = NULL) {
  center <- rowMeans(limits)
  covered <- limits[, 1L] <= truth & truth <= limits[, 2L]
  if (is.null(covered_set)) covered_set <- covered
  if (is.null(fallback)) fallback <- rep(FALSE, repeats)
  data.frame(method = method, repeats = repeats, coverage = mean(covered),
    coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / repeats),
    confidence_set_coverage = mean(covered_set), mean_length = mean(limits[, 2L] - limits[, 1L]),
    center_bias = mean(center - truth), center_rmse = sqrt(mean((center - truth)^2)),
    anchor_fallback_fraction = mean(fallback))
}

calibration_cache <- new.env(parent = emptyenv())
summaries <- list()
for (cell_index in seq_len(nrow(settings))) {
  setting <- settings[cell_index, ]
  checkpoint <- file.path(output, "cells", paste0(setting$cell_id, ".rds"))
  if (file.exists(checkpoint)) {
    entry <- readRDS(checkpoint)
    stopifnot(identical(entry$setting, setting), identical(entry$source_md5, configuration$source_md5),
              entry$repeats == repeats)
    summaries[[cell_index]] <- entry$summary
    next
  }
  count <- setting$num_sources
  valid_minimum <- ceiling(.75 * count)
  votes_required <- floor(count / 2) + 1L
  private_unit <- if (setting$variance_pattern == "equal") rep(1, count) else
    exp(seq(log(.25), log(4), length.out = count))
  private <- private_unit / setting$n_per_site
  shared <- setting$shared_variance_unit / setting$n_per_site
  anchor_private <- 1 / setting$n_per_site
  standard_errors <- sqrt(shared + c(anchor_private, private))
  key <- paste(count, setting$shared_variance_unit, setting$variance_pattern,
               setting$n_per_site, valid_minimum, votes_required, alpha, sep = ":")
  calibrations <- calibration_cache[[key]]
  if (is.null(calibrations)) {
    started <- proc.time()[["elapsed"]]
    calibrations <- list(
      search_marginal_bound = make_source_calibration(count, valid_minimum, votes_required),
      search_known_factor = make_source_calibration(count, valid_minimum, votes_required,
        method = "shared_target_gaussian", shared_variance = shared, source_private_variances = private),
      search_equicorrelation_proxy = make_source_calibration(count, valid_minimum, votes_required,
        method = "equicorrelated_gaussian", source_correlation = shared / (shared + mean(private))))
    calibration_cache[[key]] <- calibrations
    cat("Calibrated", key, "in", proc.time()[["elapsed"]] - started, "seconds\n")
  }
  shifted_count <- count - ceiling(setting$unshifted_fraction * count)
  shifted <- if (setting$invalid_location == "low_variance") seq_len(shifted_count) else
    seq.int(count - shifted_count + 1L, count)
  bias <- numeric(count + 1L)
  bias[shifted + 1L] <- setting$bias_scale * sqrt(anchor_private + private[shifted])
  set.seed(setting$seed)
  common_noise <- sqrt(shared) * rnorm(repeats)
  independent_noise <- sweep(matrix(rnorm(repeats * (count + 1L)), repeats, count + 1L),
                             2L, sqrt(c(anchor_private, private)), "*")
  estimates <- sweep(truth + common_noise + independent_noise, 2L, bias, "+")
  critical <- qnorm(1 - alpha / 2)
  interval_records <- list()
  method_rows <- list()
  add_normal <- function(name, center, se) {
    limits <- cbind(lower = center - critical * se, upper = center + critical * se)
    interval_records[[name]] <<- limits
    method_rows[[name]] <<- summarize_interval(name, limits)
  }
  add_normal("target_only", estimates[, 1L], standard_errors[1L])
  precision <- 1 / c(anchor_private, private)
  pooled <- drop(estimates %*% precision) / sum(precision)
  add_normal("pooled_gls", pooled, sqrt(shared + 1 / sum(precision)))
  valid <- which(bias == 0)
  oracle <- drop(estimates[, valid, drop = FALSE] %*% precision[valid]) / sum(precision[valid])
  add_normal("oracle_gls", oracle, sqrt(shared + 1 / sum(precision[valid])))
  differences <- abs(sweep(estimates[, -1L, drop = FALSE], 1L, estimates[, 1L], "-"))
  included <- cbind(TRUE, sweep(differences, 2L, 2 * sqrt(anchor_private + private), "<=") )
  weights <- sweep(included, 2L, precision, "*")
  naive <- rowSums(estimates * weights) / rowSums(weights)
  add_normal("naive_screen_gls", naive, sqrt(shared + 1 / rowSums(weights)))
  for (name in names(calibrations)) {
    limits <- matrix(NA_real_, repeats, 2L, dimnames = list(NULL, c("lower", "upper")))
    covered_set <- fallback <- logical(repeats)
    for (replicate in seq_len(repeats)) {
      fit <- source_search_interval(estimates[replicate, 1L], standard_errors[1L],
        estimates[replicate, -1L], standard_errors[-1L], calibrations[[name]])
      limits[replicate, ] <- c(fit$lower, fit$upper)
      covered_set[replicate] <- any(fit$intervals[, 1L] <= truth & truth <= fit$intervals[, 2L])
      fallback[replicate] <- fit$used_anchor_fallback
    }
    interval_records[[name]] <- limits
    method_rows[[name]] <- summarize_interval(name, limits, covered_set, fallback)
  }
  summary <- cbind(setting[rep(1L, length(method_rows)), ], do.call(rbind, method_rows))
  summary$true_valid_sources <- length(valid) - 1L
  summary$valid_minimum <- valid_minimum
  summary$votes_required <- votes_required
  rownames(summary) <- NULL
  entry <- list(setting = setting, source_md5 = configuration$source_md5, repeats = repeats,
                summary = summary, calibrations = calibrations, intervals = interval_records)
  temporary <- paste0(checkpoint, ".tmp")
  saveRDS(entry, temporary)
  stopifnot(file.rename(temporary, checkpoint))
  summaries[[cell_index]] <- summary
  write.csv(do.call(rbind, summaries), file.path(output, "metrics.csv"), row.names = FALSE)
  cat("Completed shared-target cell", cell_index, "of", nrow(settings), "at", format(Sys.time()), "\n")
}
write.csv(do.call(rbind, summaries), file.path(output, "metrics.csv"), row.names = FALSE)
writeLines("Known shared-target Gaussian experiment completed; full fitted-score theory remains open.", file.path(output, "COMPLETE"))
