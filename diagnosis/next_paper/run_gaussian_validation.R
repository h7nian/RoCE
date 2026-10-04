#!/usr/bin/env Rscript
# Gaussian-summary experiment, NOT the current paper's full RoCE simulation.
arguments <- commandArgs(trailingOnly = TRUE)
if (!(length(arguments) %in% c(3L, 4L))) stop("Usage: source_file output_directory repeats [profile]")
source_file <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(source_file)
runner <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
settings_file <- file.path(dirname(runner), "gaussian_validation_settings.R")
source(settings_file)
profile <- if (length(arguments) == 4L) arguments[[4L]] else "original"
output <- arguments[[2L]]
repeats <- as.integer(arguments[[3L]])
check_integer(repeats, "repeats", 100, 100000)
if (!startsWith(output, "/scratch.global/zhan9381/FACE-HD/")) stop("Output must be on FACE-HD scratch")
dir.create(output, recursive = TRUE, showWarnings = FALSE)

settings <- make_gaussian_validation_settings(profile)
truth <- .25
alpha <- .05
configuration <- list(scope = "Research Gaussian-summary benchmark, no full nuisance fitting",
  repeats = repeats, truth = truth, alpha = alpha, profile = profile,
  source_md5 = unname(tools::md5sum(source_file)),
  runner_md5 = unname(tools::md5sum(runner)), settings_md5 = unname(tools::md5sum(settings_file)),
  settings = settings)
config_path <- file.path(output, "configuration.rds")
if (file.exists(config_path)) stopifnot(identical(configuration, readRDS(config_path))) else saveRDS(configuration, config_path)
write.csv(settings, file.path(output, "manifest.csv"), row.names = FALSE)
dir.create(file.path(output, "cells"), showWarnings = FALSE)

interval_row <- function(method, lower, upper, estimate, covered_set = NULL, fallback = NULL) {
  covered <- lower <= truth & truth <= upper
  if (is.null(covered_set)) covered_set <- covered
  if (is.null(fallback)) fallback <- rep(FALSE, repeats)
  error <- estimate - truth
  data.frame(method = method, repeats = repeats, coverage = mean(covered),
    coverage_mcse = sqrt(mean(covered) * (1 - mean(covered)) / repeats),
    confidence_set_coverage = mean(covered_set), mean_length = mean(upper - lower),
    median_length = median(upper - lower), center_bias = mean(error),
    center_rmse = sqrt(mean(error^2)), anchor_fallback_fraction = mean(fallback))
}

summary_rows <- list()
for (cell_index in seq_len(nrow(settings))) {
  setting <- settings[cell_index, ]
  checkpoint <- file.path(output, "cells", paste0(setting$cell_id, ".rds"))
  if (file.exists(checkpoint)) {
    entry <- readRDS(checkpoint)
    stopifnot(identical(entry$setting, setting), entry$repeats == repeats,
              identical(entry$source_md5, configuration$source_md5))
    summary_rows[[cell_index]] <- entry$summary
    next
  }
  set.seed(setting$seed)
  count <- setting$num_sources
  valid_minimum <- setting$valid_minimum
  votes_required <- setting$votes_required
  unshifted_count <- ceiling(setting$unshifted_fraction * count)
  correlation <- setting$shared_correlation
  standard_error <- 1 / sqrt(setting$n_per_site)
  discrepancy_se <- sqrt(2 * (1 - correlation)) * standard_error
  bias <- numeric(count + 1L)
  invalid <- seq.int(unshifted_count + 1L, count)
  bias[invalid + 1L] <- setting$bias_scale * discrepancy_se
  if (setting$bias_pattern == "alternating") {
    bias[invalid + 1L] <- bias[invalid + 1L] * rep(c(1, -1), length.out = length(invalid))
  }
  common_noise <- rnorm(repeats)
  independent_noise <- matrix(rnorm(repeats * (count + 1L)), repeats, count + 1L)
  estimates <- truth + standard_error * (sqrt(correlation) * common_noise +
                  sqrt(1 - correlation) * independent_noise)
  estimates <- sweep(estimates, 2L, bias, "+")
  critical <- qnorm(1 - alpha / 2)
  method_rows <- list()
  intervals <- list()

  add_normal <- function(name, center, se) {
    lower <- center - critical * se
    upper <- center + critical * se
    method_rows[[name]] <<- interval_row(name, lower, upper, center)
    intervals[[name]] <<- cbind(lower = lower, upper = upper)
  }
  add_normal("target_only", estimates[, 1L], standard_error)
  pooled <- rowMeans(estimates)
  add_normal("pooled_gls", pooled, standard_error * sqrt(correlation + (1 - correlation) / (count + 1L)))
  true_valid <- which(bias == 0)
  oracle <- rowMeans(estimates[, true_valid, drop = FALSE])
  add_normal("oracle_gls", oracle, standard_error * sqrt(correlation + (1 - correlation) / length(true_valid)))
  included <- cbind(TRUE, abs(sweep(estimates[, -1L, drop = FALSE], 1L, estimates[, 1L], "-")) <= 2 * discrepancy_se)
  included_count <- rowSums(included)
  naive <- rowSums(estimates * included) / included_count
  add_normal("naive_screen_gls", naive, standard_error * sqrt(correlation + (1 - correlation) / included_count))

  calibrations <- list(
    search_marginal_bound = make_source_calibration(count, valid_minimum, votes_required),
    search_known_correlation = make_source_calibration(count, valid_minimum, votes_required,
      method = "equicorrelated_gaussian", source_correlation = correlation),
    search_wrong_independence = make_source_calibration(count, valid_minimum, votes_required,
      method = "equicorrelated_gaussian", source_correlation = 0))
  for (name in names(calibrations)) {
    if (name == "search_wrong_independence" && correlation == 0) {
      method_rows[[name]] <- method_rows[["search_known_correlation"]]
      method_rows[[name]]$method <- name
      intervals[[name]] <- intervals[["search_known_correlation"]]
      next
    }
    limits <- matrix(NA_real_, repeats, 2L, dimnames = list(NULL, c("lower", "upper")))
    covered_set <- fallback <- logical(repeats)
    for (repeat_index in seq_len(repeats)) {
      fit <- source_search_interval(estimates[repeat_index, 1L], standard_error,
        estimates[repeat_index, -1L], rep(standard_error, count), calibrations[[name]])
      limits[repeat_index, ] <- c(fit$lower, fit$upper)
      covered_set[repeat_index] <- any(fit$intervals[, "lower"] <= truth & truth <= fit$intervals[, "upper"])
      fallback[repeat_index] <- fit$used_anchor_fallback
    }
    method_rows[[name]] <- interval_row(name, limits[, 1L], limits[, 2L], rowMeans(limits), covered_set, fallback)
    intervals[[name]] <- limits
  }
  summary <- cbind(setting[rep(1L, length(method_rows)), ], do.call(rbind, method_rows))
  summary$true_valid_sources <- length(true_valid) - 1L
  summary$valid_minimum <- valid_minimum
  summary$votes_required <- votes_required
  rownames(summary) <- NULL
  entry <- list(setting = setting, repeats = repeats, source_md5 = configuration$source_md5,
                summary = summary, intervals = intervals, calibrations = calibrations)
  temporary <- paste0(checkpoint, ".tmp")
  saveRDS(entry, temporary)
  stopifnot(file.rename(temporary, checkpoint))
  summary_rows[[cell_index]] <- summary
  write.csv(do.call(rbind, summary_rows), file.path(output, "metrics.csv"), row.names = FALSE)
  cat("Completed Gaussian cell", cell_index, "of", nrow(settings), "K=", count,
      "correlation=", correlation, "bias_scale=", setting$bias_scale, "at", format(Sys.time()), "\n")
}
write.csv(do.call(rbind, summary_rows), file.path(output, "metrics.csv"), row.names = FALSE)
writeLines("All Gaussian-summary cells completed; this is not validation of the full RoCE extension.", file.path(output, "COMPLETE"))
