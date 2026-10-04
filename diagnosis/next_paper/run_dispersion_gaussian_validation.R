#!/usr/bin/env Rscript
# Extend the completed v6 Gaussian experiment on the same regenerated draws.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: dispersion_source previous_gaussian_run new_scratch_output")
source_file <- normalizePath(arguments[[1L]], mustWork = TRUE)
source(source_file)
previous <- arguments[[2L]]
output <- arguments[[3L]]
if (any(!startsWith(c(previous, output), "/scratch.global/zhan9381/FACE-HD/"))) stop("Use FACE-HD scratch")
if (!file.exists(file.path(previous, "COMPLETE"))) stop("The paired reference study must be complete")
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
dir.create(file.path(output, "cells"))
configuration <- readRDS(file.path(previous, "configuration.rds"))
saveRDS(list(reference = configuration, source_md5 = tools::md5sum(source_file),
  runner_md5 = tools::md5sum(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))),
  file.path(output, "configuration.rds"))
calibrations <- list()
all_metrics <- list()
truth <- configuration$mu1 - configuration$mu0
for (cell in seq_len(nrow(configuration$settings))) {
  reference <- readRDS(file.path(previous, "cells", paste0(cell, ".rds")))
  stopifnot(identical(reference$configuration, configuration), reference$cell == cell)
  setting <- configuration$settings[cell, ]
  count <- setting$num_sources
  width <- count + 1L
  covariance <- reference$covariance
  q1 <- if (setting$profile == "disjoint_quarter") count / 4 else count / 2
  q0 <- if (setting$profile == "treated_half") count else count - q1
  key <- paste(count, setting$profile)
  if (is.null(calibrations[[key]])) {
    calibrations[[key]] <- lapply(c("subset_exact", "uniform_upper"), function(method)
      make_arm_dispersion_calibration(covariance, q1, q0, method))
    names(calibrations[[key]]) <- c("exact", "upper")
  }
  set.seed(68000L + match(count, c(4L, 16L)) * 100L +
    match(setting$profile, c("treated_half", "disjoint_half", "disjoint_quarter")))
  draws <- sweep(matrix(rnorm(configuration$repeats * 2L * width), configuration$repeats, 2L * width) %*%
                   chol(covariance), 2L, reference$population_means, "+")
  reference_target <- reference$records[reference$records$method == "target", ]
  reference_target <- reference_target[order(reference_target$iteration), ]
  target_estimate <- draws[, 1L] - draws[, width + 1L]
  anchor_variance <- covariance[1L, 1L] + covariance[width + 1L, width + 1L] - 2 * covariance[1L, width + 1L]
  target_error <- max(abs(target_estimate - qnorm(.975) * sqrt(anchor_variance) - reference_target$lower),
                      abs(target_estimate + qnorm(.975) * sqrt(anchor_variance) - reference_target$upper))
  design <- cbind(c(rep(1, width), rep(0, width)), c(rep(0, width), rep(1, width)))
  precision_design <- solve(covariance, design)
  contrast <- c(1, -1)
  information <- crossprod(design, precision_design)
  pooled_weights <- drop(precision_design %*% solve(information, contrast))
  pooled_se <- sqrt(drop(crossprod(contrast, solve(information, contrast))))
  pooled <- reference$records[reference$records$method == "pooled_arm_gls", ]
  pooled <- pooled[order(pooled$iteration), ]
  pooled_error <- max(abs(drop(draws %*% pooled_weights) - qnorm(.975) * pooled_se - pooled$lower),
                      abs(drop(draws %*% pooled_weights) + qnorm(.975) * pooled_se - pooled$upper))
  stopifnot(target_error < 1e-12, pooled_error < 1e-12)
  records <- list()
  for (method in c("exact", "upper")) for (anchored in c(FALSE, TRUE)) {
    intervals <- arm_dispersion_intervals(draws, calibrations[[key]][[method]],
      alpha = configuration$alpha, anchor_fraction = if (anchored) .5 else NULL)
    label <- paste0("dispersion_", method, if (anchored) "_anchor" else "")
    records[[label]] <- data.frame(iteration = seq_len(nrow(draws)), method = label,
      lower = intervals$lower, upper = intervals$upper, pooled_estimate = intervals$pooled_estimate,
      bias_allowance = intervals$bias_allowance, anchor_fallback = intervals$used_anchor_fallback)
  }
  records <- do.call(rbind, records)
  target_length <- mean(reference_target$upper - reference_target$lower)
  metrics <- do.call(rbind, lapply(split(records, records$method), function(rows) {
    covered <- rows$lower <= truth & rows$upper >= truth
    ci <- binom.test(sum(covered), nrow(rows))$conf.int
    lengths <- rows$upper - rows$lower
    data.frame(cell = cell, setting, method = rows$method[1L], repeats = nrow(rows),
      coverage = mean(covered), coverage_lower = ci[1L], coverage_upper = ci[2L],
      mean_length = mean(lengths), mean_length_ratio = mean(lengths)/target_length,
      length_ratio_mcse = sd(lengths)/sqrt(length(lengths))/target_length,
      fallback_rate = mean(rows$anchor_fallback), row.names = NULL)
  }))
  saveRDS(list(cell = cell, configuration = configuration, records = records, metrics = metrics,
    target_pairing_error = target_error, pooled_pairing_error = pooled_error,
    data_sha256 = digest::digest(draws, algo = "sha256"), calibration_key = key),
    file.path(output, "cells", paste0(cell, ".rds")))
  all_metrics[[cell]] <- metrics
  write.csv(do.call(rbind, all_metrics), file.path(output, "metrics.csv"), row.names = FALSE)
  cat("Completed paired dispersion cell", cell, "of", nrow(configuration$settings), "\n")
}
saveRDS(calibrations, file.path(output, "calibrations.rds"))
writeLines("Completed the prespecified paired known-Gaussian reference. Fitted RoCE validity remains unproved.",
  file.path(output, "COMPLETE"))
