#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L ||
    any(!startsWith(arguments, "/scratch.global/zhan9381/FACE-HD/"))) stop("Usage: scratch_input scratch_output")
input <- arguments[[1L]]
output <- arguments[[2L]]
stopifnot(file.exists(file.path(input, "COMPLETE")))
configuration <- readRDS(file.path(input, "configuration.rds"))
settings <- configuration$settings
methods <- c("target_only", "pooled_gls", "oracle_gls", "naive_screen_gls",
             "search_marginal_bound", "search_known_correlation", "search_wrong_independence")
baselines <- methods[1:4]
seed_reference <- new.env(parent = emptyenv())
paired_cells <- 0L
metrics <- list()
for (index in seq_len(nrow(settings))) {
  entry <- readRDS(file.path(input, "cells", paste0(settings$cell_id[index], ".rds")))
  stopifnot(identical(entry$setting, settings[index, ]), entry$repeats == configuration$repeats,
            identical(entry$source_md5, configuration$source_md5),
            setequal(names(entry$intervals), methods), nrow(entry$summary) == length(methods))
  for (method in methods) {
    limits <- entry$intervals[[method]]
    row <- entry$summary[entry$summary$method == method, ]
    stopifnot(nrow(limits) == configuration$repeats, nrow(row) == 1L,
              all(is.finite(limits)), all(limits[, 1L] <= limits[, 2L]))
    covered <- sum(limits[, 1L] <= configuration$truth & configuration$truth <= limits[, 2L])
    coverage <- covered / nrow(limits)
    stopifnot(abs(row$coverage - coverage) < 1e-14,
              abs(row$mean_length - mean(limits[, 2L] - limits[, 1L])) < 1e-14,
              row$confidence_set_coverage <= row$coverage + 1e-14)
  }
  if (entry$setting$bias_scale == 0) {
    stopifnot(identical(entry$intervals$oracle_gls, entry$intervals$pooled_gls))
  }
  key <- as.character(entry$setting$seed)
  if (exists(key, seed_reference, inherits = FALSE)) {
    reference <- get(key, seed_reference, inherits = FALSE)
    stopifnot(identical(entry$intervals[baselines], reference))
    paired_cells <- paired_cells + 1L
  } else assign(key, entry$intervals[baselines], seed_reference)
  summary <- entry$summary
  target_length <- summary$mean_length[summary$method == "target_only"]
  oracle_length <- summary$mean_length[summary$method == "oracle_gls"]
  summary$length_ratio_to_target <- summary$mean_length / target_length
  summary$length_ratio_to_oracle <- summary$mean_length / oracle_length
  count <- configuration$repeats
  z <- qnorm(.975)
  denominator <- 1 + z^2 / count
  center <- (summary$coverage + z^2 / (2 * count)) / denominator
  radius <- z * sqrt((summary$coverage * (1 - summary$coverage) + z^2 / (4 * count)) / count) / denominator
  summary$coverage_lower <- center - radius
  summary$coverage_upper <- center + radius
  summary$undercoverage_p <- pbinom(round(summary$coverage * count), count, .95)
  metrics[[index]] <- summary
}
metrics <- do.call(rbind, metrics)
expected <- read.csv(file.path(input, "metrics.csv"), stringsAsFactors = FALSE)
stopifnot(isTRUE(all.equal(metrics[names(expected)], expected, check.attributes = FALSE, tolerance = 1e-12)))
validated <- metrics$assumption_valid & metrics$method %in% c("search_marginal_bound", "search_known_correlation")
metrics$undercoverage_p_holm <- NA_real_
metrics$undercoverage_p_holm[validated] <- p.adjust(metrics$undercoverage_p[validated], "holm")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(metrics, file.path(output, "metrics.csv"), row.names = FALSE)
record <- list(cells = nrow(settings), repeats_per_cell = configuration$repeats,
  metric_rows = nrow(metrics), paired_baseline_cells = paired_cells,
  source_md5 = configuration$source_md5,
  adjusted_undercoverage_flags = sum(metrics$undercoverage_p_holm[validated] < .05))
saveRDS(record, file.path(output, "validation_record.rds"))
writeLines(paste("Validated", nrow(settings), "cells and", nrow(metrics), "metric rows;",
                 paired_cells, "cells exactly match paired baseline intervals."), file.path(output, "COMPLETE"))
print(record)
