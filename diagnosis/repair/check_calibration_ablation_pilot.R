#!/usr/bin/env Rscript
# Verify the six production-grid pilot fits against each other and old data.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("Usage: ablation_root baseline_root output_directory")
if (any(!startsWith(arguments, "/scratch.global/zhan9381/FACE-HD/"))) stop("Use FACE-HD scratch")
root <- arguments[[1L]]
baseline <- arguments[[2L]]
output <- arguments[[3L]]
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
manifest <- read.csv(file.path(baseline, "manifest.csv"), stringsAsFactors = FALSE)
baseline_configuration <- jsonlite::fromJSON(file.path(baseline, "configuration.json"))
if (!"p" %in% names(manifest)) manifest$p <- baseline_configuration$p
records <- list()

for (task in seq_len(3L)) {
  directories <- lapply(c("standard_source_p10", "standard_source_target_lasso_p10"),
    function(name) file.path(root, name, "tasks", task))
  stopifnot(all(vapply(directories, function(directory)
    file.exists(file.path(directory, "COMPLETE")), logical(1L))))
  hashes <- vapply(directories, function(directory)
    trimws(readLines(file.path(directory, "data_sha256.txt"))), character(1L))
  stopifnot(hashes[[1L]] == hashes[[2L]])
  results <- lapply(directories, function(directory) readRDS(file.path(directory, "score_derivative.rds")))
  fits <- lapply(results, function(result)
    attr(result, "roce_simulation_artifacts")$direct_tate_results$one_round_crossfit)
  stopifnot(!is.null(fits[[1L]]), !is.null(fits[[2L]]))
  source_error <- max(abs(unlist(fits[[1L]]$source_estimates) - unlist(fits[[2L]]$source_estimates)))
  stopifnot(is.finite(source_error), source_error <= 1e-12,
            fits[[1L]]$crossfit_levels == 3L, fits[[2L]]$crossfit_levels == 2L)
  scenario <- paste0("C", task)
  reference_task <- manifest[manifest$config == scenario & manifest$p == 10L &
    manifest$K == 2L & manifest$rho == 0 & manifest$sim_id == 1L, , drop = FALSE]
  stopifnot(nrow(reference_task) == 1L)
  reference_directory <- file.path(baseline, "tasks", reference_task$task_id)
  stopifnot(file.exists(file.path(reference_directory, "COMPLETE")),
    identical(hashes[[1L]], trimws(readLines(file.path(reference_directory, "data_sha256.txt")))))
  reference <- read.csv(file.path(reference_directory, "results.csv"))
  reference <- reference[reference$method == "target_only_ate", c("estimate", "se", "truth")]
  stopifnot(nrow(reference) == 1L)
  target_error <- max(vapply(results, function(result) {
    target <- result[result$method == "target_only_ate", c("estimate", "se", "truth")]
    stopifnot(nrow(target) == 1L)
    max(abs(unlist(target) - unlist(reference)))
  }, numeric(1L)))
  stopifnot(target_error <= 1e-12)
  for (fit in fits) for (arm in fit$arm_results) for (outer in seq_along(arm$fold_results)) {
    for (source in arm$fold_results[[outer]]$source_results) {
      stopifnot(source$source_nuisance_method == "standard", source$n_calibrated_folds == 0L,
                identical(source$training_folds, setdiff(seq_len(10L), outer)))
      for (name in names(source$inner_fits)) {
        validation <- as.integer(sub("k2_", "", name))
        stopifnot(identical(source$inner_fits[[name]]$training_folds,
                           setdiff(seq_len(10L), c(outer, validation))))
      }
    }
  }
  records[[task]] <- data.frame(config = scenario, source_estimate_error = source_error,
    target_reference_error = target_error, matching_data_hash = hashes[[1L]])
}
write.csv(do.call(rbind, records), file.path(output, "checks.csv"), row.names = FALSE)
writeLines("Six pilot fits passed paired source, target, data and training-boundary checks.",
           file.path(output, "PILOT_CHECKS_PASSED"))
cat("CALIBRATION_ABLATION_PILOTS_PASSED\n")
