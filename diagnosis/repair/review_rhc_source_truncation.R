#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
study <- normalizePath(arguments[1L], mustWork = TRUE)
reference <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !dir.exists(output))
manifest <- read.csv(file.path(study, "manifest.csv"), stringsAsFactors = FALSE)
old_manifest <- read.csv(file.path(reference, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(manifest) > 0L, all(table(manifest$fold_seed, manifest$radius) == 1L),
  all(manifest$radius %in% 2:5), all(manifest$nuisance_rule == "min"), all(manifest$target_radius == 5))
fits <- signatures <- list()
for (index in seq_len(nrow(manifest))) {
  task <- manifest[index, ]; directory <- file.path(study, "tasks", task$task_id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")))
  fits[[task$config]] <- readRDS(file.path(directory, "checkpoints/analysis.rds"))$value
  signatures[[task$config]] <- readRDS(file.path(directory, "stage_signatures.rds"))
  fitted <- fits[[task$config]]
  stopifnot(fitted$metadata$M_tau == task$radius, fitted$metadata$M_tau_inference == task$radius,
    fitted$metadata$calibration_control$target_radius == 5)
  diagnostics <- read.csv(file.path(directory, "nuisance_fits.csv"))
  checked <- grepl("nonconverged|line_search_failures", names(diagnostics))
  stopifnot(all(as.matrix(diagnostics[, checked, drop = FALSE]) == 0))
}
equal <- function(first, second) stopifnot(isTRUE(all.equal(first, second, tolerance = 1e-10)))
summary <- list()
for (seed in unique(manifest$fold_seed)) {
  name <- paste0("RHC_source_radius5_fold", seed)
  old_task <- old_manifest[old_manifest$config == paste0("RHC_min_fold", seed), ]
  stopifnot(nrow(old_task) == 1L)
  old_directory <- file.path(reference, "tasks", old_task$task_id)
  old_fit <- readRDS(file.path(old_directory, "checkpoints/analysis.rds"))$value
  old_signature <- readRDS(file.path(old_directory, "stage_signatures.rds"))
  baseline <- old_fit
  if (name %in% names(fits)) {
    baseline <- fits[[name]]
    equal(baseline$methods, old_fit$methods)
    equal(signatures[[name]], old_signature)
  }
  for (radius in sort(unique(manifest$radius))) {
    profile <- paste0("RHC_source_radius", radius, "_fold", seed)
    fitted <- fits[[profile]]
    stopifnot(identical(fitted$metadata$data_sha256, baseline$metadata$data_sha256),
              identical(fitted$metadata$fold_sha256, baseline$metadata$fold_sha256))
    for (field in c("target", "initial_outcome_messages")) equal(signatures[[profile]][[field]], old_signature[[field]])
    for (method in c("Target-only", "Calibrated target-only")) equal(
      fitted$methods[as.character(fitted$methods$method) == method, ],
      baseline$methods[as.character(baseline$methods$method) == method, ])
    result <- fitted$tate_fit
    summary[[length(summary)+1L]] <- data.frame(profile = profile, fold_seed = seed,
      source_radius = radius, max_weight = exp(radius), estimate = result$estimate, se = result$se,
      ci_lower = result$ci_lower, ci_upper = result$ci_upper,
      gap_to_anchor = result$estimate - fitted$target_anchor$estimate,
      se_ratio_to_radius5 = result$se/baseline$tate_fit$se)
  }
}
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, summary), file.path(output, "comparisons.csv"), row.names = FALSE)
writeLines(c("SOURCE_TRUNCATION_COMPARISONS_PASSED",
  sprintf("All%d predefined fits complete. Any fitted radius5 matches archived min; target fits and initial OR messages held fixed.", nrow(manifest)),
  "Same source radius in fitting and inference; no point/SE rescaling or case deletion."), file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, summary))
