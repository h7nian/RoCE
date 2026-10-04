#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
parent <- normalizePath(configuration$parent_study, mustWork = TRUE)
stopifnot(file.exists(file.path(dirname(configuration$library), "RHC_REFERENCE_PARITY_PASSED")))
source(file.path(root, "workflow/rhc_stage_signatures.R"))
parent_manifest <- read.csv(file.path(parent, "manifest.csv"), stringsAsFactors = FALSE)
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
expected <- parent_manifest[match(manifest$task_id, parent_manifest$task_id), , drop = FALSE]
rownames(expected) <- NULL
stopifnot(identical(manifest, expected))
for (id in manifest$task_id) stopifnot(file.exists(file.path(root, "tasks", id, "COMPLETE")))
parity <- list()
for (id in configuration$control_task_ids) for (radius in configuration$source_radii) {
  old_path <- file.path(parent, "tasks", id, paste0("radius", radius), "fit.rds")
  new_path <- file.path(root, "tasks", id, paste0("radius", radius), "fit.rds")
  old <- readRDS(old_path)$value; new <- readRDS(new_path)$value
  stopifnot(isTRUE(all.equal(rhc_stage_signatures(old), rhc_stage_signatures(new), tolerance = 1e-10)),
    isTRUE(all.equal(old$fold_weights, new$fold_weights, tolerance = 1e-10)))
  difference <- max(abs(c(old$estimate - new$estimate, old$se - new$se)))
  stopifnot(difference < 1e-10)
  parity[[length(parity) + 1L]] <- data.frame(task = id, radius = radius, max_point_se_difference = difference)
}
destination <- file.path(root, "combined")
if (file.exists(file.path(destination, "CHECKS_PASSED"))) quit(status = 0L)
stopifnot(!dir.exists(destination))
combined <- tempfile("combined_pending_", tmpdir = root)
dir.create(file.path(combined, "tasks"), recursive = TRUE)
write.csv(parent_manifest, file.path(combined, "manifest.csv"), row.names = FALSE)
mapping <- list()
for (id in parent_manifest$task_id) {
  origin <- if (id %in% manifest$task_id) root else parent
  directory <- normalizePath(file.path(origin, "tasks", id), mustWork = TRUE)
  stopifnot(file.exists(file.path(directory, "COMPLETE")),
    file.symlink(directory, file.path(combined, "tasks", id)))
  metadata <- jsonlite::fromJSON(file.path(directory, "metadata.json"))
  mapping[[length(mapping) + 1L]] <- data.frame(task_id = id, source_study = origin,
    data_sha256 = metadata$data_sha256,
    disposition = if (id %in% configuration$recovery_task_ids) "recovered" else if
      (id %in% configuration$control_task_ids) "rerun_control" else "unchanged_success")
  if (id %in% configuration$control_task_ids) {
    old_metadata <- jsonlite::fromJSON(file.path(parent, "tasks", id, "metadata.json"))
    stopifnot(identical(metadata$data_sha256, old_metadata$data_sha256))
  }
}
write.csv(do.call(rbind, mapping), file.path(combined, "result_sources.csv"), row.names = FALSE)
write.csv(do.call(rbind, parity), file.path(combined, "successful_path_parity.csv"), row.names = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
  shQuote(file.path(root, "workflow/review_rhc_model_validation.R")), shQuote(combined),
  shQuote(file.path(combined, "reviews/complete_v1"))))
stopifnot(status == 0L)
writeLines("Every prescribed repeat is present once; all original failures recovered; successful-path and real-RHC parity passed.",
  file.path(combined, "CHECKS_PASSED"))
stopifnot(file.rename(combined, destination))
