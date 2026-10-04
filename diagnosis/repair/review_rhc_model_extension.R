#!/usr/bin/env Rscript
# Combine the retained pilot with disjoint new repeat IDs, without dropping failures.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 1L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
context <- jsonlite::fromJSON(file.path(root, "extension_context.json"))
prior <- normalizePath(context$retained_first_20, mustWork = TRUE)
stopifnot(file.exists(file.path(prior, "CHECKS_PASSED")),
  file.exists(file.path(dirname(configuration$library), "RHC_REFERENCE_PARITY_PASSED")))
original <- read.csv(file.path(prior, "manifest.csv"), stringsAsFactors = FALSE)
additional <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
stopifnot(identical(names(original), names(additional)))
original$input_study <- prior; original$input_task_id <- original$task_id
additional$input_study <- root; additional$input_task_id <- additional$task_id
stopifnot(setequal(original$sim_id, 1:20), setequal(additional$sim_id, 21:200))
index <- rbind(original, additional)
stopifnot(nrow(index) == 800L, !anyDuplicated(index[c("config", "sim_id")]),
  setequal(index$config, c("O_case_mix", "W_overlap", "W_tail", "W_departure")))
for (scenario in unique(index$config)) stopifnot(setequal(index$sim_id[index$config == scenario], 1:200))
index$task_id <- seq_len(nrow(index))
destination <- file.path(root, "combined_mc200")
if (file.exists(file.path(destination, "CHECKS_PASSED"))) quit(status = 0L)
stopifnot(!dir.exists(destination))
combined <- tempfile("combined_pending_", tmpdir = root)
dir.create(file.path(combined, "tasks"), recursive = TRUE)
for (i in seq_len(nrow(index))) {
  directory <- normalizePath(file.path(index$input_study[i], "tasks", index$input_task_id[i]), mustWork = TRUE)
  stopifnot(file.exists(file.path(directory, "COMPLETE")),
    file.exists(file.path(directory, "PAIRED_TARGET_FITS_PASSED")))
  stopifnot(file.symlink(directory, file.path(combined, "tasks", index$task_id[i])))
  index$input_task_directory[i] <- directory
}
write.csv(index, file.path(combined, "result_sources.csv"), row.names = FALSE)
write.csv(index[setdiff(names(index), c("input_study", "input_task_id", "input_task_directory"))],
  file.path(combined, "manifest.csv"), row.names = FALSE)
status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
  shQuote(file.path(root, "review_workflow_v1/review_rhc_model_validation.R")), shQuote(combined),
  shQuote(file.path(combined, "reviews/complete_v1"))))
stopifnot(status == 0L)
writeLines("All four scenarios have all200 original repeat IDs; no failed-fit or seed exclusion; paired radii and convergence checked.",
  file.path(combined, "CHECKS_PASSED"))
stopifnot(file.rename(combined, destination))
