#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/manifest_main.csv"
}
output_directory <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  dirname(input_path)
}

manifest <- read.csv(input_path, stringsAsFactors = FALSE)
required_columns <- c("task_id", "K")
missing_columns <- setdiff(required_columns, names(manifest))
if (length(missing_columns) > 0L) {
  stop(
    "manifest is missing required column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}
if (anyDuplicated(manifest$task_id)) {
  stop("manifest task_id values must be unique.")
}
if (nrow(manifest) == 0L) {
  stop("manifest contains no tasks.")
}

source_counts <- sort(unique(as.integer(manifest$K)))
if (any(is.na(source_counts)) || any(source_counts < 1L)) {
  stop("manifest K values must be positive integers.")
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_paths <- file.path(
  output_directory, sprintf("manifest_main_K%d.csv", source_counts)
)
if (any(file.exists(output_paths))) {
  stop(
    "refusing to overwrite existing K-specific manifest(s): ",
    paste(output_paths[file.exists(output_paths)], collapse = ", "),
    call. = FALSE
  )
}
for (source_count in source_counts) {
  subset_manifest <- manifest[manifest$K == source_count, , drop = FALSE]
  output_path <- file.path(
    output_directory,
    sprintf("manifest_main_K%d.csv", source_count)
  )
  temporary_path <- tempfile(
    pattern = sprintf(".manifest_main_K%d_", source_count),
    tmpdir = output_directory, fileext = ".csv"
  )
  write.csv(subset_manifest, temporary_path, row.names = FALSE)
  if (!file.rename(temporary_path, output_path)) {
    unlink(temporary_path, force = TRUE)
    stop("failed to atomically commit K-specific manifest: ", output_path,
         call. = FALSE)
  }
  message(
    "wrote ", output_path, " with ", nrow(subset_manifest),
    " tasks (K=", source_count, " source sites)"
  )
}
