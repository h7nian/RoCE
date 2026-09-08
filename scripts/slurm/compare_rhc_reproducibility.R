#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "usage: compare_rhc_reproducibility.R REFERENCE_DIR CANDIDATE_DIR OUTPUT_DIR",
    call. = FALSE
  )
}

reference_directory <- normalizePath(args[[1L]], mustWork = TRUE)
candidate_directory <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
source(file.path("scripts", "slurm", "resource_topology.R"))
resource_columns <- roce_resource_metadata_columns()
planned_outputs <- file.path(
  output_directory,
  c(
    "rhc_reproducibility_summary.csv",
    "rhc_reproducibility_mismatches.csv",
    "exact_reproducibility_passed.txt"
  )
)
if (any(file.exists(planned_outputs))) {
  stop(
    "refusing to overwrite an existing RHC reproducibility audit: ",
    paste(planned_outputs[file.exists(planned_outputs)], collapse = ", "),
    call. = FALSE
  )
}
artifact_names <- c(
  "rhc_direct_tate_methods.csv",
  "rhc_direct_tate_sources.csv",
  "rhc_direct_tate_metadata.csv",
  "rhc_direct_tate_diagnostics.csv",
  "rhc_direct_tate_comparison_diagnostics.csv",
  "rhc_direct_tate_nuisance_fits.csv"
)

read_artifact <- function(directory, artifact) {
  path <- file.path(directory, artifact)
  if (!file.exists(path) || is.na(file.info(path)$size) ||
      file.info(path)$size <= 0) {
    stop("missing or empty RHC artifact: ", path, call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

mismatch_parts <- list()
comparison_rows <- list()
for (artifact in artifact_names) {
  reference <- read_artifact(reference_directory, artifact)
  candidate <- read_artifact(candidate_directory, artifact)
  if (!identical(names(reference), names(candidate))) {
    stop(
      "RHC artifact schema mismatch: ", artifact,
      "; reference columns={", paste(names(reference), collapse = ","),
      "}; candidate columns={", paste(names(candidate), collapse = ","), "}.",
      call. = FALSE
    )
  }
  common_columns <- intersect(names(reference), names(candidate))
  compared_columns <- common_columns[
    !grepl("seconds$", common_columns) &
      !(common_columns %in% c("package_library", resource_columns))
  ]
  if (length(compared_columns) == 0L || nrow(reference) != nrow(candidate)) {
    stop(
      "RHC artifact has no comparable columns or a row-count mismatch: ",
      artifact,
      call. = FALSE
    )
  }
  unequal_rows <- logical(nrow(reference))
  for (column in compared_columns) {
    unequal <- !mapply(
      identical,
      as.list(reference[[column]]),
      as.list(candidate[[column]])
    )
    unequal_rows <- unequal_rows | unequal
    if (any(unequal)) {
      mismatch_parts[[length(mismatch_parts) + 1L]] <- data.frame(
        artifact = artifact,
        row = which(unequal),
        column = column,
        reference = as.character(reference[[column]][unequal]),
        candidate = as.character(candidate[[column]][unequal]),
        stringsAsFactors = FALSE
      )
    }
  }
  comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(
    artifact = artifact,
    rows = nrow(reference),
    columns_compared = length(compared_columns),
    mismatched_rows = sum(unequal_rows),
    stringsAsFactors = FALSE
  )
}

mismatches <- if (length(mismatch_parts) == 0L) {
  data.frame(
    artifact = character(0), row = integer(0), column = character(0),
    reference = character(0), candidate = character(0),
    stringsAsFactors = FALSE
  )
} else {
  do.call(rbind, mismatch_parts)
}
summary <- do.call(rbind, comparison_rows)
summary$exact_reproducibility_passed <- nrow(mismatches) == 0L

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
write_atomic_csv <- function(value, path) {
  temporary <- tempfile(
    pattern = paste0(".", basename(path), "_"),
    tmpdir = dirname(path), fileext = ".tmp"
  )
  write.csv(value, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) {
    unlink(temporary)
    stop("failed to atomically write ", path, call. = FALSE)
  }
}
write_atomic_csv(
  summary, file.path(output_directory, "rhc_reproducibility_summary.csv")
)
write_atomic_csv(
  mismatches, file.path(output_directory, "rhc_reproducibility_mismatches.csv")
)
if (nrow(mismatches) > 0L) {
  stop(
    "RHC candidate changed scientific result fields; inspect mismatch output.",
    call. = FALSE
  )
}

sentinel <- file.path(output_directory, "exact_reproducibility_passed.txt")
temporary_sentinel <- tempfile(
  pattern = ".rhc_exact_reproducibility_",
  tmpdir = output_directory, fileext = ".tmp"
)
writeLines(
  c(
    "exact_reproducibility=passed",
    paste0("artifacts_compared=", nrow(summary)),
    paste0("fields_compared=", sum(summary$rows * summary$columns_compared))
  ),
  temporary_sentinel
)
if (!file.rename(temporary_sentinel, sentinel)) {
  unlink(temporary_sentinel)
  stop("failed to atomically write RHC reproducibility sentinel.", call. = FALSE)
}

print(summary, row.names = FALSE)
message("[pass] RHC candidate exactly reproduces all common scientific fields.")
