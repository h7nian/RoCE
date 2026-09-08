#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "usage: compare_parallel_reproducibility.R",
      "REFERENCE.csv CANDIDATE.csv OUTPUT_DIR"
    ),
    call. = FALSE
  )
}

reference_path <- normalizePath(args[[1L]], mustWork = TRUE)
candidate_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
source(file.path("scripts", "slurm", "resource_topology.R"))

reference <- read.csv(reference_path, stringsAsFactors = FALSE)
candidate <- read.csv(candidate_path, stringsAsFactors = FALSE)
required_columns <- c(
  "task_id", "sim_id", "method", "estimate", "se", "bias", "coverage",
  "workflow_fingerprint", "manifest_fingerprint", "package_library",
  "package_fingerprint", "task_elapsed_seconds",
  roce_resource_metadata_columns()
)
for (label in c("reference", "candidate")) {
  value <- get(label)
  missing_columns <- setdiff(required_columns, names(value))
  if (length(missing_columns) > 0L) {
    stop(
      label, " result is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(value[c("task_id", "sim_id", "method")])) {
    stop(label, " result has duplicated task/seed/method rows.", call. = FALSE)
  }
  roce_validate_result_resource_metadata(value)
}
if (!identical(names(reference), names(candidate))) {
  stop("reference and candidate result schemas differ.", call. = FALSE)
}

row_order <- function(value) {
  order(value$task_id, value$sim_id, value$method)
}
reference <- reference[row_order(reference), , drop = FALSE]
candidate <- candidate[row_order(candidate), , drop = FALSE]
rownames(reference) <- NULL
rownames(candidate) <- NULL
key_columns <- c("task_id", "sim_id", "method")
if (!identical(reference[key_columns], candidate[key_columns])) {
  stop("reference and candidate task/seed/method keys differ.", call. = FALSE)
}

single_value <- function(value, column, label) {
  observed <- unique(trimws(tolower(as.character(value[[column]]))))
  if (length(observed) != 1L || !nzchar(observed)) {
    stop(label, " must contain one nonempty ", column, ".", call. = FALSE)
  }
  observed
}
reference_workflow <- single_value(
  reference, "workflow_fingerprint", "reference"
)
candidate_workflow <- single_value(
  candidate, "workflow_fingerprint", "candidate"
)
reference_manifest <- single_value(
  reference, "manifest_fingerprint", "reference"
)
candidate_manifest <- single_value(
  candidate, "manifest_fingerprint", "candidate"
)
if (!identical(reference_workflow, candidate_workflow) ||
    !identical(reference_manifest, candidate_manifest)) {
  stop(
    "reference and candidate do not share workflow and manifest fingerprints.",
    call. = FALSE
  )
}
for (fingerprint in c(
    single_value(reference, "package_fingerprint", "reference"),
    single_value(candidate, "package_fingerprint", "candidate"),
    reference_workflow, reference_manifest)) {
  if (!grepl("^[0-9a-f]{64}$", fingerprint)) {
    stop("an audited fingerprint is not a lowercase SHA-256 value.",
         call. = FALSE)
  }
}

# Parallel execution is allowed to change only wall/CPU timing fields, the
# declared resource topology, and the immutable package installation used for
# the implementation comparison. Each topology is validated above against its
# own K/fold metadata. All statistical values and diagnostics must retain exact
# double representations.
excluded_columns <- unique(c(
  "package_library", "package_fingerprint",
  grep("^(slurm_|compute_hostname$)", names(reference), value = TRUE),
  roce_resource_metadata_columns(),
  names(reference)[grepl("seconds$", names(reference))]
))
compared_columns <- setdiff(names(reference), excluded_columns)
mismatches <- do.call(rbind, lapply(compared_columns, function(column) {
  unequal <- !mapply(
    identical,
    as.list(reference[[column]]),
    as.list(candidate[[column]])
  )
  if (!any(unequal)) {
    return(NULL)
  }
  data.frame(
    task_id = reference$task_id[unequal],
    sim_id = reference$sim_id[unequal],
    method = reference$method[unequal],
    column = column,
    reference = as.character(reference[[column]][unequal]),
    candidate = as.character(candidate[[column]][unequal]),
    stringsAsFactors = FALSE
  )
}))
if (is.null(mismatches)) {
  mismatches <- data.frame(
    task_id = integer(0), sim_id = integer(0), method = character(0),
    column = character(0), reference = character(0), candidate = character(0),
    stringsAsFactors = FALSE
  )
}

reference_runtime <- unique(reference$task_elapsed_seconds)
candidate_runtime <- unique(candidate$task_elapsed_seconds)
if (length(reference_runtime) != 1L || length(candidate_runtime) != 1L ||
    !is.finite(reference_runtime) || !is.finite(candidate_runtime) ||
    reference_runtime <= 0 || candidate_runtime <= 0) {
  stop("each result must contain one positive finite task runtime.",
       call. = FALSE)
}
summary <- data.frame(
  reference_path = reference_path,
  candidate_path = candidate_path,
  reference_package_fingerprint = single_value(
    reference, "package_fingerprint", "reference"
  ),
  candidate_package_fingerprint = single_value(
    candidate, "package_fingerprint", "candidate"
  ),
  workflow_fingerprint = reference_workflow,
  manifest_fingerprint = reference_manifest,
  rows_compared = nrow(reference),
  columns_compared = length(compared_columns),
  excluded_timing_package_or_resource_columns = length(excluded_columns),
  exact_mismatches = nrow(mismatches),
  reference_runtime_seconds = reference_runtime,
  candidate_runtime_seconds = candidate_runtime,
  candidate_speedup_over_reference = reference_runtime / candidate_runtime,
  parallel_speedup_over_sequential = if (
    isTRUE(unique(reference$parallel_treatment_arms)) &&
      identical(unique(candidate$parallel_treatment_arms), FALSE)
  ) {
    candidate_runtime / reference_runtime
  } else if (
    identical(unique(reference$parallel_treatment_arms), FALSE) &&
      isTRUE(unique(candidate$parallel_treatment_arms))
  ) {
    reference_runtime / candidate_runtime
  } else {
    NA_real_
  },
  exact_reproducibility_passed = nrow(mismatches) == 0L,
  stringsAsFactors = FALSE
)

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
  mismatches, file.path(output_directory, "exact_mismatches.csv")
)
write_atomic_csv(
  summary, file.path(output_directory, "parallel_reproducibility_summary.csv")
)
if (nrow(mismatches) > 0L) {
  stop(
    "parallel candidate changed non-timing result fields; inspect exact_mismatches.csv.",
    call. = FALSE
  )
}

sentinel <- file.path(output_directory, "exact_reproducibility_passed.txt")
temporary_sentinel <- tempfile(
  pattern = ".exact_reproducibility_passed_",
  tmpdir = output_directory, fileext = ".tmp"
)
writeLines(
  c(
    "exact_reproducibility=passed",
    paste0("rows_compared=", nrow(reference)),
    paste0("columns_compared=", length(compared_columns)),
    paste0(
      "candidate_speedup_over_reference=",
      format(summary$candidate_speedup_over_reference, digits = 10)
    ),
    paste0(
      "parallel_speedup_over_sequential=",
      format(summary$parallel_speedup_over_sequential, digits = 10)
    ),
    paste0("workflow_fingerprint=", reference_workflow),
    paste0("manifest_fingerprint=", reference_manifest)
  ),
  temporary_sentinel
)
if (!file.rename(temporary_sentinel, sentinel)) {
  unlink(temporary_sentinel)
  stop("failed to atomically write reproducibility sentinel.", call. = FALSE)
}

print(summary, row.names = FALSE)
message("[pass] scheduling topologies exactly reproduce every non-timing field.")
