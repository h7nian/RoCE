#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    paste(
      "usage: audit_rho_reuse_equivalence.R",
      "INDEPENDENT.csv GROUP_REUSED.csv OUTPUT_DIR"
    ),
    call. = FALSE
  )
}

independent_path <- normalizePath(args[[1L]], mustWork = TRUE)
grouped_path <- normalizePath(args[[2L]], mustWork = TRUE)
output_directory <- args[[3L]]
if (file.exists(output_directory)) {
  stop("refusing to overwrite an existing reuse-audit directory: ",
       output_directory, call. = FALSE)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

independent <- read.csv(independent_path, stringsAsFactors = FALSE)
grouped <- read.csv(grouped_path, stringsAsFactors = FALSE)
required <- c(
  "task_id", "sim_id", "method", "estimate", "se", "bias", "coverage",
  "config", "p", "K", "rho", "nlambda_init", "n_bootstrap",
  "cutoff", "aggregation_cutoff", "aggregation_lambda",
  "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
  "rho_reuse_enabled", "rho_reuse_changed_sources",
  "task_elapsed_seconds"
)
for (label in c("independent", "grouped")) {
  value <- get(label)
  missing <- setdiff(required, names(value))
  if (length(missing) > 0L) {
    stop(label, " result is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(value[c("task_id", "sim_id", "method")])) {
    stop(label, " result contains duplicate task/seed/method rows.",
         call. = FALSE)
  }
}
group_required <- c(
  "rho_group_task_id", "rho_group_size",
  "rho_group_manifest_fingerprint", "rho_group_positive_rho_workers",
  "rho_group_positive_rho_backend", "rho_group_source_workers_per_fit"
)
if (!all(group_required %in% names(grouped))) {
  stop("grouped result lacks rho-group provenance fields.", call. = FALSE)
}

row_order <- function(value) order(value$task_id, value$sim_id, value$method)
independent <- independent[row_order(independent), , drop = FALSE]
grouped <- grouped[row_order(grouped), , drop = FALSE]
rownames(independent) <- rownames(grouped) <- NULL
keys <- c("task_id", "sim_id", "method")
if (!identical(independent[keys], grouped[keys])) {
  stop("independent and grouped result keys differ.", call. = FALSE)
}

single_value <- function(value, column, label) {
  observed <- unique(value[[column]])
  if (length(observed) != 1L || is.na(observed[[1L]])) {
    stop(label, " must contain one nonmissing ", column, ".", call. = FALSE)
  }
  observed[[1L]]
}
expected <- list(
  config = "C3", p = 100L, K = 4L, rho = 2.5,
  nlambda_init = 100L, n_bootstrap = 5000L
)
for (field in names(expected)) {
  for (label in c("independent", "grouped")) {
    if (!identical(
      as.character(single_value(get(label), field, label)),
      as.character(expected[[field]])
    )) {
      stop(label, " does not use the locked ", field, " setting.",
           call. = FALSE)
    }
  }
}
primary_cutoff <- suppressWarnings(as.numeric(single_value(
  independent, "cutoff", "independent"
)))
if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
    primary_cutoff <= 0) {
  stop("independent result has no valid primary cutoff.", call. = FALSE)
}
for (label in c("independent", "grouped")) {
  value <- get(label)
  cutoff_consistent <-
    abs(as.numeric(single_value(value, "cutoff", label)) -
          primary_cutoff) <= 1e-12 &&
    abs(as.numeric(single_value(value, "aggregation_cutoff", label)) -
          primary_cutoff) <= 1e-12 &&
    abs(as.numeric(single_value(value, "aggregation_lambda", label)) -
          1 / primary_cutoff) <= 1e-12
  if (!isTRUE(cutoff_consistent)) {
    stop(label, " aggregation cutoff/lambda metadata are inconsistent.",
         call. = FALSE)
  }
}
# This audit compares the positive-rho endpoint file. Its update deliberately
# uses one source worker per PSOCK process; the separate rho=0 file records the
# four-source reference topology and is covered by the all-or-none group commit.
if (any(as.logical(independent$rho_reuse_enabled)) ||
    !all(as.logical(grouped$rho_reuse_enabled)) ||
    !identical(unique(grouped$rho_reuse_changed_sources), "s1") ||
    !identical(single_value(grouped, "rho_group_size", "grouped"), 6L) ||
    !identical(single_value(
      grouped, "rho_group_positive_rho_workers", "grouped"
    ), 5L) ||
    !identical(single_value(
      grouped, "rho_group_positive_rho_backend", "grouped"
    ), "psock") ||
    !identical(single_value(
      grouped, "rho_group_source_workers_per_fit", "grouped"
    ), 1L)) {
  stop("independent/grouped reuse metadata is inconsistent.", call. = FALSE)
}

fingerprint_columns <- c(
  "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint"
)
fingerprints <- setNames(character(length(fingerprint_columns)),
                         fingerprint_columns)
for (column in fingerprint_columns) {
  independent_value <- tolower(trimws(as.character(
    single_value(independent, column, "independent")
  )))
  grouped_value <- tolower(trimws(as.character(
    single_value(grouped, column, "grouped")
  )))
  if (!identical(independent_value, grouped_value) ||
      !grepl("^[0-9a-f]{64}$", independent_value)) {
    stop("invalid or mismatched ", column, ".", call. = FALSE)
  }
  fingerprints[[column]] <- independent_value
}
group_manifest_fingerprint <- tolower(trimws(as.character(
  single_value(grouped, "rho_group_manifest_fingerprint", "grouped")
)))
if (!grepl("^[0-9a-f]{64}$", group_manifest_fingerprint)) {
  stop("invalid rho-group manifest fingerprint.", call. = FALSE)
}

group_only_columns <- setdiff(names(grouped), names(independent))
if (!setequal(group_only_columns, group_required)) {
  stop(
    "unexpected grouped-only columns: ",
    paste(setdiff(group_only_columns, group_required), collapse = ", "),
    call. = FALSE
  )
}
independent_only_columns <- setdiff(names(independent), names(grouped))
if (length(independent_only_columns) > 0L) {
  stop("grouped result lost independent columns: ",
       paste(independent_only_columns, collapse = ", "), call. = FALSE)
}
excluded <- unique(c(
  "package_library", "rho_reuse_enabled", "rho_reuse_changed_sources",
  grep("^(slurm_|compute_hostname$)", names(independent), value = TRUE),
  names(independent)[grepl("seconds$", names(independent))]
))
compared_columns <- setdiff(intersect(names(independent), names(grouped)), excluded)
mismatches <- do.call(rbind, lapply(compared_columns, function(column) {
  unequal <- !mapply(
    identical,
    as.list(independent[[column]]),
    as.list(grouped[[column]])
  )
  if (!any(unequal)) return(NULL)
  data.frame(
    task_id = independent$task_id[unequal],
    method = independent$method[unequal],
    column = column,
    independent = as.character(independent[[column]][unequal]),
    grouped = as.character(grouped[[column]][unequal]),
    stringsAsFactors = FALSE
  )
}))
if (is.null(mismatches)) {
  mismatches <- data.frame(
    task_id = integer(0), method = character(0), column = character(0),
    independent = character(0), grouped = character(0),
    stringsAsFactors = FALSE
  )
}

independent_seconds <- as.numeric(single_value(
  independent, "task_elapsed_seconds", "independent"
))
group_seconds <- as.numeric(single_value(
  grouped, "task_elapsed_seconds", "grouped"
))
if (!is.finite(independent_seconds) || independent_seconds <= 0 ||
    !is.finite(group_seconds) || group_seconds <= 0) {
  stop("audit runtimes must be positive and finite.", call. = FALSE)
}
summary <- data.frame(
  independent_path = independent_path,
  grouped_path = grouped_path,
  task_id = single_value(grouped, "task_id", "grouped"),
  group_task_id = single_value(
    grouped, "rho_group_task_id", "grouped"
  ),
  rows_compared = nrow(grouped),
  columns_compared = length(compared_columns),
  exact_mismatches = nrow(mismatches),
  independent_task_seconds = independent_seconds,
  six_rho_group_seconds = group_seconds,
  approximate_six_independent_over_group_speedup =
    6 * independent_seconds / group_seconds,
  exact_equivalence_passed = nrow(mismatches) == 0L,
  stringsAsFactors = FALSE
)

write.csv(mismatches, file.path(output_directory, "exact_mismatches.csv"),
          row.names = FALSE)
write.csv(summary, file.path(output_directory, "reuse_equivalence_summary.csv"),
          row.names = FALSE)
if (nrow(mismatches) > 0L) {
  stop("rho reuse changed non-timing statistical fields.", call. = FALSE)
}

project_root <- normalizePath(".", mustWork = TRUE)
md5_line <- function(name, relative_path) {
  paste0(name, "=", unname(tools::md5sum(file.path(
    project_root, relative_path
  ))))
}
gate_lines <- c(
  "rho_reuse_equivalence=passed",
  paste0("task_id=", summary$task_id),
  paste0("group_task_id=", summary$group_task_id),
  "config=C3", "p=100", "K=4", "rho=2.5",
  "nlambda_init=100", "n_bootstrap=5000",
  paste0("primary_cutoff=", primary_cutoff),
  "positive_rho_workers=5", "positive_rho_backend=psock",
  paste0("package_fingerprint=", fingerprints[["package_fingerprint"]]),
  paste0("workflow_fingerprint=", fingerprints[["workflow_fingerprint"]]),
  paste0("manifest_fingerprint=", fingerprints[["manifest_fingerprint"]]),
  paste0("group_manifest_fingerprint=", group_manifest_fingerprint),
  paste0("independent_result_md5=", unname(tools::md5sum(independent_path))),
  paste0("grouped_result_md5=", unname(tools::md5sum(grouped_path))),
  md5_line(
    "reuse_audit_driver_md5",
    "scripts/slurm/audit_rho_reuse_equivalence.R"
  ),
  md5_line(
    "rho_group_task_driver_md5",
    "scripts/slurm/run_direct_tate_rho_group_task.R"
  ),
  md5_line(
    "rho_group_array_driver_md5",
    "scripts/slurm/run_direct_tate_rho_group_array.sh"
  ),
  md5_line(
    "task_helpers_md5",
    "scripts/slurm/direct_tate_task_helpers.R"
  ),
  md5_line(
    "checkpoint_audit_driver_md5",
    "scripts/slurm/audit_direct_tate_checkpoint.R"
  ),
  md5_line(
    "simulation_qc_policy_md5",
    "scripts/slurm/simulation_qc_policy.R"
  ),
  md5_line(
    "group_checkpoint_driver_md5",
    "scripts/slurm/run_rho_group_checkpoint_audit.sh"
  ),
  md5_line(
    "rho_group_submit_driver_md5",
    "scripts/slurm/submit_rho_group_direct_tate.sh"
  )
)
gate_path <- file.path(output_directory, "rho_reuse_equivalence_passed.txt")
writeLines(gate_lines, gate_path)
print(summary, row.names = FALSE)
message("[pass] p=100 grouped rho reuse exactly matches the independent fit.")
