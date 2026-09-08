#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop(
    paste(
      "usage: audit_cv_tail_equivalence.R",
      "BASELINE_PRIMARY.csv CANDIDATE_PRIMARY.csv",
      "BASELINE_SENSITIVITY.csv CANDIDATE_SENSITIVITY.csv OUTPUT_DIR"
    ),
    call. = FALSE
  )
}

input_paths <- setNames(vapply(args[1:4], normalizePath, character(1L),
                               mustWork = TRUE), c(
  "baseline_primary", "candidate_primary",
  "baseline_sensitivity", "candidate_sensitivity"
))
output_directory <- args[[5L]]
if (file.exists(output_directory)) {
  stop("refusing to overwrite CV-tail audit directory: ", output_directory,
       call. = FALSE)
}
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

results <- lapply(input_paths, utils::read.csv, stringsAsFactors = FALSE)
keys_by_scope <- list(
  primary = c("task_id", "sim_id", "method"),
  sensitivity = c(
    "task_id", "sim_id", "method", "sensitivity_kind",
    "aggregation_cutoff", "M_tau_inference"
  )
)
required <- c(
  "config", "p", "K", "rho", "nlambda_init", "n_bootstrap",
  "package_library", "package_fingerprint", "workflow_fingerprint",
  "manifest_fingerprint", "task_elapsed_seconds"
)

single_value <- function(value, column, label) {
  observed <- unique(value[[column]])
  observed <- observed[!is.na(observed)]
  if (length(observed) != 1L) {
    stop(label, " must contain exactly one nonmissing ", column, ".",
         call. = FALSE)
  }
  observed[[1L]]
}

validate_frame <- function(value, label, keys) {
  missing <- setdiff(c(keys, required), names(value))
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(value[keys])) {
    stop(label, " has duplicate result keys.", call. = FALSE)
  }
  expected <- list(
    config = "C3", p = 100L, K = 4L, rho = 2.5,
    nlambda_init = 100L, n_bootstrap = 5000L
  )
  for (field in names(expected)) {
    if (!identical(
      as.character(single_value(value, field, label)),
      as.character(expected[[field]])
    )) {
      stop(label, " does not use the locked ", field, " setting.",
           call. = FALSE)
    }
  }
  value[do.call(order, unname(value[keys])), , drop = FALSE]
}

for (scope in names(keys_by_scope)) {
  keys <- keys_by_scope[[scope]]
  for (variant in c("baseline", "candidate")) {
    label <- paste(variant, scope, sep = "_")
    results[[label]] <- validate_frame(results[[label]], label, keys)
    rownames(results[[label]]) <- NULL
  }
  if (!identical(
    results[[paste0("baseline_", scope)]][keys],
    results[[paste0("candidate_", scope)]][keys]
  )) {
    stop("baseline and candidate ", scope, " keys differ.", call. = FALSE)
  }
}

digest_columns <- c(
  "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint"
)
read_digests <- function(value, label) {
  setNames(vapply(digest_columns, function(column) {
    digest <- tolower(trimws(as.character(single_value(value, column, label))))
    if (!grepl("^[0-9a-f]{64}$", digest)) {
      stop(label, " has an invalid ", column, ".", call. = FALSE)
    }
    digest
  }, character(1L)), digest_columns)
}
digests <- lapply(names(results), function(label) {
  read_digests(results[[label]], label)
})
names(digests) <- names(results)
for (scope in names(keys_by_scope)) {
  baseline_digest <- digests[[paste0("baseline_", scope)]]
  candidate_digest <- digests[[paste0("candidate_", scope)]]
  if (!identical(
    baseline_digest[["workflow_fingerprint"]],
    candidate_digest[["workflow_fingerprint"]]
  ) || !identical(
    baseline_digest[["manifest_fingerprint"]],
    candidate_digest[["manifest_fingerprint"]]
  )) {
    stop(scope, " workflow or manifest fingerprints differ.", call. = FALSE)
  }
}
if (!identical(digests$baseline_primary, digests$baseline_sensitivity) ||
    !identical(digests$candidate_primary, digests$candidate_sensitivity)) {
  stop("primary and sensitivity provenance differ within a package.",
       call. = FALSE)
}
if (identical(
  digests$baseline_primary[["package_fingerprint"]],
  digests$candidate_primary[["package_fingerprint"]]
)) {
  stop("candidate must use a distinct package fingerprint.", call. = FALSE)
}

is_tail_column <- function(columns) {
  grepl("cv_path_tail_skipped_fold_fits$", columns)
}
is_intentional_count <- function(columns) {
  grepl("cv_invalid_(fold_fits|lambdas)$", columns)
}
is_runtime_or_provenance <- function(columns) {
  columns %in% c("package_library", "package_fingerprint") |
    grepl("^(slurm_|compute_hostname$)", columns) |
    grepl("seconds$", columns)
}

make_mismatches <- function(baseline, candidate, columns, keys) {
  pieces <- lapply(columns, function(column) {
    unequal <- !mapply(
      identical, as.list(baseline[[column]]), as.list(candidate[[column]])
    )
    if (!any(unequal)) return(NULL)
    data.frame(
      baseline[unequal, keys, drop = FALSE],
      column = column,
      baseline = as.character(baseline[[column]][unequal]),
      candidate = as.character(candidate[[column]][unequal]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  pieces <- pieces[!vapply(pieces, is.null, logical(1L))]
  if (length(pieces) == 0L) {
    empty <- as.data.frame(setNames(
      replicate(length(keys), logical(0), simplify = FALSE), keys
    ))
    empty$column <- empty$baseline <- empty$candidate <- character(0)
    return(empty)
  }
  do.call(rbind, pieces)
}

scope_audits <- lapply(names(keys_by_scope), function(scope) {
  keys <- keys_by_scope[[scope]]
  baseline <- results[[paste0("baseline_", scope)]]
  candidate <- results[[paste0("candidate_", scope)]]
  candidate_only <- setdiff(names(candidate), names(baseline))
  expected_candidate_only <- names(candidate)[is_tail_column(names(candidate))]
  if (!setequal(candidate_only, expected_candidate_only) ||
      length(expected_candidate_only) == 0L) {
    stop(scope, " has unexpected candidate-only fields: ",
         paste(setdiff(candidate_only, expected_candidate_only), collapse = ", "),
         call. = FALSE)
  }
  baseline_only <- setdiff(names(baseline), names(candidate))
  if (length(baseline_only) > 0L) {
    stop(scope, " candidate lost baseline fields: ",
         paste(baseline_only, collapse = ", "), call. = FALSE)
  }

  shared <- intersect(names(baseline), names(candidate))
  scientific_columns <- shared[
    !is_runtime_or_provenance(shared) & !is_intentional_count(shared)
  ]
  intentional_columns <- shared[is_intentional_count(shared)]
  forbidden <- make_mismatches(
    baseline, candidate, scientific_columns, keys
  )
  intentional <- make_mismatches(
    baseline, candidate, intentional_columns, keys
  )

  # A skipped candidate is deliberately left at +Inf and therefore can only
  # increase the legacy invalid-count diagnostics. Require every increase to
  # be bounded by the corresponding explicitly reported skipped-fold count.
  for (index in seq_len(nrow(intentional))) {
    column <- intentional$column[[index]]
    row_key <- intentional[index, keys, drop = FALSE]
    row_match <- rep(TRUE, nrow(candidate))
    for (key in keys) row_match <- row_match & candidate[[key]] == row_key[[key]]
    if (sum(row_match) != 1L) {
      stop("could not match an intentional diagnostic row.", call. = FALSE)
    }
    baseline_value <- suppressWarnings(as.numeric(intentional$baseline[[index]]))
    candidate_value <- suppressWarnings(as.numeric(intentional$candidate[[index]]))
    tail_column <- sub(
      "cv_invalid_(fold_fits|lambdas)$",
      "cv_path_tail_skipped_fold_fits",
      column
    )
    tail_value <- suppressWarnings(as.numeric(candidate[[tail_column]][row_match]))
    if (any(!is.finite(c(baseline_value, candidate_value, tail_value))) ||
        candidate_value < baseline_value ||
        candidate_value - baseline_value > tail_value) {
      stop("unreconciled CV diagnostic change in ", scope, ": ", column,
           call. = FALSE)
    }
  }

  list(
    scientific_columns = scientific_columns,
    forbidden = forbidden,
    intentional = intentional,
    tail_columns = expected_candidate_only
  )
})
names(scope_audits) <- names(keys_by_scope)

all_result_keys <- unique(unlist(keys_by_scope, use.names = FALSE))
add_scope <- function(value, scope) {
  missing_keys <- setdiff(all_result_keys, names(value))
  for (key in missing_keys) value[[key]] <- rep(NA, nrow(value))
  value$result_scope <- rep(scope, nrow(value))
  value[c(
    "result_scope", all_result_keys,
    "column", "baseline", "candidate"
  )]
}
forbidden <- do.call(rbind, lapply(names(scope_audits), function(scope) {
  add_scope(scope_audits[[scope]]$forbidden, scope)
}))
intentional <- do.call(rbind, lapply(names(scope_audits), function(scope) {
  add_scope(scope_audits[[scope]]$intentional, scope)
}))

primary_candidate <- results$candidate_primary
aggregate_tail_columns <- grep(
  "^face_(initial_dr|calibrated_dr|calibrated_outcome)_cv_path_tail_skipped_fold_fits$",
  names(primary_candidate), value = TRUE
)
tail_values <- suppressWarnings(as.numeric(unlist(
  primary_candidate[aggregate_tail_columns], use.names = FALSE
)))
tail_values <- tail_values[!is.na(tail_values)]
if (length(tail_values) == 0L || any(!is.finite(tail_values)) ||
    any(tail_values < 0) || any(tail_values != floor(tail_values)) ||
    sum(tail_values) <= 0) {
  stop("candidate did not report valid, exercised path-tail diagnostics.",
       call. = FALSE)
}

baseline_seconds <- as.numeric(single_value(
  results$baseline_primary, "task_elapsed_seconds", "baseline_primary"
))
candidate_seconds <- as.numeric(single_value(
  results$candidate_primary, "task_elapsed_seconds", "candidate_primary"
))
if (any(!is.finite(c(baseline_seconds, candidate_seconds))) ||
    any(c(baseline_seconds, candidate_seconds) <= 0)) {
  stop("task runtimes must be positive and finite.", call. = FALSE)
}
summary <- data.frame(
  primary_rows_compared = nrow(results$candidate_primary),
  sensitivity_rows_compared = nrow(results$candidate_sensitivity),
  primary_scientific_columns_compared =
    length(scope_audits$primary$scientific_columns),
  sensitivity_scientific_columns_compared =
    length(scope_audits$sensitivity$scientific_columns),
  forbidden_mismatches = nrow(forbidden),
  intentional_diagnostic_differences = nrow(intentional),
  path_tail_skipped_fold_fits_total = sum(tail_values),
  baseline_task_seconds = baseline_seconds,
  candidate_task_seconds = candidate_seconds,
  candidate_over_baseline_runtime = candidate_seconds / baseline_seconds,
  speedup = baseline_seconds / candidate_seconds,
  scientific_equivalence_passed = nrow(forbidden) == 0L,
  stringsAsFactors = FALSE
)
write.csv(forbidden, file.path(output_directory, "forbidden_mismatches.csv"),
          row.names = FALSE)
write.csv(
  intentional,
  file.path(output_directory, "intentional_diagnostic_differences.csv"),
  row.names = FALSE
)
write.csv(
  summary,
  file.path(output_directory, "cv_tail_equivalence_summary.csv"),
  row.names = FALSE
)
if (nrow(forbidden) > 0L) {
  stop("CV-tail candidate changed a scientific or final-fit output field.",
       call. = FALSE)
}

gate <- c(
  "cv_tail_scientific_equivalence=passed",
  "config=C3", "p=100", "K=4", "rho=2.5",
  "nlambda_init=100", "n_bootstrap=5000",
  paste0("baseline_package_fingerprint=",
         digests$baseline_primary[["package_fingerprint"]]),
  paste0("candidate_package_fingerprint=",
         digests$candidate_primary[["package_fingerprint"]]),
  paste0("workflow_fingerprint=",
         digests$candidate_primary[["workflow_fingerprint"]]),
  paste0("manifest_fingerprint=",
         digests$candidate_primary[["manifest_fingerprint"]]),
  paste0("baseline_primary_md5=", unname(tools::md5sum(
    input_paths[["baseline_primary"]]
  ))),
  paste0("candidate_primary_md5=", unname(tools::md5sum(
    input_paths[["candidate_primary"]]
  ))),
  paste0("baseline_sensitivity_md5=", unname(tools::md5sum(
    input_paths[["baseline_sensitivity"]]
  ))),
  paste0("candidate_sensitivity_md5=", unname(tools::md5sum(
    input_paths[["candidate_sensitivity"]]
  ))),
  paste0("audit_driver_md5=", unname(tools::md5sum(
    file.path("scripts", "slurm", "audit_cv_tail_equivalence.R")
  )))
)
writeLines(
  gate,
  file.path(output_directory, "cv_tail_scientific_equivalence_passed.txt")
)
print(summary, row.names = FALSE)
message(
  "[pass] CV-tail acceleration preserves every audited scientific and ",
  "final-fit field; intentional invalid-count changes are reconciled."
)
