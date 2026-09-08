#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 9L) {
  stop(
    paste(
      "usage: summarize_tate_source_diagnostics.R RAW_DIR OUTPUT.csv",
      "PRIMARY_MANIFEST.csv GROUP_MANIFEST.csv CONFIG K",
      "N_REPLICATIONS CUTOFF EXPECT_HARD_THRESHOLD"
    ),
    call. = FALSE
  )
}

raw_directory <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]
primary_manifest_path <- normalizePath(args[[3L]], mustWork = TRUE)
group_manifest_path <- normalizePath(args[[4L]], mustWork = TRUE)
config <- args[[5L]]
source_count <- suppressWarnings(as.integer(args[[6L]]))
expected_replications <- suppressWarnings(as.integer(args[[7L]]))
cutoff <- suppressWarnings(as.numeric(args[[8L]]))
hard_text <- tolower(trimws(args[[9L]]))
if (!config %in% c("C1", "C2", "C3") ||
    is.na(source_count) || !source_count %in% c(2L, 4L, 8L) ||
    is.na(expected_replications) || expected_replications < 1L ||
    expected_replications > 500L || !is.finite(cutoff) || cutoff <= 0 ||
    !hard_text %in% c("true", "false")) {
  stop("invalid diagnostic summary scope.", call. = FALSE)
}
expect_hard <- identical(hard_text, "true")
if (file.exists(output_path)) {
  stop("refusing to overwrite existing summary: ", output_path,
       call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))

runtime <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment("ROCE_WORKFLOW_FINGERPRINT")
manifest_fingerprint <- roce_sha256_file(primary_manifest_path)
group_manifest_fingerprint <- roce_sha256_file(group_manifest_path)
expected_manifest_fingerprint <- roce_sha256_environment(
  "ROCE_MANIFEST_FINGERPRINT"
)
expected_group_manifest_fingerprint <- roce_sha256_environment(
  "ROCE_GROUP_MANIFEST_FINGERPRINT"
)
if (!identical(manifest_fingerprint, expected_manifest_fingerprint) ||
    !identical(
      group_manifest_fingerprint, expected_group_manifest_fingerprint
    )) {
  stop("manifest fingerprints do not match the supplied files.", call. = FALSE)
}

primary <- read.csv(primary_manifest_path, stringsAsFactors = FALSE)
group_manifest <- read.csv(group_manifest_path, stringsAsFactors = FALSE)
required_group_columns <- c(
  "group_task_id", "sim_id", "config", "p", "K", "cutoff",
  "rho_values", "primary_task_ids"
)
if (length(setdiff(required_group_columns, names(group_manifest))) > 0L) {
  stop("group manifest is missing required columns.", call. = FALSE)
}
groups <- group_manifest[
  group_manifest$config == config & group_manifest$p == 100L &
    group_manifest$K == source_count & group_manifest$cutoff == cutoff &
    group_manifest$sim_id %in% seq_len(expected_replications),
  , drop = FALSE
]
groups <- groups[order(groups$sim_id), , drop = FALSE]
if (nrow(groups) != expected_replications ||
    !identical(as.integer(groups$sim_id), seq_len(expected_replications)) ||
    anyDuplicated(groups$group_task_id)) {
  stop("group manifest does not contain the exact requested replications.",
       call. = FALSE)
}
rho_grid <- c(0, 0.5, 1, 1.5, 2, 2.5)
parse_semicolon <- function(value) {
  strsplit(as.character(value), ";", fixed = TRUE)[[1L]]
}
task_ids_by_group <- lapply(seq_len(nrow(groups)), function(index) {
  rhos <- suppressWarnings(as.numeric(parse_semicolon(groups$rho_values[[index]])))
  ids <- suppressWarnings(as.integer(parse_semicolon(
    groups$primary_task_ids[[index]]
  )))
  if (!identical(rhos, rho_grid) || length(ids) != length(rho_grid) ||
      anyNA(ids) || anyDuplicated(ids)) {
    stop("group manifest has an invalid rho/task partition.", call. = FALSE)
  }
  ids
})
expected_task_ids <- unlist(task_ids_by_group, use.names = FALSE)
if (anyDuplicated(expected_task_ids)) {
  stop("group manifest reuses a primary task id.", call. = FALSE)
}

required_primary_columns <- c(
  "task_id", "sim_id", "config", "p", "K", "rho", "cutoff"
)
if (length(setdiff(required_primary_columns, names(primary))) > 0L) {
  stop("primary manifest is missing required columns.", call. = FALSE)
}
expected_primary <- primary[match(expected_task_ids, primary$task_id), , drop = FALSE]
if (anyNA(expected_primary$task_id) || anyDuplicated(primary$task_id) ||
    any(expected_primary$config != config) ||
    any(expected_primary$p != 100L) ||
    any(expected_primary$K != source_count) ||
    any(expected_primary$cutoff != cutoff) ||
    !identical(as.numeric(expected_primary$rho), rep(rho_grid, expected_replications)) ||
    !identical(
      as.integer(expected_primary$sim_id),
      rep(seq_len(expected_replications), each = length(rho_grid))
    )) {
  stop("primary and group manifests do not define the same task partition.",
       call. = FALSE)
}

parse_sentinel <- function(path) {
  if (!file.exists(path)) stop("missing group commit sentinel: ", path,
                              call. = FALSE)
  lines <- readLines(path, warn = FALSE)
  split <- strsplit(lines, "=", fixed = TRUE)
  valid <- lengths(split) >= 2L
  if (!all(valid)) stop("malformed group commit sentinel: ", path,
                        call. = FALSE)
  keys <- vapply(split, `[[`, character(1L), 1L)
  values <- vapply(split, function(parts) paste(parts[-1L], collapse = "="),
                   character(1L))
  if (anyDuplicated(keys)) stop("duplicate sentinel key: ", path,
                                call. = FALSE)
  stats::setNames(values, keys)
}
required_sentinel <- c(
  "rho_group_commit", "group_task_id", "primary_task_ids",
  "package_fingerprint", "workflow_fingerprint", "manifest_fingerprint",
  "group_manifest_fingerprint"
)
for (index in seq_len(nrow(groups))) {
  sentinel <- parse_sentinel(file.path(
    raw_directory, "rho_group_commits",
    sprintf("group_%06d_committed.txt", groups$group_task_id[[index]])
  ))
  if (!all(required_sentinel %in% names(sentinel)) ||
      sentinel[["rho_group_commit"]] != "complete" ||
      suppressWarnings(as.integer(sentinel[["group_task_id"]])) !=
        groups$group_task_id[[index]] ||
      sentinel[["primary_task_ids"]] != groups$primary_task_ids[[index]] ||
      tolower(sentinel[["package_fingerprint"]]) != runtime$fingerprint ||
      tolower(sentinel[["workflow_fingerprint"]]) != workflow_fingerprint ||
      tolower(sentinel[["manifest_fingerprint"]]) != manifest_fingerprint ||
      tolower(sentinel[["group_manifest_fingerprint"]]) !=
        group_manifest_fingerprint) {
    stop("group commit sentinel does not match audited provenance.",
         call. = FALSE)
  }
}

task_files <- file.path(
  raw_directory, sprintf("task_%06d.csv", expected_task_ids)
)
if (any(!file.exists(task_files))) {
  stop("one or more expected task CSVs are missing.", call. = FALSE)
}
results <- RoCE:::.read_simulation_result_files(task_files)
required_result_columns <- c(
  "task_id", "sim_id", "method", "config", "p", "K", "rho", "cutoff",
  "package_library", "package_fingerprint", "workflow_fingerprint",
  "manifest_fingerprint", "rho_group_manifest_fingerprint",
  "rho_group_task_id", "rho_group_size",
  "hard_threshold_diagnostic_requested"
)
if (length(setdiff(required_result_columns, names(results))) > 0L) {
  stop("task results are missing audit metadata.", call. = FALSE)
}
base_methods <- c(
  "one_round_crossfit", "target_only", "sample_size", "inverse_variance",
  "federated_dr", "pooled_dr", "one_round_crossfit_ate_armwise",
  "one_round_crossfit_ate", "one_round_crossfit_ate_quadratic_bias",
  "target_only_ate", "sample_size_ate",
  "inverse_variance_ate", "federated_dr_ate", "pooled_dr_ate"
)
expected_methods <- c(base_methods, if (expect_hard) {
  "one_round_crossfit_ate_hard_threshold"
})
rows_by_task <- split(results, results$task_id)
method_sets_ok <- length(rows_by_task) == length(expected_task_ids) &&
  all(vapply(rows_by_task, function(rows) {
    nrow(rows) == length(expected_methods) && !anyDuplicated(rows$method) &&
      identical(sort(as.character(rows$method)), sort(expected_methods))
  }, logical(1L)))
result_key <- paste(results$sim_id, results$rho, results$method, sep = "|")
provenance <- roce_result_provenance(results, runtime$library)
hard_values <- as.character(results$hard_threshold_diagnostic_requested)
expected_hard_text <- if (expect_hard) "TRUE" else "FALSE"
expected_map <- data.frame(
  task_id = expected_task_ids,
  sim_id = rep(seq_len(expected_replications), each = length(rho_grid)),
  rho = rep(rho_grid, expected_replications),
  rho_group_task_id = rep(
    as.integer(groups$group_task_id), each = length(rho_grid)
  )
)
observed_map <- unique(results[c(
  "task_id", "sim_id", "rho", "rho_group_task_id"
)])
observed_map <- observed_map[match(expected_task_ids, observed_map$task_id), ]
row.names(observed_map) <- NULL
map_ok <- nrow(observed_map) == nrow(expected_map) &&
  identical(as.integer(observed_map$task_id), expected_map$task_id) &&
  identical(as.integer(observed_map$sim_id), expected_map$sim_id) &&
  identical(as.numeric(observed_map$rho), expected_map$rho) &&
  identical(
    as.integer(observed_map$rho_group_task_id),
    expected_map$rho_group_task_id
  )
metadata_ok <- method_sets_ok && !anyDuplicated(result_key) &&
  map_ok &&
  setequal(as.integer(results$task_id), expected_task_ids) &&
  all(results$config == config) && all(results$p == 100L) &&
  all(results$K == source_count) && all(results$cutoff == cutoff) &&
  all(results$rho_group_size == length(rho_grid)) &&
  all(hard_values == expected_hard_text) && provenance$passed &&
  all(tolower(results$workflow_fingerprint) == workflow_fingerprint) &&
  all(tolower(results$manifest_fingerprint) == manifest_fingerprint) &&
  all(tolower(results$rho_group_manifest_fingerprint) ==
        group_manifest_fingerprint)
if (!isTRUE(metadata_ok)) {
  stop("task results fail completeness or provenance validation: ",
       provenance$detail, call. = FALSE)
}

extract_system <- function(data, method, prefix, weight_system) {
  selected <- data[data$method == method, , drop = FALSE]
  if (nrow(selected) == 0L) return(NULL)
  weight_pattern <- paste0("^", prefix, "source_([^_]+)_weight$")
  weight_columns <- grep(weight_pattern, names(selected), value = TRUE)
  if (length(weight_columns) == 0L) {
    stop("missing source-specific weights for ", weight_system,
         call. = FALSE)
  }
  rows <- list()
  for (weight_column in weight_columns) {
    source <- sub(weight_pattern, "\\1", weight_column)
    key <- paste0(prefix, "source_", source, "_")
    required <- paste0(key, c(
      "target_estimate", "source_estimate", "weight", "fold_weight_sd",
      "wald_mean", "wald_max",
      "penalty_activation_fraction", "screening_discrepancy_mean_abs",
      "screening_discrepancy_max_abs",
      "evaluation_discrepancy_mean_abs",
      "evaluation_discrepancy_max_abs", "screening_discrepancy_se_mean",
      "wald_identity_error_max"
    ))
    missing <- setdiff(required, names(selected))
    if (length(missing) > 0L) {
      stop(
        "missing source diagnostics for ", weight_system, "/", source,
        ": ", paste(missing, collapse = ", "), call. = FALSE
      )
    }
    numeric_values <- suppressWarnings(as.numeric(unlist(
      selected[required], use.names = FALSE
    )))
    if (length(numeric_values) != nrow(selected) * length(required) ||
        any(!is.finite(numeric_values))) {
      stop("non-finite source diagnostics for ", weight_system, "/", source,
           call. = FALSE)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      sim_id = selected$sim_id,
      config = selected$config,
      K = selected$K,
      rho = selected$rho,
      source = source,
      weight_system = weight_system,
      target_estimate = selected[[paste0(key, "target_estimate")]],
      source_estimate = selected[[paste0(key, "source_estimate")]],
      weight = selected[[paste0(key, "weight")]],
      fold_weight_sd = selected[[paste0(key, "fold_weight_sd")]],
      wald_mean = selected[[paste0(key, "wald_mean")]],
      wald_max = selected[[paste0(key, "wald_max")]],
      penalty_activation_fraction =
        selected[[paste0(key, "penalty_activation_fraction")]],
      screening_discrepancy_mean_abs = selected[[paste0(
        key, "screening_discrepancy_mean_abs"
      )]],
      screening_discrepancy_max_abs = selected[[paste0(
        key, "screening_discrepancy_max_abs"
      )]],
      evaluation_discrepancy_mean_abs = selected[[paste0(
        key, "evaluation_discrepancy_mean_abs"
      )]],
      evaluation_discrepancy_max_abs = selected[[paste0(
        key, "evaluation_discrepancy_max_abs"
      )]],
      screening_discrepancy_se_mean = selected[[paste0(
        key, "screening_discrepancy_se_mean"
      )]],
      wald_identity_error_max = selected[[paste0(
        key, "wald_identity_error_max"
      )]],
      inclusion_fraction = if (
        paste0(key, "inclusion_fraction") %in% names(selected)
      ) selected[[paste0(key, "inclusion_fraction")]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

long <- do.call(rbind, Filter(Negate(is.null), list(
  extract_system(results, "one_round_crossfit_ate", "", "common_tate"),
  extract_system(
    results, "one_round_crossfit_ate_hard_threshold", "", "hard_tate"
  ),
  extract_system(
    results, "one_round_crossfit_ate_armwise", "mu1_", "arm_mu1"
  ),
  extract_system(
    results, "one_round_crossfit_ate_armwise", "mu0_", "arm_mu0"
  )
)))
if (is.null(long) || nrow(long) == 0L) {
  stop("no source diagnostic rows were extracted.", call. = FALSE)
}
expected_systems <- c("common_tate", "arm_mu1", "arm_mu0", if (expect_hard) {
  "hard_tate"
})
expected_sources <- paste0("s", seq_len(source_count))
observed_system_sources <- split(long$source, long$weight_system)
if (!setequal(names(observed_system_sources), expected_systems) ||
    any(vapply(observed_system_sources, function(values) {
      !setequal(unique(values), expected_sources)
    }, logical(1L))) ||
    any(long$fold_weight_sd < 0) || any(long$wald_mean < 0) ||
    any(long$wald_max < long$wald_mean) ||
    any(long$penalty_activation_fraction < 0 |
          long$penalty_activation_fraction > 1) ||
    any(!is.na(long$inclusion_fraction) &
          (long$inclusion_fraction < 0 | long$inclusion_fraction > 1))) {
  stop("source diagnostic schema or value constraints failed.", call. = FALSE)
}

group_columns <- c("config", "K", "rho", "source", "weight_system")
group_key <- interaction(long[group_columns], drop = TRUE, lex.order = TRUE)
summary_rows <- lapply(split(long, group_key), function(group) {
  data.frame(
    config = group$config[[1L]],
    K = group$K[[1L]],
    rho = group$rho[[1L]],
    source = group$source[[1L]],
    weight_system = group$weight_system[[1L]],
    n_replications = length(unique(group$sim_id)),
    mean_target_estimate = mean(group$target_estimate),
    mean_source_estimate = mean(group$source_estimate),
    mean_weight = mean(group$weight),
    sd_weight = stats::sd(group$weight),
    mean_abs_weight = mean(abs(group$weight)),
    mean_fold_weight_sd = mean(group$fold_weight_sd),
    mean_wald = mean(group$wald_mean),
    max_wald = max(group$wald_max),
    mean_penalty_activation = mean(group$penalty_activation_fraction),
    mean_abs_screening_discrepancy =
      mean(group$screening_discrepancy_mean_abs),
    max_abs_screening_discrepancy =
      max(group$screening_discrepancy_max_abs),
    mean_abs_evaluation_discrepancy =
      mean(group$evaluation_discrepancy_mean_abs),
    max_abs_evaluation_discrepancy =
      max(group$evaluation_discrepancy_max_abs),
    mean_screening_discrepancy_se =
      mean(group$screening_discrepancy_se_mean),
    max_wald_identity_error = max(group$wald_identity_error_max),
    mean_inclusion_fraction = if (all(is.na(group$inclusion_fraction))) {
      NA_real_
    } else {
      mean(group$inclusion_fraction, na.rm = TRUE)
    },
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[do.call(order, summary[group_columns]), , drop = FALSE]
rownames(summary) <- NULL
if (any(summary$n_replications != expected_replications)) {
  stop("one or more diagnostic cells are missing replications.", call. = FALSE)
}
summary$expected_replications <- expected_replications
summary$audited_package_library <- runtime$library
summary$package_fingerprint <- runtime$fingerprint
summary$workflow_fingerprint <- workflow_fingerprint
summary$manifest_fingerprint <- manifest_fingerprint
summary$group_manifest_fingerprint <- group_manifest_fingerprint
summary$hard_threshold_diagnostic_requested <- expect_hard

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile(
  pattern = ".tate_source_diagnostics_", tmpdir = dirname(output_path)
)
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
utils::write.csv(summary, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  stop("failed to commit source diagnostic summary: ", output_path,
       call. = FALSE)
}
message("[done] wrote source diagnostic summary: ", output_path)
