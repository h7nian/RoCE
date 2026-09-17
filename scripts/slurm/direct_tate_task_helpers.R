roce_parse_method_list <- function(methods) {
  values <- strsplit(as.character(methods), ",", fixed = TRUE)[[1L]]
  trimws(values[nzchar(trimws(values))])
}

# The manifest's deviation mechanism for one task/group row (HISTORY #0009):
# required, and restricted to the mechanisms the FACE DGP implements.
roce_task_deviation_mechanism <- function(task) {
  mechanism <- as.character(task$deviation_mechanism)
  if (length(mechanism) != 1L || is.na(mechanism) ||
      !mechanism %in% RoCE:::.face_deviation_mechanisms()) {
    stop("manifest row lacks a valid deviation_mechanism.", call. = FALSE)
  }
  mechanism
}

roce_parse_boolean <- function(value, name) {
  normalized <- tolower(trimws(as.character(value)))
  if (length(normalized) != 1L ||
      !normalized %in% c("true", "false", "1", "0")) {
    stop(name, " must be TRUE/FALSE or 1/0.", call. = FALSE)
  }
  normalized %in% c("true", "1")
}

roce_claim_task_lock <- function(lock_path, claim_lines) {
  if (length(lock_path) != 1L || is.na(lock_path) || !nzchar(lock_path)) {
    stop("lock_path must be one nonempty path.", call. = FALSE)
  }
  if (length(claim_lines) == 0L || anyNA(claim_lines)) {
    stop("claim_lines must contain at least one nonmissing value.",
         call. = FALSE)
  }
  dir.create(dirname(lock_path), recursive = TRUE, showWarnings = FALSE)
  if (!dir.create(lock_path, showWarnings = FALSE)) {
    stop("task is already claimed by another worker: ", lock_path,
         call. = FALSE)
  }
  claim_path <- file.path(lock_path, "claim.txt")
  tryCatch(
    writeLines(as.character(claim_lines), claim_path),
    error = function(error) {
      unlink(lock_path, recursive = TRUE, force = TRUE)
      stop(error)
    }
  )
  function() {
    unlink(lock_path, recursive = TRUE, force = TRUE)
    invisible(NULL)
  }
}

roce_scheduler_provenance <- function() {
  value_or_local <- function(name) {
    value <- trimws(Sys.getenv(name, ""))
    if (nzchar(value)) value else "local"
  }
  list(
    slurm_job_id = value_or_local("SLURM_JOB_ID"),
    slurm_array_job_id = value_or_local("SLURM_ARRAY_JOB_ID"),
    slurm_array_task_id = value_or_local("SLURM_ARRAY_TASK_ID"),
    slurm_cluster_name = value_or_local("SLURM_CLUSTER_NAME"),
    slurm_job_partition = value_or_local("SLURM_JOB_PARTITION"),
    slurm_node_list = value_or_local("SLURM_JOB_NODELIST"),
    compute_hostname = value_or_local("HOSTNAME")
  )
}

roce_scheduler_provenance_columns <- function() {
  c(
    "slurm_job_id", "slurm_array_job_id", "slurm_array_task_id",
    "slurm_cluster_name", "slurm_job_partition", "slurm_node_list",
    "compute_hostname"
  )
}

roce_validate_scheduler_provenance <- function(
    result, require_slurm = FALSE) {
  if (!is.data.frame(result) || nrow(result) < 1L) {
    stop("result must be a nonempty data frame.", call. = FALSE)
  }
  columns <- roce_scheduler_provenance_columns()
  missing <- setdiff(columns, names(result))
  if (length(missing) > 0L) {
    stop("result is missing scheduler provenance: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  normalized <- result[columns]
  normalized[] <- lapply(normalized, function(values) {
    trimws(as.character(values))
  })
  if (anyNA(normalized) || any(vapply(normalized, function(values) {
    any(!nzchar(values))
  }, logical(1L)))) {
    stop("scheduler provenance must be nonmissing and nonempty.",
         call. = FALSE)
  }
  if (isTRUE(require_slurm) && any(vapply(normalized, function(values) {
    any(values == "local")
  }, logical(1L)))) {
    stop("production scheduler provenance cannot use local placeholders.",
         call. = FALSE)
  }
  invisible(TRUE)
}

roce_default_sensitivity_grid <- function(
    reference_cutoff = 1 / RoCE:::AGG_WALD_LAMBDA) {
  reference_cutoff <- suppressWarnings(as.numeric(reference_cutoff))
  if (length(reference_cutoff) != 1L || !is.finite(reference_cutoff) ||
      reference_cutoff <= 0) {
    stop("reference_cutoff must be one positive finite number.",
         call. = FALSE)
  }
  cutoff_grid <- unique(c(1, 1.5, 2, 2.5, 3, reference_cutoff))
  rbind(
    data.frame(
      cutoff = cutoff_grid,
      M_tau_inference = 5
    ),
    data.frame(
      cutoff = reference_cutoff,
      M_tau_inference = c(4, 6, Inf)
    )
  )
}

roce_is_reused_sensitivity_task <- function(
    task, primary_cutoff = 1 / RoCE:::AGG_WALD_LAMBDA) {
  required <- c(
    "experiment", "config", "K", "rho", "cutoff", "M_tau",
    "M_tau_inference"
  )
  if (!is.data.frame(task) || nrow(task) != 1L ||
      length(setdiff(required, names(task))) > 0L) {
    stop("task must be one row with the complete sensitivity setting.",
         call. = FALSE)
  }
  primary_cutoff <- suppressWarnings(as.numeric(primary_cutoff))
  if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
      primary_cutoff <= 0) {
    stop("primary_cutoff must be one positive finite number.",
         call. = FALSE)
  }
  as.character(task$experiment) %in%
      c("negative_transfer", "single_task_smoke") &&
    identical(as.character(task$deviation_mechanism), "treated_arm") &&
    identical(as.character(task$config), "C3") &&
    as.integer(task$K) == 4L &&
    as.numeric(task$rho) %in% c(0, 2.5) &&
    abs(as.numeric(task$cutoff) - primary_cutoff) <= 1e-12 &&
    as.numeric(task$M_tau) == 5 &&
    as.numeric(task$M_tau_inference) == 5
}

roce_validate_sensitivity_identity <- function(
    sensitivity_result, primary_result, tolerance = 1e-12) {
  reference <- sensitivity_result[
    sensitivity_result$sensitivity_kind == "reference_identity",
    , drop = FALSE
  ]
  primary_reference <- primary_result[
    primary_result$method %in% c(
      "one_round_crossfit_ate", "target_only_ate"
    ),
    , drop = FALSE
  ]
  reference <- reference[order(reference$method), , drop = FALSE]
  primary_reference <- primary_reference[
    order(primary_reference$method), , drop = FALSE
  ]
  identity_columns <- c(
    "method", "estimate", "se", "bias", "coverage", "ci_width"
  )
  if (nrow(reference) != 2L || nrow(primary_reference) != 2L ||
      !all(identity_columns %in% names(reference)) ||
      !all(identity_columns %in% names(primary_reference)) ||
      !isTRUE(all.equal(
        reference[identity_columns],
        primary_reference[identity_columns],
        tolerance = tolerance,
        check.attributes = FALSE
      ))) {
    stop(
      "default reused sensitivity did not reproduce the primary TATE rows.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

roce_annotate_direct_tate_rows <- function(
    result, task, resource_plan, task_elapsed_seconds,
    package_provenance, workflow_fingerprint, manifest_fingerprint,
    metadata_policy = c("from_task", "preserve")) {
  metadata_policy <- match.arg(metadata_policy)
  if (!is.data.frame(result) || nrow(result) < 1L ||
      !is.data.frame(task) || nrow(task) != 1L) {
    stop("result must be nonempty and task must contain one row.",
         call. = FALSE)
  }
  assert_matches_task <- function(field, expected_field = field) {
    if (!field %in% names(result) || !expected_field %in% names(task)) {
      stop("result/task metadata field is missing: ", field,
           call. = FALSE)
    }
    observed <- result[[field]]
    expected <- task[[expected_field]][[1L]]
    matches <- if (is.numeric(expected) || is.integer(expected)) {
      observed_numeric <- suppressWarnings(as.numeric(observed))
      expected_numeric <- suppressWarnings(as.numeric(expected))
      all(is.finite(observed_numeric)) && is.finite(expected_numeric) &&
        all(abs(observed_numeric - expected_numeric) <= 1e-12)
    } else {
      all(!is.na(observed)) &&
        all(as.character(observed) == as.character(expected))
    }
    if (!isTRUE(matches)) {
      stop("result metadata '", field, "' does not match task.",
           call. = FALSE)
    }
    invisible(TRUE)
  }
  for (field in c("sim_id", "p", "K", "config", "deviation_mechanism")) {
    assert_matches_task(field)
  }

  if (identical(metadata_policy, "from_task")) {
    if ("rho" %in% names(result)) assert_matches_task("rho")
    task_cutoff <- suppressWarnings(as.numeric(task$cutoff))
    if (length(task_cutoff) != 1L || !is.finite(task_cutoff) ||
        task_cutoff <= 0) {
      stop("task cutoff must be one positive finite number.", call. = FALSE)
    }
    if ("cutoff" %in% names(result)) {
      observed <- suppressWarnings(as.numeric(result$cutoff))
      if (anyNA(observed) || any(!is.finite(observed)) ||
          any(abs(observed - task_cutoff) > 1e-12)) {
        stop("result metadata 'cutoff' does not match task.", call. = FALSE)
      }
    }
    if ("aggregation_cutoff" %in% names(result)) {
      observed <- suppressWarnings(as.numeric(result$aggregation_cutoff))
      if (anyNA(observed) || any(!is.finite(observed)) ||
          any(abs(observed - task_cutoff) > 1e-12)) {
        stop("result aggregation cutoff does not match task.", call. = FALSE)
      }
    }
    if ("aggregation_lambda" %in% names(result)) {
      observed <- suppressWarnings(as.numeric(result$aggregation_lambda))
      if (anyNA(observed) || any(!is.finite(observed)) ||
          any(abs(observed - 1 / task_cutoff) > 1e-12)) {
        stop("result aggregation lambda does not match task.", call. = FALSE)
      }
    }
    result$task_id <- as.integer(task$task_id)
    result$experiment <- as.character(task$experiment)
    result$rho <- as.numeric(task$rho)
    result$cutoff <- task_cutoff
    result$aggregation_cutoff <- task_cutoff
    result$aggregation_lambda <- 1 / task_cutoff
    result$primary_cutoff <- task_cutoff
    result$n_site <- as.integer(task$n_site)
    result$n_folds <- as.integer(task$n_folds)
  } else {
    required_preserved <- c(
      "task_id", "experiment", "primary_experiment", "rho", "cutoff",
      "aggregation_cutoff", "aggregation_lambda", "primary_cutoff",
      "n_site", "n_folds"
    )
    missing_preserved <- setdiff(required_preserved, names(result))
    if (length(missing_preserved) > 0L) {
      stop("preserved sidecar metadata is missing: ",
           paste(missing_preserved, collapse = ", "), call. = FALSE)
    }
    for (field in c("task_id", "rho", "n_site", "n_folds")) {
      assert_matches_task(field)
    }
    assert_matches_task("primary_experiment", "experiment")
    sidecar_cutoff <- suppressWarnings(as.numeric(result$cutoff))
    aggregation_cutoff <- suppressWarnings(as.numeric(
      result$aggregation_cutoff
    ))
    aggregation_lambda <- suppressWarnings(as.numeric(
      result$aggregation_lambda
    ))
    preserved_primary_cutoff <- suppressWarnings(as.numeric(
      result$primary_cutoff
    ))
    task_primary_cutoff <- suppressWarnings(as.numeric(task$cutoff))
    preserved_numeric <- c(
      sidecar_cutoff, aggregation_cutoff, aggregation_lambda,
      preserved_primary_cutoff, task_primary_cutoff
    )
    preserved_labels <- c(
      as.character(result$experiment),
      as.character(result$primary_experiment)
    )
    if (anyNA(preserved_labels) ||
        any(result$experiment != "c3_reused_sensitivity") ||
        anyNA(preserved_numeric) || any(!is.finite(preserved_numeric)) ||
        any(sidecar_cutoff <= 0) || task_primary_cutoff <= 0 ||
        any(abs(preserved_primary_cutoff - task_primary_cutoff) >
              1e-12) ||
        any(abs(aggregation_cutoff - sidecar_cutoff) > 1e-12) ||
        any(abs(aggregation_lambda - 1 / sidecar_cutoff) > 1e-12)) {
      stop("preserved sensitivity metadata is internally inconsistent.",
           call. = FALSE)
    }
  }

  if (length(task_elapsed_seconds) != 1L ||
      !is.finite(task_elapsed_seconds) || task_elapsed_seconds < 0) {
    stop("task_elapsed_seconds must be one finite nonnegative number.",
         call. = FALSE)
  }
  result$allocated_cores <- resource_plan$allocated_cores
  result$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
  result$nuisance_solver <- RoCE:::nuisance_solver_cpp()
  result$source_workers_per_arm <- resource_plan$source_workers
  result$parallel_treatment_arms <- resource_plan$parallel_treatment_arms
  result$fully_parallel_cores <- resource_plan$fully_parallel_cores
  result$target_anchor_nlambda <- RoCE:::LAMBDA_GRID_SIZE_STANDARD
  result$task_elapsed_seconds <- as.numeric(task_elapsed_seconds)
  result$package_library <- package_provenance$library
  result$package_fingerprint <- package_provenance$fingerprint
  result$workflow_fingerprint <- workflow_fingerprint
  result$manifest_fingerprint <- manifest_fingerprint
  scheduler <- roce_scheduler_provenance()
  for (field in names(scheduler)) {
    result[[field]] <- scheduler[[field]]
  }
  result
}

roce_commit_csv_bundle <- function(
    values, output_paths, sentinel_path = NULL, sentinel_lines = NULL) {
  output_paths <- as.character(output_paths)
  if (!is.list(values) || length(values) == 0L ||
      length(values) != length(output_paths) ||
      any(!vapply(values, is.data.frame, logical(1L)))) {
    stop("CSV bundle values must be a nonempty list of data frames.",
         call. = FALSE)
  }
  if (anyNA(output_paths) || any(!nzchar(trimws(output_paths))) ||
      anyDuplicated(output_paths)) {
    stop("CSV bundle values and unique output paths must align.", call. = FALSE)
  }
  if (is.null(sentinel_path) != is.null(sentinel_lines)) {
    stop("sentinel_path and sentinel_lines must be supplied together.",
         call. = FALSE)
  }
  if (!is.null(sentinel_path)) {
    sentinel_path <- as.character(sentinel_path)
    sentinel_lines <- as.character(sentinel_lines)
    if (length(sentinel_path) != 1L || is.na(sentinel_path) ||
        !nzchar(trimws(sentinel_path)) || sentinel_path %in% output_paths) {
      stop("sentinel path must be nonempty and distinct from CSV outputs.",
           call. = FALSE)
    }
    if (length(sentinel_lines) == 0L || anyNA(sentinel_lines) ||
        any(!nzchar(trimws(sentinel_lines)))) {
      stop("sentinel lines must be nonempty and nonmissing.",
           call. = FALSE)
    }
  }
  final_paths <- c(output_paths, sentinel_path)
  if (any(file.exists(final_paths))) {
    stop(
      "refusing to overwrite existing bundle path(s): ",
      paste(final_paths[file.exists(final_paths)], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(lapply(unique(dirname(final_paths)), dir.create,
                   recursive = TRUE, showWarnings = FALSE))
  staging_root <- tempfile(
    pattern = ".roce_csv_bundle_", tmpdir = dirname(output_paths[[1L]])
  )
  if (!dir.create(staging_root)) {
    stop("failed to create CSV bundle staging directory.", call. = FALSE)
  }
  staged_paths <- file.path(
    staging_root, sprintf("value_%04d.csv", seq_along(values))
  )
  committed_paths <- character(0)
  completed <- FALSE
  on.exit({
    if (!completed && length(committed_paths) > 0L) {
      unlink(committed_paths, force = TRUE)
    }
    unlink(staging_root, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  for (index in seq_along(values)) {
    utils::write.csv(
      values[[index]], staged_paths[[index]], row.names = FALSE
    )
  }
  staged_sentinel <- NULL
  if (!is.null(sentinel_path)) {
    staged_sentinel <- file.path(staging_root, "commit_sentinel.tmp")
    writeLines(as.character(sentinel_lines), staged_sentinel)
  }
  for (index in seq_along(staged_paths)) {
    if (!file.rename(staged_paths[[index]], output_paths[[index]])) {
      stop("failed to commit CSV bundle output: ", output_paths[[index]],
           call. = FALSE)
    }
    committed_paths <- c(committed_paths, output_paths[[index]])
  }
  if (!is.null(sentinel_path)) {
    if (!file.rename(staged_sentinel, sentinel_path)) {
      stop("failed to commit CSV bundle sentinel: ", sentinel_path,
           call. = FALSE)
    }
    committed_paths <- c(committed_paths, sentinel_path)
  }
  completed <- TRUE
  invisible(final_paths)
}
