roce_read_positive_integer_env <- function(name, default, maximum = Inf) {
  raw_value <- trimws(Sys.getenv(name, as.character(default)))
  value <- suppressWarnings(as.integer(raw_value))
  if (length(value) != 1L || is.na(value) || value < 1L ||
      !identical(as.character(value), raw_value) || value > maximum) {
    maximum_text <- if (is.finite(maximum)) {
      paste0(" and at most ", as.integer(maximum))
    } else {
      ""
    }
    stop(
      name, " must be one positive integer", maximum_text, ".",
      call. = FALSE
    )
  }
  value
}

roce_slurm_resource_plan <- function(
    source_count, n_folds, allocated_cores, nuisance_cv_threads) {
  values <- c(
    source_count = source_count,
    n_folds = n_folds,
    allocated_cores = allocated_cores,
    nuisance_cv_threads = nuisance_cv_threads
  )
  integer_values <- suppressWarnings(as.integer(values))
  names(integer_values) <- names(values)
  if (anyNA(integer_values) || any(integer_values < 1L) ||
      !isTRUE(all.equal(
        as.numeric(values), as.numeric(integer_values), tolerance = 0
      ))) {
    stop("resource-plan inputs must be positive integers.", call. = FALSE)
  }
  values <- integer_values
  if (values[["nuisance_cv_threads"]] > values[["n_folds"]]) {
    stop(
      "ROCE_NUISANCE_CV_THREADS cannot exceed the nuisance CV fold count.",
      call. = FALSE
    )
  }
  if (values[["allocated_cores"]] < values[["nuisance_cv_threads"]]) {
    stop(
      "SLURM_CPUS_PER_TASK is smaller than ROCE_NUISANCE_CV_THREADS.",
      call. = FALSE
    )
  }

  fully_parallel_cores <-
    2L * values[["source_count"]] * values[["nuisance_cv_threads"]]
  parallel_treatment_arms <-
    values[["allocated_cores"]] >= fully_parallel_cores
  arm_divisor <- if (parallel_treatment_arms) 2L else 1L
  source_workers <- min(
    values[["source_count"]],
    floor(
      values[["allocated_cores"]] /
        (arm_divisor * values[["nuisance_cv_threads"]])
    )
  )
  if (source_workers < 1L) {
    stop("the resolved resource plan has no source worker.", call. = FALSE)
  }

  list(
    allocated_cores = values[["allocated_cores"]],
    nuisance_cv_threads = values[["nuisance_cv_threads"]],
    source_workers = as.integer(source_workers),
    parallel_treatment_arms = parallel_treatment_arms,
    fully_parallel_cores = as.integer(fully_parallel_cores)
  )
}

roce_resource_metadata_columns <- function() {
  c(
    "allocated_cores", "nuisance_cv_threads", "source_workers_per_arm",
    "parallel_treatment_arms", "fully_parallel_cores"
  )
}

roce_as_strict_logical <- function(values) {
  if (is.logical(values)) {
    return(values)
  }
  normalized <- toupper(trimws(as.character(values)))
  result <- rep(NA, length(normalized))
  result[normalized %in% c("TRUE", "T", "1")] <- TRUE
  result[normalized %in% c("FALSE", "F", "0")] <- FALSE
  result
}

roce_validate_result_resource_metadata <- function(
    result, expected_threads = NULL) {
  if (!is.data.frame(result) || nrow(result) < 1L) {
    stop("result must be a nonempty data frame.", call. = FALSE)
  }
  required_columns <- c(
    "K", "n_folds", roce_resource_metadata_columns()
  )
  missing_columns <- setdiff(required_columns, names(result))
  if (length(missing_columns) > 0L) {
    stop(
      "result is missing resource metadata: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!is.null(expected_threads)) {
    expected_threads <- suppressWarnings(as.integer(expected_threads))
    if (length(expected_threads) != 1L || is.na(expected_threads) ||
        expected_threads < 1L) {
      stop("expected_threads must be one positive integer.", call. = FALSE)
    }
  }

  integer_columns <- c(
    "K", "n_folds", "allocated_cores", "nuisance_cv_threads",
    "source_workers_per_arm", "fully_parallel_cores"
  )
  normalized <- result[integer_columns]
  normalized[] <- lapply(normalized, function(values) {
    numeric_values <- suppressWarnings(as.numeric(values))
    integer_values <- suppressWarnings(as.integer(numeric_values))
    if (anyNA(numeric_values) || any(!is.finite(numeric_values)) ||
        any(numeric_values != integer_values)) {
      return(rep(NA_integer_, length(values)))
    }
    integer_values
  })
  normalized$parallel_treatment_arms <-
    roce_as_strict_logical(result$parallel_treatment_arms)
  if (anyNA(normalized)) {
    stop("resource metadata contains an invalid value.", call. = FALSE)
  }
  if (!is.null(expected_threads) &&
      any(normalized$nuisance_cv_threads != expected_threads)) {
    stop(
      "resource metadata does not use the expected nuisance-CV thread count.",
      call. = FALSE
    )
  }

  combinations <- unique(normalized)
  for (index in seq_len(nrow(combinations))) {
    values <- combinations[index, , drop = FALSE]
    plan <- roce_slurm_resource_plan(
      source_count = values$K,
      n_folds = values$n_folds,
      allocated_cores = values$allocated_cores,
      nuisance_cv_threads = values$nuisance_cv_threads
    )
    observed <- list(
      source_workers = values$source_workers_per_arm,
      parallel_treatment_arms = values$parallel_treatment_arms,
      fully_parallel_cores = values$fully_parallel_cores
    )
    expected <- plan[c(
      "source_workers", "parallel_treatment_arms", "fully_parallel_cores"
    )]
    if (!identical(observed, expected)) {
      stop(
        "resource metadata is inconsistent with its Slurm resource plan.",
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

roce_expected_production_cv_threads <- function(source_count) {
  source_count <- suppressWarnings(as.integer(source_count))
  if (anyNA(source_count) || any(!source_count %in% c(2L, 4L, 8L))) {
    stop("production source count must be 2, 4, or 8.", call. = FALSE)
  }
  ifelse(source_count == 8L, 2L, 5L)
}

roce_validate_production_resource_metadata <- function(result) {
  roce_validate_result_resource_metadata(result)
  expected <- roce_expected_production_cv_threads(result$K)
  observed <- suppressWarnings(as.integer(result$nuisance_cv_threads))
  if (anyNA(observed) || !identical(observed, as.integer(expected))) {
    stop(
      paste0(
        "production resource metadata must use five nuisance-CV threads ",
        "for K=2/4 and two for K=8."
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}
