roce_sha256_environment <- function(name) {
  value <- tolower(trimws(Sys.getenv(name, "")))
  if (!grepl("^[0-9a-f]{64}$", value)) {
    stop(name, " must be a 64-character SHA-256 digest.", call. = FALSE)
  }
  value
}

roce_sha256_file <- function(path) {
  if (length(path) != 1L || is.na(path) || !nzchar(trimws(path)) ||
      !file.exists(path) || dir.exists(path)) {
    stop("path must identify one existing regular file.", call. = FALSE)
  }
  normalized_path <- normalizePath(path, mustWork = TRUE)
  output <- system2(
    "sha256sum", args = shQuote(normalized_path),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  digest <- if (length(output) == 1L && status == 0L) {
    strsplit(trimws(output[[1L]]), "[[:space:]]+")[[1L]][[1L]]
  } else {
    ""
  }
  digest <- tolower(digest)
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    stop("could not compute a valid SHA-256 file fingerprint.",
         call. = FALSE)
  }
  digest
}

roce_result_provenance <- function(results, expected_library = "") {
  required_columns <- c("package_library", "package_fingerprint")
  missing_columns <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    return(list(
      passed = FALSE,
      library = NA_character_,
      fingerprint = NA_character_,
      detail = paste0("missing=", paste(missing_columns, collapse = ","))
    ))
  }

  libraries <- unique(trimws(as.character(results$package_library)))
  fingerprints <- unique(tolower(trimws(
    as.character(results$package_fingerprint)
  )))
  libraries <- libraries[!is.na(libraries) & nzchar(libraries)]
  fingerprints <- fingerprints[
    !is.na(fingerprints) & nzchar(fingerprints)
  ]

  observed_library <- if (length(libraries) == 1L) libraries[[1L]] else NA_character_
  observed_fingerprint <- if (length(fingerprints) == 1L) {
    fingerprints[[1L]]
  } else {
    NA_character_
  }
  library_matches <- length(libraries) == 1L
  if (nzchar(expected_library) && library_matches) {
    expected_normalized <- normalizePath(expected_library, mustWork = FALSE)
    observed_normalized <- normalizePath(observed_library, mustWork = FALSE)
    library_matches <- identical(observed_normalized, expected_normalized)
  }
  fingerprint_valid <- length(fingerprints) == 1L &&
    grepl("^[0-9a-f]{64}$", observed_fingerprint)
  expected_fingerprint <- tolower(trimws(
    Sys.getenv("ROCE_PACKAGE_FINGERPRINT", "")
  ))
  fingerprint_matches <- !nzchar(expected_fingerprint) ||
    identical(observed_fingerprint, expected_fingerprint)

  list(
    passed = library_matches && fingerprint_valid && fingerprint_matches,
    library = observed_library,
    fingerprint = observed_fingerprint,
    detail = paste0(
      "libraries=", length(libraries),
      ";fingerprints=", length(fingerprints),
      ";library=", ifelse(is.na(observed_library), "mixed_or_missing", observed_library),
      ";fingerprint=", ifelse(
        is.na(observed_fingerprint), "mixed_or_missing", observed_fingerprint
      )
    )
  )
}

roce_runtime_package_provenance <- function(
    project_library, package = "RoCE") {
  fingerprint <- roce_sha256_environment("ROCE_PACKAGE_FINGERPRINT")
  loaded_library <- dirname(normalizePath(
    find.package(package), mustWork = TRUE
  ))
  if (!nzchar(project_library)) {
    stop("ROCE_PROJECT_LIB must identify the tested package library.",
         call. = FALSE)
  }
  expected_library <- normalizePath(project_library, mustWork = TRUE)
  if (!identical(loaded_library, expected_library)) {
    stop(package, " was not loaded from ROCE_PROJECT_LIB.", call. = FALSE)
  }
  list(library = loaded_library, fingerprint = fingerprint)
}
