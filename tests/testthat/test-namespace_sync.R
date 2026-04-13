# test-namespace_sync.R - Guard against drift between roxygen exports and NAMESPACE

library(testthat)

# Package loaded by helper-load.R (all functions available via FACEC namespace)

extract_roxygen_exports <- function(r_files) {
  exported <- character(0)

  for (file in r_files) {
    lines <- readLines(file, warn = FALSE)
    if (!length(lines)) next

    export_idx <- grep("^\\s*#'\\s*@export\\b", lines)
    if (!length(export_idx)) next

    for (idx in export_idx) {
      # Find the first non-empty, non-roxygen line after @export
      j <- idx + 1L
      while (j <= length(lines) && (grepl("^\\s*$", lines[j]) || grepl("^\\s*#'", lines[j]))) {
        j <- j + 1L
      }
      if (j > length(lines)) next

      line <- lines[j]
      # Match standard assignment forms: name <- function(...), name <- value
      m <- regexec("^\\s*([A-Za-z][A-Za-z0-9._]*)\\s*<-", line)
      g <- regmatches(line, m)[[1]]
      if (length(g) >= 2L) {
        exported <- c(exported, g[2])
      }
    }
  }

  unique(exported)
}

extract_namespace_exports <- function(namespace_file) {
  lines <- readLines(namespace_file, warn = FALSE)
  export_lines <- grep("^\\s*export\\(", lines, value = TRUE)
  if (!length(export_lines)) return(character(0))

  # Handles one export(name) per line (current project style)
  sub("^\\s*export\\(([^)]+)\\)\\s*$", "\\1", export_lines)
}

test_that("roxygen @export tags and NAMESPACE exports stay synchronized", {
  root <- normalizePath(file.path("..", ".."), mustWork = FALSE)
  ns_file <- file.path(root, "NAMESPACE")
  r_dir <- file.path(root, "R")

  skip_if_not(file.exists(ns_file), "NAMESPACE file not found")
  skip_if_not(dir.exists(r_dir), "R directory not found")

  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  roxy_exports <- sort(extract_roxygen_exports(r_files))
  ns_exports <- sort(extract_namespace_exports(ns_file))

  missing_in_ns <- setdiff(roxy_exports, ns_exports)
  extra_in_ns <- setdiff(ns_exports, roxy_exports)

  expect(
    length(missing_in_ns) == 0,
    paste0(
      "Roxygen @export symbols missing from NAMESPACE: ",
      paste(missing_in_ns, collapse = ", ")
    )
  )

  expect(
    length(extra_in_ns) == 0,
    paste0(
      "NAMESPACE exports without matching roxygen @export: ",
      paste(extra_in_ns, collapse = ", ")
    )
  )
})
