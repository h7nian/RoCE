#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "usage: render_rhc_direct_tate_figure.R METHODS.csv OUTPUT.pdf [OUTPUT.pdf ...]",
    call. = FALSE
  )
}

input_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_paths <- args[-1L]
normalized_outputs <- file.path(
  normalizePath(dirname(output_paths), mustWork = FALSE),
  basename(output_paths)
)
if (anyDuplicated(normalized_outputs)) {
  stop("output paths must be unique.", call. = FALSE)
}
if (any(tolower(tools::file_ext(output_paths)) != "pdf")) {
  stop("every output path must end in .pdf.", call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))

methods <- read.csv(input_path, stringsAsFactors = FALSE)
required_columns <- c("method", "estimate", "se", "ci_lower", "ci_upper")
missing_columns <- setdiff(required_columns, names(methods))
if (length(missing_columns) > 0L) {
  stop(
    "RHC methods file is missing column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}
expected_methods <- c(
  "Target-only", "SS", "IVW", "Federated-DR", "Pooled-DR", "RoCE"
)
if (!identical(as.character(methods$method), expected_methods) ||
    any(!is.finite(unlist(methods[setdiff(required_columns, "method")]))) ||
    any(methods$se <= 0)) {
  stop("RHC methods file failed its method-order or finiteness check.",
       call. = FALSE)
}

allow_overwrite <- identical(
  tolower(trimws(Sys.getenv("ROCE_ALLOW_FIGURE_OVERWRITE", "false"))),
  "true"
)
existing <- output_paths[file.exists(output_paths)]
if (length(existing) > 0L && !allow_overwrite) {
  stop(
    "refusing to overwrite existing figure(s): ",
    paste(existing, collapse = ", "),
    call. = FALSE
  )
}

forest <- plot_forest_methods(
  methods,
  highlight = "RoCE",
  xlab = "Target average treatment effect (risk difference)"
)
forest <- forest + ggplot2::theme(
  axis.title = ggplot2::element_text(size = 21),
  axis.text = ggplot2::element_text(size = 20)
)

render_and_install <- function() {
  temporary_figure <- tempfile(fileext = ".pdf")
  on.exit(unlink(temporary_figure), add = TRUE)
  save_plot(forest, temporary_figure, width = 8.2, height = 4.2)

  for (output_path in output_paths) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    temporary_output <- tempfile(
      pattern = paste0(".", basename(output_path), "_"),
      tmpdir = dirname(output_path),
      fileext = ".pdf"
    )
    if (!file.copy(temporary_figure, temporary_output, overwrite = TRUE) ||
        !file.rename(temporary_output, output_path)) {
      unlink(temporary_output)
      stop("failed to atomically install ", output_path, call. = FALSE)
    }
    message("[done] wrote ", output_path)
  }
}
render_and_install()
