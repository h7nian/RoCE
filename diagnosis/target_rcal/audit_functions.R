#!/usr/bin/env Rscript
# Static inventory, not a proof that an unreferenced function is dead. Exported
# APIs and dynamic calls require manual review before deletion or renaming.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: audit_functions.R OUTPUT_DIRECTORY")
  output <- args[[1L]]
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
  namespace <- asNamespace("RoCE")
  exports <- getNamespaceExports("RoCE")
  source_files <- list.files("R", "[.]R$", full.names = TRUE)
  runner_files <- c("main.R", "realdata.R",
                    list.files("scripts", "[.]R$", recursive = TRUE, full.names = TRUE))
  test_files <- list.files("tests/testthat", "[.]R$", full.names = TRUE)
  functions <- list()
  parse_rows <- lapply(c(source_files, runner_files, test_files), function(path) {
    error <- ""
    expressions <- tryCatch(parse(path, keep.source = TRUE), error = function(condition) {
      error <<- conditionMessage(condition)
      expression()
    })
    if (path %in% source_files) for (index in seq_along(expressions)) {
      statement <- expressions[[index]]
      if (!is.call(statement) || !is.symbol(statement[[1L]]) ||
          !(as.character(statement[[1L]]) %in% c("<-", "=")) ||
          !is.symbol(statement[[2L]]) || !is.call(statement[[3L]]) ||
          !identical(statement[[3L]][[1L]], as.name("function"))) next
      name <- as.character(statement[[2L]])
      reference <- attr(expressions, "srcref")[[index]]
      functions[[name]] <<- list(name = name, path = path, line = as.integer(reference)[1L],
                                definition = statement[[3L]])
    }
    data.frame(path, parsed = !nzchar(error), error)
  })
  write.csv(do.call(rbind, parse_rows), file.path(output, "parse_checks.csv"), row.names = FALSE)
  test_text <- paste(unlist(lapply(test_files, readLines, warn = FALSE)), collapse = "\n")
  runner_text <- paste(unlist(lapply(runner_files, readLines, warn = FALSE)), collapse = "\n")
  rows <- lapply(functions, function(item) {
    fn <- get(item$name, envir = namespace, inherits = FALSE)
    parameter_names <- names(formals(fn))
    symbols <- all.names(body(fn), functions = TRUE, unique = TRUE)
    calls <- codetools::findGlobals(fn, merge = FALSE)$functions
    data.frame(name = item$name, path = item$path, line = item$line,
               exported = item$name %in% exports,
               n_parameters = length(parameter_names),
               parameters = paste(parameter_names, collapse = ";"),
               parameters_absent_from_body = paste(setdiff(parameter_names, c(symbols, "...")),
                                                   collapse = ";"),
               package_calls = paste(intersect(calls, names(functions)), collapse = ";"),
               mentioned_in_tests = grepl(item$name, test_text, fixed = TRUE),
               mentioned_in_runners = grepl(item$name, runner_text, fixed = TRUE))
  })
  inventory <- do.call(rbind, rows)
  rownames(inventory) <- NULL
  write.csv(inventory, file.path(output, "function_inventory.csv"), row.names = FALSE)
  usages <- unlist(lapply(functions, function(item) {
    fn <- get(item$name, envir = namespace, inherits = FALSE)
    capture.output(codetools::checkUsage(fn, name = item$name))
  }), use.names = FALSE)
  writeLines(usages, file.path(output, "codetools_findings.txt"))
  write.csv(inventory[nzchar(inventory$parameters_absent_from_body), ],
            file.path(output, "unused_parameter_candidates.csv"), row.names = FALSE)
  cat(sprintf("Parsed %d files; inventoried %d functions (%d exported).\n",
              length(parse_rows), nrow(inventory), sum(inventory$exported)))
  cat(sprintf("Unused-parameter candidates: %d functions. Codetools messages: %d.\n",
              sum(nzchar(inventory$parameters_absent_from_body)), length(usages)))
  stopifnot(all(vapply(parse_rows, function(row) row$parsed, logical(1L))))
}
if (sys.nframe() == 0L) main()
