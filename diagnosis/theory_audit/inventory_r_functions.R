#!/usr/bin/env Rscript
# Parse functions, formals and assignments without sourcing project code.
# Presence in this inventory is not a semantic review or evidence of dead code.

inventory_file <- function(path) {
  expressions <- parse(path, keep.source = TRUE)
  tokens <- getParseData(expressions)
  function_tokens <- tokens[tokens$token == "FUNCTION", , drop = FALSE]
  function_tokens <- function_tokens[order(function_tokens$line1, function_tokens$col1), , drop = FALSE]
  function_index <- 0L
  functions <- parameters <- assignments <- list()

  describe <- function(value) paste(deparse(value, width.cutoff = 120L), collapse = " ")
  head_name <- function(node) {
    if (is.call(node) && is.symbol(node[[1L]])) as.character(node[[1L]]) else ""
  }
  collect_assignments <- function(node) {
    if (!is.call(node) || identical(head_name(node), "function")) return(character())
    found <- if (head_name(node) %in% c("<-", "=", "<<-") && is.symbol(node[[2L]])) {
      as.character(node[[2L]])
    } else character()
    for (index in seq_along(node)[-1L]) {
      if (!identical(node[[index]], quote(expr = ))) {
        found <- c(found, collect_assignments(node[[index]]))
      }
    }
    unique(found)
  }
  visit <- function(node, parent = "", assigned_name = "") {
    if (is.expression(node) || is.pairlist(node)) {
      for (index in seq_along(node)) {
        if (!identical(node[[index]], quote(expr = ))) visit(node[[index]], parent)
      }
      return(invisible(NULL))
    }
    if (!is.call(node)) return(invisible(NULL))
    if (head_name(node) %in% c("<-", "=", "<<-") && length(node) == 3L) {
      visit(node[[2L]], parent)
      visit(node[[3L]], parent, describe(node[[2L]]))
      return(invisible(NULL))
    }
    if (identical(head_name(node), "function")) {
      function_index <<- function_index + 1L
      location <- function_tokens[function_index, , drop = FALSE]
      label <- if (nzchar(assigned_name)) assigned_name else "<anonymous>"
      identifier <- paste(path, location$line1, location$col1, sep = ":")
      formals <- node[[2L]]
      body <- node[[3L]]
      body_names <- all.names(body, functions = TRUE, unique = TRUE)
      functions[[length(functions) + 1L]] <<- data.frame(
        function_id = identifier, path = path, line = location$line1,
        column = location$col1, name = label, parent_function = parent,
        parameter_count = length(formals), semantic_review = "pending",
        stringsAsFactors = FALSE
      )
      for (index in seq_along(formals)) {
        name <- names(formals)[index]
        has_default <- !identical(formals[[index]], quote(expr = ))
        parameters[[length(parameters) + 1L]] <<- data.frame(
          function_id = identifier, parameter = name, has_default = has_default,
          default = if (has_default) describe(formals[[index]]) else "",
          mentioned_in_body = name %in% body_names,
          stringsAsFactors = FALSE
        )
      }
      for (name in collect_assignments(body)) {
        assignments[[length(assignments) + 1L]] <<- data.frame(
          function_id = identifier, assigned_symbol = name,
          also_a_parameter = name %in% names(formals), stringsAsFactors = FALSE
        )
      }
      # Defaults can themselves contain function definitions.
      visit(formals, identifier)
      visit(body, identifier)
      return(invisible(NULL))
    }
    for (index in seq_along(node)) {
      if (!identical(node[[index]], quote(expr = ))) visit(node[[index]], parent)
    }
    invisible(NULL)
  }
  visit(expressions)
  if (function_index != nrow(function_tokens)) {
    stop(sprintf("Function traversal mismatch in %s: AST=%d, tokens=%d",
                 path, function_index, nrow(function_tokens)))
  }
  list(functions = functions, parameters = parameters, assignments = assignments)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: inventory_r_functions.R FILE_MANIFEST OUTPUT_DIRECTORY")
  paths <- read.csv(args[[1L]], stringsAsFactors = FALSE)$path
  paths <- paths[grepl("\\.[Rr]$", paths) | basename(paths) == ".Rprofile"]
  output <- args[[2L]]
  parent <- normalizePath(dirname(output), mustWork = TRUE)
  if (!startsWith(parent, "/scratch.global/")) stop("Output parent must be in scratch.global")
  if (file.exists(output)) stop("Output already exists; choose a new directory")
  dir.create(output)
  records <- list(functions = list(), parameters = list(), assignments = list())
  parse_checks <- list()
  for (path in paths) {
    result <- tryCatch(inventory_file(path), error = identity)
    error <- if (inherits(result, "error")) conditionMessage(result) else ""
    parse_checks[[length(parse_checks) + 1L]] <- data.frame(
      path = path, parsed_and_indexed = !nzchar(error), error = error
    )
    if (!nzchar(error)) {
      for (category in names(records)) records[[category]] <- c(records[[category]], result[[category]])
    }
  }
  for (category in names(records)) {
    write.csv(do.call(rbind, records[[category]]), file.path(output, paste0(category, ".csv")),
              row.names = FALSE)
  }
  parse_checks <- do.call(rbind, parse_checks)
  write.csv(parse_checks, file.path(output, "parse_checks.csv"), row.names = FALSE)
  cat(sprintf("Indexed %d functions in %d R files; %d indexing errors.\n",
              length(records$functions), length(paths), sum(!parse_checks$parsed_and_indexed)))
  if (any(!parse_checks$parsed_and_indexed)) stop("Review parse_checks.csv; inventory is incomplete")
}

if (sys.nframe() == 0L) main()
