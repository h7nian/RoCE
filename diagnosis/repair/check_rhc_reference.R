#!/usr/bin/env Rscript
# Release a new RHC tuning study only after its default reproduces saved results.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 3L)
study <- normalizePath(arguments[1L], mustWork = TRUE)
reference_directory <- normalizePath(arguments[2L], mustWork = TRUE)
output <- arguments[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"))
candidate_directory <- file.path(study, "tasks/1")
stopifnot(file.exists(file.path(candidate_directory, "COMPLETE")),
          file.exists(file.path(reference_directory, "COMPLETE")))
source(file.path(study, "workflow/rhc_stage_signatures.R"))
first <- readRDS(file.path(reference_directory, "checkpoints/analysis.rds"))$value
second <- readRDS(file.path(candidate_directory, "checkpoints/analysis.rds"))$value
stopifnot(identical(readLines(file.path(reference_directory, "data_sha256.txt")),
                    readLines(file.path(candidate_directory, "data_sha256.txt"))))
checks <- list()
compare <- function(a, b, label, tolerance = 1e-10) {
  if (is.list(a) || is.list(b)) {
    stopifnot(is.list(a), is.list(b), identical(names(a), names(b)), length(a) == length(b))
    for (i in seq_along(a)) compare(a[[i]], b[[i]], paste(label, i, sep = ":"), tolerance)
  } else if (is.numeric(a) && is.numeric(b)) {
    stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)))
    difference <- if (length(a)) max(abs(a-b)) else 0
    if (difference > tolerance) stop("Default RHC parity failed: ", label)
    checks[[length(checks)+1L]] <<- data.frame(component = label, max_difference = difference)
  } else stopifnot(identical(a, b))
}
for (method in as.character(first$methods$method)) {
  a <- first$methods[as.character(first$methods$method) == method, c("estimate", "se")]
  b <- second$methods[as.character(second$methods$method) == method, c("estimate", "se")]
  compare(as.numeric(a[1L, ]), as.numeric(b[1L, ]), method)
}
compare(rhc_stage_signatures(first$tate_fit), rhc_stage_signatures(second$tate_fit), "fitted_stages")
compare(as.numeric(first$tate_fit$fold_weights), as.numeric(second$tate_fit$fold_weights), "aggregation_weights")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, checks), file.path(output, "default_parity.csv"), row.names = FALSE)
writeLines("RHC_DEFAULT_REFERENCE_PARITY_PASSED", file.path(output, "DEFAULT_PARITY_PASSED"))
cat("Default RHC parity passed; maximum numerical difference:",
    max(vapply(checks, function(row) row$max_difference, numeric(1L))), "\n")
