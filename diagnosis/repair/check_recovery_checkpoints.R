#!/usr/bin/env Rscript
# Read every finalized checkpoint before recovering an interrupted campaign.
# Estimator cache identity/model validation still runs when a worker resumes.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
inventory <- read.csv(arguments[1L], stringsAsFactors = FALSE)
output <- arguments[2L]
prefix <- "/scratch.global/zhan9381/FACE-HD/"
stopifnot(startsWith(output, prefix), !file.exists(output), nrow(inventory) > 0L,
          all(startsWith(inventory$path, prefix)), !anyDuplicated(inventory$path))
groups <- split(seq_len(nrow(inventory)), rep(seq_len(4L), length.out = nrow(inventory)))
results <- parallel::mclapply(groups, function(indices) {
  lapply(indices, function(index) {
    path <- inventory$path[index]
    error <- tryCatch({ readRDS(path); "" }, error = function(error) conditionMessage(error))
    data.frame(path = path, readable = !nzchar(error), error = error, stringsAsFactors = FALSE)
  })
}, mc.cores = 4L)
stopifnot(!any(vapply(results, inherits, logical(1L), "try-error")))
result <- do.call(rbind, unlist(results, recursive = FALSE))
stopifnot(nrow(result) == nrow(inventory), setequal(result$path, inventory$path))
write.csv(result, output, row.names = FALSE)
cat("Read", nrow(result), "checkpoints; unreadable:", sum(!result$readable), "\n")
if (any(!result$readable)) print(result[!result$readable, ], row.names = FALSE)
