#!/usr/bin/env Rscript
root <- Sys.getenv("ROCE_CHECK_ROOT")
stopifnot(nzchar(root), startsWith(normalizePath(root), "/scratch.global/zhan9381/FACE-HD/"))
.libPaths(c(file.path(root, "Rlib"), .libPaths()))
stopifnot(identical(normalizePath(find.package("RoCE")),
                    normalizePath(file.path(root, "Rlib", "RoCE"))))
filter <- Sys.getenv("ROCE_TEST_FILTER")
output <- Sys.getenv("ROCE_TEST_OUTPUT", root)
stopifnot(startsWith(normalizePath(output), "/scratch.global/zhan9381/FACE-HD/"))
result <- testthat::test_dir(
  "tests/testthat", filter = if (nzchar(filter)) filter else NULL,
  reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE
)
saveRDS(result, file.path(output, "test_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
summary <- as.data.frame(result)
write.csv(summary[setdiff(names(summary), "result")], file.path(output, "test_summary.csv"), row.names = FALSE)
if (any(summary$failed > 0 | summary$error)) stop("Test failures; structured results retained.")
