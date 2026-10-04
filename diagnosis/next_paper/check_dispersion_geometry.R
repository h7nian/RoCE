#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !startsWith(arguments[[1L]], "/scratch.global/zhan9381/FACE-HD/")) {
  stop("Usage: new_FACE_HD_scratch_output_directory")
}
output <- arguments[[1L]]
if (dir.exists(output) || !dir.create(output, recursive = TRUE)) stop("Use a new output directory")
set.seed(68123)
records <- list()
for (count in 2:8) for (trial in 1:10) {
  loading <- matrix(rnorm(count * count), count)
  covariance <- tcrossprod(loading) + diag(.5, count)
  inverse <- solve(covariance)
  precision <- drop(inverse %*% rep(1, count))
  variance <- 1 / sum(precision)
  weights <- variance * precision
  projection <- inverse - variance * tcrossprod(precision)
  for (valid_count in seq_len(count - 1L)) {
    valid <- sort(sample(count, valid_count))
    invalid <- setdiff(seq_len(count), valid)
    subset_variance <- 1 / sum(solve(covariance[valid, valid, drop = FALSE]))
    direct <- drop(crossprod(weights[invalid],
      solve(projection[invalid, invalid, drop = FALSE], weights[invalid])))
    identity_error <- abs(direct - (subset_variance - variance))
    upper <- max(diag(covariance)) / valid_count + (valid_count - 1) / valid_count *
      max(covariance[row(covariance) != col(covariance)]) - variance
    bias <- numeric(count)
    bias[invalid] <- rnorm(count - valid_count)
    bias_slack <- direct * drop(crossprod(bias, projection %*% bias)) - sum(weights * bias)^2
    stopifnot(identity_error < 1e-9, upper + 1e-10 >= direct, bias_slack >= -1e-10)
    records[[length(records) + 1L]] <- data.frame(count = count, trial = trial, valid_count = valid_count,
      identity_error = identity_error, covariance_bound_slack = upper - direct, bias_bound_slack = bias_slack)
  }
}
records <- do.call(rbind, records)
write.csv(records, file.path(output, "geometry_checks.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
cat("DISPERSION_GEOMETRY_PASSED", nrow(records), "subset cases; maximum identity error",
    max(records$identity_error), "\n")
