#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) %in% c(2L, 3L))
root <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(startsWith(output, paste0(root, "/reviews/")), !dir.exists(output))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
completion <- do.call(rbind, lapply(split(manifest, manifest$config), function(group) {
  data.frame(scenario = group$config[1L], prescribed = nrow(group),
    complete = sum(vapply(group$task_id, function(id)
      file.exists(file.path(root, "tasks", id, "COMPLETE")), logical(1L))),
    failed = sum(vapply(group$task_id, function(id)
      file.exists(file.path(root, "tasks", id, "analysis_FAILED.txt")), logical(1L))))
}))
scope <- "All prescribed repeats"
if (length(arguments) == 3L) {
  if (arguments[3L] == "--complete-scenarios") {
    completed_scenarios <- completion$scenario[completion$complete == completion$prescribed]
    manifest <- manifest[manifest$config %in% completed_scenarios, , drop = FALSE]
    stopifnot(nrow(manifest) > 0L)
    scope <- paste("Interim review of complete scenarios only. Omitted scenarios:",
      paste(setdiff(completion$scenario, completed_scenarios), collapse = ", "))
  } else {
    count <- suppressWarnings(as.integer(arguments[3L]))
    stopifnot(length(count) == 1L, !is.na(count), count > 0L, count <= nrow(manifest),
      as.character(count) == arguments[3L])
    manifest <- head(manifest, count)
    scope <- paste("Interface smoke review of the first", count, "tasks; not a performance study")
  }
}
records <- lapply(manifest$task_id, function(id) {
  directory <- file.path(root, "tasks", id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")),
    file.exists(file.path(directory, "PAIRED_TARGET_FITS_PASSED")))
  rows <- read.csv(file.path(directory, "methods.csv"), stringsAsFactors = FALSE)
  task <- manifest[manifest$task_id == id, ]
  stopifnot(nrow(rows) == 5L, all(is.finite(rows$estimate)), all(is.finite(rows$se)), all(rows$se > 0),
    all(rows$scenario == task$config), all(rows$repeat_id == task$sim_id), length(unique(rows$truth)) == 1L,
    identical(sort(rows$source_radius[rows$method == "RoCE"]), c(2L, 3L, 5L)),
    all(rows$covered == (rows$ci_lower <= rows$truth & rows$ci_upper >= rows$truth)))
  for (radius in c(2, 3, 5)) {
    fits <- read.csv(file.path(directory, paste0("radius", radius), "nuisance_fits.csv"))
    failure_fields <- grep("nonconverged|line_search_failures", names(fits), value = TRUE)
    stopifnot(length(failure_fields) == 6L, all(is.finite(as.matrix(fits[failure_fields]))),
      all(as.matrix(fits[failure_fields]) == 0))
  }
  rows
})
rows <- do.call(rbind, records)
stopifnot(!anyDuplicated(rows[c("scenario", "repeat_id", "method", "source_radius")]))
rows$label <- ifelse(rows$method == "RoCE", paste0("RoCE-r", rows$source_radius), rows$method)
summaries <- lapply(split(rows, list(rows$scenario, rows$label), drop = TRUE), function(group) {
  error <- group$estimate - group$truth; n <- nrow(group)
  coverage <- mean(group$covered); z <- qnorm(.975)
  center <- (coverage + z^2 / (2 * n)) / (1 + z^2 / n)
  half_width <- z * sqrt(coverage * (1 - coverage) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  data.frame(scenario = group$scenario[1L], method = group$label[1L], repeats = n,
    truth = group$truth[1L], bias = mean(error), bias_mcse = sd(error) / sqrt(n),
    rmse = sqrt(mean(error^2)), empirical_sd = sd(group$estimate), mean_se = mean(group$se),
    coverage = coverage, coverage_lower = center - half_width, coverage_upper = center + half_width)
})
paired <- lapply(split(rows, rows$scenario), function(group) {
  reference <- group[group$method == "Target-only", ]
  do.call(rbind, lapply(c(2, 3, 5), function(radius) {
    candidate <- group[group$method == "RoCE" & !is.na(group$source_radius) & group$source_radius == radius, ]
    target <- reference[match(candidate$repeat_id, reference$repeat_id), ]
    difference <- (candidate$estimate - candidate$truth)^2 - (target$estimate - target$truth)^2
    data.frame(scenario = group$scenario[1L], radius = radius, repeats = nrow(candidate),
      paired_mse_difference = mean(difference), paired_mse_difference_mcse = sd(difference) / sqrt(length(difference)))
  }))
})
dir.create(output, recursive = TRUE)
writeLines(scope, file.path(output, "scope.txt"))
write.csv(completion, file.path(output, "scenario_completion.csv"), row.names = FALSE)
write.csv(rows, file.path(output, "replicates.csv"), row.names = FALSE)
write.csv(do.call(rbind, summaries), file.path(output, "summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, paired), file.path(output, "paired_mse.csv"), row.names = FALSE)
writeLines(sprintf("Checked all %d prescribed tasks in this review; paired inputs and convergence passed.", nrow(manifest)),
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, summaries))
