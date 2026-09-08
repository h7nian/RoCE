#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/bootstrap_benchmark"
}
benchmark_repeats <- if (length(args) >= 2L) {
  as.integer(args[[2L]])
} else {
  10L
}
if (is.na(benchmark_repeats) || benchmark_repeats < 2L) {
  stop("benchmark_repeats must be an integer >= 2.", call. = FALSE)
}

project_library <- Sys.getenv(
  "ROCE_PROJECT_LIB",
  "results/direct_tate_mc500_b5000/Rlib_current"
)
if (!dir.exists(project_library)) {
  stop("tested RoCE library not found: ", project_library, call. = FALSE)
}
.libPaths(c(normalizePath(project_library), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
package_provenance <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment(
  "ROCE_BOOTSTRAP_WORKFLOW_FINGERPRINT"
)
output_paths <- file.path(
  output_root,
  c("bootstrap_benchmark_raw.csv", "bootstrap_benchmark_summary.csv")
)
if (any(file.exists(output_paths))) {
  stop(
    "refusing to overwrite existing bootstrap benchmark output(s): ",
    paste(output_paths[file.exists(output_paths)], collapse = ", "),
    call. = FALSE
  )
}

# Match the primary K=4 design: one target and four source influence blocks,
# each with 1,000 observations. Fixed blocks isolate multiplier Monte Carlo
# error and runtime from nuisance-fitting variability.
set.seed(20260812L)
block_sizes <- rep(1000L, 5L)
block_weights <- c(0.36, 0.22, 0.18, 0.14, 0.10)
blocks <- Map(function(block_size, block_weight, block_index) {
  influence <- stats::rnorm(
    block_size,
    mean = 0,
    sd = 0.8 + 0.2 * block_index
  )
  list(influence = influence, weight = block_weight)
}, block_sizes, block_weights, seq_along(block_sizes))

analytic_se <- sqrt(sum(vapply(blocks, function(block) {
  centered <- block$influence - mean(block$influence)
  block$weight^2 * mean(centered^2) / length(centered)
}, numeric(1L))))

bootstrap_counts <- c(1000L, 2000L, 5000L)
records <- vector("list", length(bootstrap_counts) * benchmark_repeats)
record_index <- 0L
for (n_bootstrap in bootstrap_counts) {
  for (replication in seq_len(benchmark_repeats)) {
    set.seed(910000L + 100L * n_bootstrap + replication)
    elapsed <- system.time({
      bootstrap_se <- RoCE:::.multiplier_bootstrap_se(
        blocks, n_bootstrap = n_bootstrap
      )
    })[["elapsed"]]
    record_index <- record_index + 1L
    records[[record_index]] <- data.frame(
      n_bootstrap = n_bootstrap,
      benchmark_replication = replication,
      analytic_se = analytic_se,
      bootstrap_se = bootstrap_se,
      bootstrap_to_analytic = bootstrap_se / analytic_se,
      absolute_relative_error = abs(bootstrap_se / analytic_se - 1),
      elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  }
}
raw <- do.call(rbind, records)
summary <- do.call(rbind, lapply(split(raw, raw$n_bootstrap), function(group) {
  data.frame(
    n_bootstrap = group$n_bootstrap[[1L]],
    benchmark_repeats = nrow(group),
    mean_bootstrap_se = mean(group$bootstrap_se),
    sd_bootstrap_se = stats::sd(group$bootstrap_se),
    relative_sd_bootstrap_se = stats::sd(group$bootstrap_se) /
      mean(group$bootstrap_se),
    mean_absolute_relative_error = mean(group$absolute_relative_error),
    max_absolute_relative_error = max(group$absolute_relative_error),
    mean_elapsed_seconds = mean(group$elapsed_seconds),
    max_elapsed_seconds = max(group$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}))
rownames(summary) <- NULL
raw$package_library <- package_provenance$library
raw$package_fingerprint <- package_provenance$fingerprint
raw$workflow_fingerprint <- workflow_fingerprint
summary$package_library <- package_provenance$library
summary$package_fingerprint <- package_provenance$fingerprint
summary$workflow_fingerprint <- workflow_fingerprint

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
roce_commit_csv_bundle(
  values = list(raw, summary),
  output_paths = output_paths
)
print(summary, row.names = FALSE)
