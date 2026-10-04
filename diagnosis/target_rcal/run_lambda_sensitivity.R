#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_lambda_sensitivity.R TASK_ID OUTPUT_ROOT")
task_id <- as.integer(args[1L])
stopifnot(!is.na(task_id), task_id >= 1L, task_id <= 120L)
configuration <- c("C1", "C2", "C3")[(task_id - 1L) %/% 40L + 1L]
seed <- (task_id - 1L) %% 40L + 1L
.libPaths(c(Sys.getenv("ROCE_PROJECT_LIB"), .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source("scripts/slurm/result_provenance.R")
output <- file.path(args[2L], "lambda_sensitivity", configuration,
                    sprintf("seed_%04d", seed))
if (dir.exists(output)) stop("Refusing to overwrite a diagnostic seed.")
dir.create(output, recursive = TRUE)
set.seed(seed)
data <- RoCE:::generate_simulation_data(
  n_total = 3000L, K = 2L, p = 100L, config = configuration,
  estimand_type = "superpopulation", outcome_type = "binary", dgp_type = "face",
  ate_deviation = 0, n_deviated_sites = 0L, warn_ignored = FALSE
)
split <- RoCE:::split_data_by_site(data)
stopifnot(all(vapply(split, function(site) site$n, integer(1L)) == 1000L))
folds <- RoCE:::build_crossfit_folds(split, 10L)$target_folds
rows <- lapply(c("min", "1se"), function(rule) {
  arm_values <- lapply(c(1L, 0L), function(A_val) {
    unlist(lapply(seq_len(10L), function(k1) {
      fit <- RoCE::estimate_target_only_from_complement(
        folds, k1, 10L, A_val = A_val, nuisance_lambda_rule = rule
      )
      fit$estimate + fit$varphi_ot
    }), use.names = FALSE)
  })
  contrast <- arm_values[[1L]] - arm_values[[2L]]
  estimate <- mean(contrast)
  data.frame(configuration, seed, rule, n_site = 1000L, p = 100L, n_folds = 10L,
             estimate, bias = estimate - (data$mu1_true - data$mu0_true),
             se = sqrt(mean((contrast - estimate)^2) / length(contrast)))
})
write.csv(do.call(rbind, rows), file.path(output, "results.csv"), row.names = FALSE)
writeLines(c(paste0("library=", getNamespaceInfo("RoCE", "path")),
             paste0("script_sha256=", roce_sha256_file("diagnosis/target_rcal/run_lambda_sensitivity.R"))),
           file.path(output, "metadata.txt"))
writeLines("completed", file.path(output, "COMPLETED"))
