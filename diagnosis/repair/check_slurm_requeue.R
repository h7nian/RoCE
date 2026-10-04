#!/usr/bin/env Rscript
# Real Slurm integration check: claim, save a fit, requeue, restore, reject duplicates.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1L], mustWork = TRUE)
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
.libPaths(c(configuration$library, .libPaths()))
library(RoCE)
library(testthat)
output <- file.path(root, "tasks", args[2L])
job_id <- Sys.getenv("SLURM_JOB_ID")
restart_count <- as.integer(Sys.getenv("SLURM_RESTART_COUNT", "0"))
stopifnot(grepl("^[0-9]+$", job_id), restart_count %in% 0:1)
owner <- jsonlite::fromJSON(file.path(output, "owner.json"))
stopifnot(owner$job_id == job_id, owner$restart == restart_count)
set.seed(8011)
fit_args <- list(W_outcome = matrix(rnorm(800), 200, 4),
                 Y = rbinom(200, 1, .5), A = rep(1L, 200), A_val = 1L, nlambda = 4L)
cache <- RoCE:::.new_nuisance_cache(TRUE, checkpoint_dir = file.path(output, "checkpoints"))
fit_once <- function() RoCE:::.fit_nuisance_training_subset(
  "fit_initial_outcome", fit_args, "s1", 2:4, cache)
if (restart_count == 0L) {
  fit <- fit_once()
  saveRDS(list(job_id = job_id, restart_count = restart_count,
               coefficients = as.numeric(fit), lambda = attr(fit, "lambda_used")),
          file.path(output, "before_requeue.rds"))
  writeLines("CHECKPOINTED", file.path(output, "status.txt"))
  stopifnot(system2("scontrol", c("requeue", job_id)) == 0L)
  Sys.sleep(30)
  stop("The requested Slurm requeue did not stop the original attempt.")
}
read_cached_fit <- function() {
  testthat::local_mocked_bindings(fit_initial_outcome = function(...) stop("Unexpected cache miss"),
                                 .package = "RoCE")
  fit_once()
}
fit <- read_cached_fit()
before <- readRDS(file.path(output, "before_requeue.rds"))
stopifnot(before$job_id == job_id, before$restart_count == 0L,
          identical(before$coefficients, as.numeric(fit)), identical(before$lambda, attr(fit, "lambda_used")))
duplicate_log <- file.path(output, "duplicate_guard_expected.txt")
duplicate_status <- suppressWarnings(system2("python3", c("-B",
  shQuote(file.path(root, "workflow", "submit_repeat_pilot.py")), "validate-run",
  shQuote(root), shQuote(args[2L])), stdout = duplicate_log, stderr = duplicate_log))
stopifnot(duplicate_status != 0L, any(grepl("already claimed", readLines(duplicate_log))))
jsonlite::write_json(list(job_id = job_id, initial_restart_count = 0L,
  resumed_restart_count = restart_count, exact_fit_restored = TRUE,
  duplicate_attempt_rejected = TRUE), file.path(output, "probe_result.json"), pretty = TRUE)
writeLines("COMPLETE", file.path(output, "status.txt"))
writeLines("Real Slurm requeue and exact-fit restoration passed.", file.path(output, "COMPLETE"))
