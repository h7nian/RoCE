#!/usr/bin/env Rscript
# One p100 baseline pipeline fixture, not a coverage or RMSE experiment.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: check_baseline_pipeline.R LIBRARY INPUT_RDS NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
data <- readRDS(args[2L])$data
stopifnot(all(vapply(data, function(site) site$n == 1000L && ncol(site$W_outcome) == 200L,
                    logical(1L))))
writeLines(c("C3 fixture: 1000/site, raw p100, 200 working features, two sources",
  "Baseline TATE pipeline with ten target folds and fixed sample sizes",
  "Conditional derivative validation only; not coverage or general post-selection inference"),
  file.path(output, "scope.txt"))
invisible(RoCE:::set_nuisance_solver_cpp("proximal_newton"))
set.seed(1159)
started <- proc.time()[["elapsed"]]
result <- tryCatch(run_all_comparisons_tate(data,
  methods = c("target_only", "federated_dr", "pooled_dr"), family = "binomial",
  n_folds = 10L, n_cores = 1L, variance_method = "analytic"), error = identity)
saveRDS(result, file.path(output, "results.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines(as.character(proc.time()[["elapsed"]] - started), file.path(output, "elapsed_seconds.txt"))
if (inherits(result, "error")) {
  writeLines(conditionMessage(result), file.path(output, "FAILED.txt"))
  stop(conditionMessage(result))
}
rows <- do.call(rbind, lapply(names(result), function(method) {
  fit <- result[[method]]
  reconstructed <- fit$components$mu1_influence_se^2 + fit$components$mu0_influence_se^2 -
    2 * fit$components$cross_arm_covariance
  stopifnot(is.finite(fit$estimate), is.finite(fit$se), fit$se > 0,
            abs(reconstructed - fit$variance) < 1e-10)
  data.frame(method = method, estimate = fit$estimate, se = fit$se,
    inference_scope = fit$components$inference_scope,
    cross_arm_covariance = fit$components$cross_arm_covariance)
}))
write.csv(rows, file.path(output, "estimates.csv"), row.names = FALSE)
writeLines("PASSED", file.path(output, "status.txt"))
