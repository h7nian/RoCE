#!/usr/bin/env Rscript
# One complete p100 inner calibration task using already saved 1000/site data.
# The output is a computational/estimator-component check, not MC coverage.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("usage: check_calibration_cost.R LIBRARY INPUT_RDS NEW_OUTPUT [legacy|score_derivative]")
}
recipe <- if (length(args) == 4L) match.arg(args[4L], c("legacy", "score_derivative")) else "legacy"
control <- if (recipe == "score_derivative") {
  list(recipe = recipe, target_propensity_initialization = "calibrated", target_radius = log(9))
} else list(recipe = "legacy")
target_radius <- if (recipe == "score_derivative") log(9) else 5
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
saved <- readRDS(args[2L])
output <- args[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
RoCE:::set_nuisance_solver_cpp("proximal_newton")
folds <- saved$folds
cache <- new.env(parent = emptyenv())
writeLines(c("1000/site, raw p100, 200 working features, ten original folds",
  "outer fold 1 and validation fold 2 excluded; calibration folds 3:10",
  paste("one source and target, treated arm only; calibration recipe", recipe),
  "100 lambdas, nuisance tolerance 1e-10, compact matrices, one CPU"), file.path(output, "scope.txt"))
started <- proc.time()[["elapsed"]]
messages <- RoCE:::.prepare_source_calibration_messages(
  folds$target_folds, folds$source_folds$s1, "s1", 3:10, "one_round",
  1L, "binomial", 5, 100L, "min", cache, calibration_recipe = recipe)
message_seconds <- proc.time()[["elapsed"]] - started
saveRDS(messages, file.path(output, "target_messages.rds"))
cat("Target messages ready after", message_seconds, "seconds\n")
started <- proc.time()[["elapsed"]]
source_fit <- RoCE:::.fit_source_calibration(
  folds$source_folds$s1, "s1", 3:10, function(site, fold) messages[[paste0("k2_", fold)]],
  1L, 1L, 1L, 5, 100L, 10000L, "min", cache, layout = "compact", tol = 1e-10,
  calibration_recipe = recipe)
source_seconds <- proc.time()[["elapsed"]] - started
saveRDS(source_fit, file.path(output, "source_fit.rds"))
cat("Source calibration finished after", source_seconds, "seconds\n")
started <- proc.time()[["elapsed"]]
target_fit <- RoCE:::.fit_target_calibration(folds$target_folds, 3:10,
  1L, "binomial", 5, 100L, lambda_rule = "min", fit_cache = cache,
  layout = "compact", tol = 1e-10, calibration_control = control)
target_seconds <- proc.time()[["elapsed"]] - started
saveRDS(target_fit, file.path(output, "target_fit.rds"))
target_fold <- RoCE:::materialize_fold(folds$target_folds, 2L)
source_fold <- RoCE:::materialize_fold(folds$source_folds$s1, 2L)
source_evaluation <- RoCE:::.compute_inner_source_arm_info(source_fold, target_fold,
  source_fit$weight, source_fit$outcome, 5, 1L, 1L, 1L)
target_evaluation <- RoCE:::.evaluate_target_calibration(target_fit, target_fold, 1L, "binomial", target_radius)
write.csv(data.frame(component = c("source_assisted_mu1", "target_anchor_mu1"),
  estimate = c(source_evaluation$estimate, target_evaluation$estimate)),
  file.path(output, "held_out_estimates.csv"), row.names = FALSE)
timing <- rbind(
  data.frame(component = "target_messages", seconds = message_seconds),
  data.frame(component = "source_total", seconds = source_seconds),
  data.frame(component = "target_total", seconds = target_seconds),
  data.frame(component = paste0("source_", names(source_fit$timing)), seconds = as.numeric(source_fit$timing)),
  data.frame(component = paste0("target_", names(target_fit$timing)), seconds = as.numeric(target_fit$timing)))
write.csv(timing, file.path(output, "timing.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("PASSED", file.path(output, "status.txt"))
