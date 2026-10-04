#!/usr/bin/env Rscript
# Independently solve the saved four-coordinate aggregation problems by
# enumerating active signs. No nuisance models or simulation draws are refit.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) stop("Usage: campaign_directory output_directory")
root <- arguments[[1L]]
output <- arguments[[2L]]
if (!all(startsWith(c(root, output), "/scratch.global/zhan9381/FACE-HD/"))) {
  stop("Campaign and output must be on FACE-HD scratch")
}
configuration <- jsonlite::read_json(file.path(root, "configuration.json"), simplifyVector = TRUE)
.libPaths(c(configuration$library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
dir.create(output, recursive = TRUE, showWarnings = FALSE)

enumerate_optimum <- function(moments, penalty, fixed_zero = integer()) {
  dimension <- length(penalty)
  if (dimension > 4L) stop("This exhaustive diagnostic is limited to four coordinates")
  free <- setdiff(seq_len(dimension), fixed_zero)
  patterns <- as.matrix(expand.grid(rep(list(-1:1), length(free))))
  candidates <- list()
  for (index in seq_len(nrow(patterns))) {
    signs <- numeric(dimension)
    signs[free] <- patterns[index, ]
    active <- which(signs != 0)
    weights <- numeric(dimension)
    if (length(active)) {
      weights[active] <- solve(2 * moments$Q[active, active, drop = FALSE],
        -2 * moments$l[active] - penalty[active] * signs[active] / moments$N_all)
      if (any(weights[active] * signs[active] <= 0)) next
    }
    gradient <- 2 * drop(moments$Q %*% weights + moments$l)
    inactive <- setdiff(free, active)
    if (any(abs(gradient[inactive]) > penalty[inactive] / moments$N_all + 1e-12)) next
    objective <- moments$N_all * (moments$constant + 2 * sum(moments$l * weights) +
      drop(crossprod(weights, moments$Q %*% weights))) + sum(penalty * abs(weights))
    candidates[[length(candidates) + 1L]] <- list(weights = weights, objective = objective)
  }
  if (!length(candidates)) stop("No feasible optimality pattern found")
  candidates[[which.min(vapply(candidates, `[[`, numeric(1L), "objective"))]]
}

manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
tasks <- manifest[manifest$sim_id == configuration$first_seed, ]
rows <- list()
for (task_index in seq_len(nrow(tasks))) {
  task <- tasks[task_index, ]
  directory <- file.path(root, "tasks", task$task_id)
  stopifnot(file.exists(file.path(directory, "COMPLETE")), task$K == 2L)
  result <- readRDS(file.path(directory, "score_derivative.rds"))
  fitted <- attr(result, "roce_simulation_artifacts")$direct_tate_results[[paste0(task$protocol, "_crossfit")]]
  joint <- result[result$method == paste0(task$protocol, "_crossfit_ate_joint_tate"), ]
  stopifnot(nrow(joint) == 1L, !is.null(fitted))
  for (fold in seq_len(configuration$n_folds)) {
    records <- RoCE:::.joint_inner_records(
      fitted$arm_results$mu1$intermediates$inner_fold_info[[fold]],
      fitted$arm_results$mu0$intermediates$inner_fold_info[[fold]])
    moments <- RoCE:::.joint_inner_moments(records)
    wald <- abs(moments$discrepancy) / sqrt(diag(moments$Q))
    penalty <- pmax(configuration$aggregation_lambda * wald - 1, 0)
    columns <- unlist(lapply(c("mu1", "mu0"), function(arm) {
      paste0(arm, "_source_s", 1:2, "_fold_", fold, "_weight")
    }))
    stored <- as.numeric(joint[1L, columns])
    exact <- enumerate_optimum(moments, penalty)
    restricted <- enumerate_optimum(moments, penalty, fixed_zero = 1L)
    difference <- max(abs(stored - exact$weights))
    stopifnot(difference < 1e-6)
    gradient <- 2 * moments$N_all * drop(moments$Q %*% restricted$weights + moments$l)
    stored_gradient <- 2 * drop(moments$Q %*% stored + moments$l)
    residual <- ifelse(stored != 0, abs(stored_gradient + penalty / moments$N_all * sign(stored)),
                       pmax(abs(stored_gradient) - penalty / moments$N_all, 0))
    rows[[length(rows) + 1L]] <- data.frame(task_id = task$task_id, config = task$config,
      p = task$p, rho = task$rho, fold = fold, stored_source_weight = stored[1L],
      exact_source_weight = exact$weights[1L], wald = wald[1L],
      penalty = penalty[1L], penalty_required_for_zero = abs(gradient[1L]),
      objective_increase_if_forced_zero = restricted$objective - exact$objective,
      derivative_from_zero_towards_positive = gradient[1L] + penalty[1L],
      max_weight_difference = difference, max_stored_kkt_residual = max(residual))
  }
}
diagnostics <- do.call(rbind, rows)
write.csv(diagnostics, file.path(output, "independent_penalty_check.csv"), row.names = FALSE)
writeLines(c(
  paste("Independently checked outer-fold problems:", nrow(diagnostics)),
  sprintf("Maximum weight difference versus exhaustive sign enumeration: %.6g", max(diagnostics$max_weight_difference)),
  sprintf("Maximum saved-weight KKT residual on variance scale: %.6g", max(diagnostics$max_stored_kkt_residual)),
  "Other coordinates were reoptimized when forcing the shifted treated-arm weight to zero.",
  "All penalties and objective differences are on the N_all * variance scale."
), file.path(output, "checks.txt"))
print(diagnostics[diagnostics$task_id %in% 1:3 & diagnostics$fold == 1L, ])
cat(readLines(file.path(output, "checks.txt")), sep = "\n")
