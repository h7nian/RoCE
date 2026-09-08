#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) stop("usage: decompose_nested_target_projection.R LIBRARY SEED ORIGINAL_STATES PROJECTION OUTPUT")
  .libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  suppressPackageStartupMessages(library(glmnet))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_target_nuisance_system.R", equations)
  checked_rds <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  target <- checked_rds(args[2], "artifacts.rds")$group_result$artifacts[["0"]]$data_split$t
  fits <- checked_rds(args[3], "target_cv_fit_states.rds")
  states <- checked_rds(args[4], "inner_projection_states.rds")
  limits <- list(ps = c(RoCE:::PROP_SCORE_LOWER, RoCE:::PROP_SCORE_UPPER),
    outcome = c(RoCE:::OUTCOME_PRED_LOWER, RoCE:::OUTCOME_PRED_UPPER))
  summaries <- block_rows <- observation_rows <- list()
  for (key in names(states)) {
    state <- states[[key]]
    pair <- as.integer(strsplit(key, ":", fixed = TRUE)[[1]])
    theta <- unlist(lapply(fits[state$final_nuisance_keys], function(x)
      as.numeric(stats::coef(x, s = "lambda.min"))), use.names = FALSE)
    index <- state$evaluation_ids
    block <- list(W = cbind(1, target$W_outcome[index, , drop = FALSE]),
      Z = cbind(1, target$Z_site[index, , drop = FALSE]), A = target$A[index], Y = target$Y[index])
    system <- equations$.target_tate_system(block, theta, limits)
    coefficients <- state$final_fit$coefficients
    contributions <- vapply(system$index, function(j)
      -drop(system$moment_rows[, j, drop = FALSE] %*% coefficients[j]), numeric(length(index)))
    adjustment <- rowSums(contributions)
    stopifnot(max(abs(system$score-state$original_score)) < 1e-12,
      max(abs(system$score+adjustment-state$adjusted_score)) < 1e-12,
      identical(system$clipping, state$evaluation_clipping))
    variance_change <- var(adjustment)+2*cov(system$score, adjustment)
    stopifnot(abs(variance_change-var(state$adjusted_score)+var(system$score)) < 1e-12)
    centered <- adjustment-mean(adjustment)
    ss <- sum(centered^2)
    summaries[[key]] <- data.frame(outer_fold = pair[1], evaluation_fold = pair[2],
      mean_adjustment = mean(adjustment), adjustment_variance = var(adjustment),
      twice_score_adjustment_covariance = 2*cov(system$score, adjustment),
      score_variance_change = variance_change,
      maximum_absolute_adjustment = max(abs(adjustment)),
      maximum_adjustment_ss_share = if (ss == 0) NA_real_ else max(centered^2)/ss,
      propensity_clipped = system$clipping[["propensity"]],
      outcome_treated_clipped = system$clipping[["outcome_treated"]],
      outcome_control_clipped = system$clipping[["outcome_control"]])
    for (name in names(system$index)) {
      j <- system$index[[name]]
      block_rows[[paste(key, name)]] <- data.frame(outer_fold = pair[1], evaluation_fold = pair[2],
        nuisance_block = name, mean_adjustment = mean(contributions[, name]),
        coefficient_l1 = sum(abs(coefficients[j])), active_coefficients = sum(coefficients[j] != 0),
        contribution_sd = sd(contributions[, name]))
    }
    observation_rows[[key]] <- data.frame(outer_fold = pair[1], evaluation_fold = pair[2],
      observation_id = index, A = block$A, Y = block$Y, propensity = system$predictions$p,
      original_score = system$score, contributions, adjustment, adjusted_score = state$adjusted_score)
  }
  summaries <- do.call(rbind, summaries)
  block_rows <- do.call(rbind, block_rows)
  observation_rows <- do.call(rbind, observation_rows)
  stopifnot(nrow(summaries) == 20L, nrow(block_rows) == 60L, nrow(observation_rows) == 4000L)
  roce_write_atomic_directory(args[5], function(stage) {
    write.csv(summaries, file.path(stage, "pair_decomposition.csv"), row.names = FALSE)
    write.csv(block_rows, file.path(stage, "block_decomposition.csv"), row.names = FALSE)
    write.csv(observation_rows, file.path(stage, "observation_decomposition.csv"), row.names = FALSE)
    writeLines(c("inference_validated=FALSE", "observations_removed=0", "nuisance_refits=0",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/decompose_nested_target_projection.R")),
      paste0("input_", 2:4, "_manifest_sha256=", vapply(args[2:4], function(root)
        roce_sha256_file(file.path(root, "sha256.txt")), ""))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "nested target projection decomposition")
  selected <- paste(summaries$outer_fold, summaries$evaluation_fold, sep = ":") %in% c("2:4", "4:2")
  print(summaries[selected, ], row.names = FALSE)
  print(block_rows[paste(block_rows$outer_fold, block_rows$evaluation_fold, sep = ":") %in% c("2:4", "4:2"), ], row.names = FALSE)
}
if (sys.nframe() == 0L) main()
