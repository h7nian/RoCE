#!/usr/bin/env Rscript

# Prespecified single structural candidate, applied to all saved n=100 C1/K2
# data. Exploratory reweighting only: no nuisance refits and no new MC sample.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: reweight_n100_quadratic_bias.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[1L], mustWork = TRUE)
  output <- args[2L]
  if (file.exists(output)) stop("candidate reweighting output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/quadratic_bias_weights.R")
  scripts <- c("diagnosis/tate_common_weight/quadratic_bias_weights.R",
    "diagnosis/tate_common_weight/reweight_n100_quadratic_bias.R",
    "diagnosis/tate_common_weight/QUADRATIC_BIAS_WEIGHT_CANDIDATE.md")
  before <- vapply(scripts, roce_sha256_file, "")
  summary_root <- file.path(root, "summaries/n100")
  stopifnot(roce_sha256_file(file.path(summary_root, "sha256.txt")) ==
    "bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe")
  references <- read.csv(file.path(summary_root, "input_references.csv"))
  stopifnot(identical(references$sim_id, 10001:10100))
  rows <- weights <- list()
  for (i in 1:100) {
    directory <- file.path(root, sprintf("seed_%06d", 10000L+i))
    manifest <- file.path(directory, "sha256.txt")
    stopifnot(roce_sha256_file(manifest) == references$bundle_checksum_manifest_hash[i])
    lines <- readLines(manifest)
    for (name in c("artifacts.rds", "results.csv")) {
      expected <- substr(lines[substring(lines, 67L) == name], 1L, 64L)
      stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, name)) == expected)
    }
    bundle <- readRDS(file.path(directory, "artifacts.rds"))
    original <- read.csv(file.path(directory, "results.csv"))
    for (rho in c(0, .5, 1, 1.5, 2, 2.5)) {
      fit <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
      baseline <- .common_tate_from_weights(fit, fit$fold_weights)
      stopifnot(abs(baseline$estimate-fit$estimate) < 1e-12,
                abs(baseline$variance-fit$se^2) < 1e-12)
      soft <- original[original$rho == rho & original$method == "one_round_crossfit_ate", ]
      target <- original[original$rho == rho & original$method == "target_only_ate", ]
      stopifnot(nrow(soft) == 1L, nrow(target) == 1L,
                abs(soft$estimate-fit$estimate) < 1e-12)
      candidate_weights <- matrix(NA_real_, 5L, 2L)
      fold <- 0L
      candidate <- tryCatch({
        for (k in 1:5) {
          fold <- k
          inner <- fit$intermediates$inner_fold_info[[k]]
          answer <- .quadratic_bias_weights(inner, power = .75)
          candidate_weights[k, ] <- answer$weights
          weights[[length(weights)+1L]] <- data.frame(sim_id = 10000L+i,
            rho, fold = k, source = c("s1", "s2"), weight = answer$weights,
            baseline_weight = fit$fold_weights[k, ], discrepancy = answer$discrepancy,
            penalty = answer$penalty, target_training_n = answer$target_training_n,
            relative_normal_equation_error = answer$relative_normal_equation_error)
        }
        .common_tate_from_weights(fit, candidate_weights)
      }, error = identity)
      failed <- inherits(candidate, "error")
      estimate <- if (failed) NA_real_ else candidate$estimate
      se <- if (failed) NA_real_ else sqrt(candidate$variance)
      rows[[length(rows)+1L]] <- data.frame(sim_id = 10000L+i, rho,
        method = "quadratic_bias_weight_candidate", power = .75,
        estimate, se, truth = soft$truth,
        coverage = if (failed) NA else abs(estimate-soft$truth) <= qnorm(.975)*se,
        baseline_estimate = soft$estimate, baseline_se = soft$se,
        baseline_coverage = soft$coverage, target_estimate = target$estimate,
        target_se = target$se, target_coverage = target$coverage,
        failed, failed_fold = if (failed) fold else NA_integer_,
        failure = if (failed) conditionMessage(candidate) else "")
    }
    rm(bundle, fit)
    if (i %% 10L == 0L) cat("reweighted retained seeds:", i, "\n")
  }
  rows <- do.call(rbind, rows)
  weights <- do.call(rbind, weights)
  summaries <- lapply(split(rows, rows$rho), function(x) {
    complete <- nrow(x) == 100L && !any(x$failed)
    error <- x$estimate-x$truth
    baseline_error <- x$baseline_estimate-x$truth
    target_error <- x$target_estimate-x$truth
    data.frame(rho = x$rho[1L], n_attempts = nrow(x), failures = sum(x$failed),
      bias = if (complete) mean(error) else NA_real_,
      bias_mcse = if (complete) sd(error)/10 else NA_real_,
      empirical_sd = if (complete) sd(error) else NA_real_,
      mean_se = if (complete) mean(x$se) else NA_real_,
      rmse = if (complete) sqrt(mean(error^2)) else NA_real_,
      coverage = if (complete) mean(x$coverage) else NA_real_,
      baseline_rmse = sqrt(mean(baseline_error^2)), baseline_coverage = mean(x$baseline_coverage),
      target_rmse = sqrt(mean(target_error^2)), target_coverage = mean(x$target_coverage),
      mse_minus_baseline = if (complete) mean(error^2-baseline_error^2) else NA_real_,
      mse_minus_baseline_mcse = if (complete) sd(error^2-baseline_error^2)/10 else NA_real_,
      mse_minus_target = if (complete) mean(error^2-target_error^2) else NA_real_,
      mse_minus_target_mcse = if (complete) sd(error^2-target_error^2)/10 else NA_real_)
  })
  summaries <- do.call(rbind, summaries)
  stopifnot(identical(before, vapply(scripts, roce_sha256_file, "")))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "candidate_rows.csv"), row.names = FALSE)
    write.csv(weights, file.path(stage, "fold_source_weights.csv"), row.names = FALSE)
    write.csv(summaries, file.path(stage, "rho_summary.csv"), row.names = FALSE)
    write.csv(data.frame(path = scripts, sha256 = before), file.path(stage, "input_code_hashes.csv"), row.names = FALSE)
    writeLines(c("quadratic_bias_weight_reanalysis=complete", "scope=exploratory_saved_C1_K2_n100",
      "power=0.75", "power_selected_from_results=FALSE", "new_mc_replications=0",
      "nuisance_refits=0", "bootstrap_used=FALSE", "SE_multiplier_used=FALSE",
      paste0("failed_attempts=", sum(rows$failed)),
      "baseline_reconstruction_passed=TRUE", "original_results_changed=FALSE",
      "independent_validation=FALSE", "inference_validated=FALSE"), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "quadratic bias-weight exploratory reanalysis")
  print(summaries, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
