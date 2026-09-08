#!/usr/bin/env Rscript
# Paired held-out score check, not a fitted-estimator variance estimate.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 5L) stop("usage: audit_source_tate_projection.R LIBRARY SEED TREATED_PROBE CONTROL_PROBE OUTPUT")
  .libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", equations)
  holdout <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R", holdout)
  read_checked <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  bundle <- read_checked(args[2], "artifacts.rds")
  probes <- list(mu1 = read_checked(args[3], "matrix_probe.rds"),
                 mu0 = read_checked(args[4], "matrix_probe.rds"))
  artifact <- bundle$group_result$artifacts[["0"]]
  fit <- artifact$direct_tate_results$one_round_crossfit
  info <- fit$intermediates$fold_info[[1L]]
  site_ids <- list(target = info$target_idx, source = info$source_idx[[1L]])
  for (arm in 0:1) stopifnot(identical(probes[[paste0("mu", arm)]]$state,
    equations$.saved_projection_state(bundle, arm)))
  results <- summaries <- covariance_rows <- list()
  for (attempt_index in 1:3) {
    fractions <- vapply(probes, function(p) p$attempts[[attempt_index]]$summary$fraction, numeric(1L))
    stopifnot(length(unique(fractions)) == 1L)
    fraction <- unname(fractions[1])
    scores <- lapply(c("target", "source"), function(site) {
      arm_scores <- lapply(1:0, function(arm) {
        probe <- probes[[paste0("mu", arm)]]
        attempt <- probe$attempts[[attempt_index]]
        stopifnot(is.null(attempt$failure))
        data <- artifact$data_split[[if (site == "target") "t" else "s1"]]
        block <- equations$.saved_projection_block(data, site_ids[[site]], arm)
        stopifnot(!any(block$index %in% probe$state[[site]]$index))
        for (piece in probe$state$pieces) for (part in c("training", "calibration"))
          stopifnot(!any(block$index %in% piece[[paste(site, part, sep = "_")]]$index))
        holdout$.projection_holdout_scores(probe$state, block, site, attempt$coefficients)
      })
      names(arm_scores) <- c("mu1", "mu0")
      for (version in c("original", "corrected")) {
        treated <- arm_scores$mu1[[version]]; control <- arm_scores$mu0[[version]]
        contrast <- treated-control
        independent_arm_sum <- var(treated)+var(control)
        cross_arm_term <- -2*cov(treated, control)
        stopifnot(abs(var(contrast)-independent_arm_sum-cross_arm_term) < 1e-10)
        covariance_rows[[length(covariance_rows)+1L]] <<- data.frame(fraction, site, version,
          n = length(contrast), contrast_mean = mean(contrast),
          contrast_score_variance = var(contrast), independent_arm_sum, cross_arm_term)
      }
      data.frame(observation_id = site_ids[[site]],
        original_mu1 = arm_scores$mu1$original, original_mu0 = arm_scores$mu0$original,
        corrected_mu1 = arm_scores$mu1$corrected, corrected_mu0 = arm_scores$mu0$corrected,
        original_tate = arm_scores$mu1$original-arm_scores$mu0$original,
        corrected_tate = arm_scores$mu1$corrected-arm_scores$mu0$corrected)
    })
    names(scores) <- names(site_ids)
    original <- sum(vapply(scores, function(x) mean(x$original_tate), numeric(1L)))
    corrected <- sum(vapply(scores, function(x) mean(x$corrected_tate), numeric(1L)))
    saved <- fit$arm_results$mu1$fold_results[[1L]]$source_results$s1$mu_ts-
      fit$arm_results$mu0$fold_results[[1L]]$source_results$s1$mu_ts
    stopifnot(abs(original-saved) < 1e-12)
    summaries[[attempt_index]] <- data.frame(fraction, original_source_assisted_tate = original,
      corrected_source_assisted_tate = corrected, inference_validated = FALSE)
    results[[attempt_index]] <- scores
  }
  summaries <- do.call(rbind, summaries)
  covariance_rows <- do.call(rbind, covariance_rows)
  roce_write_atomic_directory(args[5], function(stage) {
    write.csv(summaries, file.path(stage, "tate_point_checks.csv"), row.names = FALSE)
    write.csv(covariance_rows, file.path(stage, "cross_arm_covariance.csv"), row.names = FALSE)
    saveRDS(results, file.path(stage, "paired_site_scores.rds"))
    writeLines(c("inference_validated=FALSE", "penalty_policy_selected=FALSE", "observations_removed=0",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/audit_source_tate_projection.R")),
      paste0("input_", 2:4, "_manifest_sha256=", vapply(args[2:4], function(root)
        roce_sha256_file(file.path(root, "sha256.txt")), ""))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "paired source TATE projection audit")
  print(summaries, row.names = FALSE); print(covariance_rows, row.names = FALSE)
}
if (sys.nframe() == 0L) main()
