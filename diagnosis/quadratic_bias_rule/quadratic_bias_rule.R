#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0007 Land the smooth quadratic-bias weight rule as the sensitivity estimator
# Task:    Replay the 100 saved v19 seeds (C1, K = 2, six rho values) with the
#          in-package quadratic-bias rule: re-solve the fold weights from the stored
#          inner-fold moments with RoCE:::.quadratic_bias_weights(), rebuild the
#          pseudo-values, and compute the weight-layer variance with the quadratic
#          derivative. Compare with the Stage-1 candidate rows. No refits.
#          Writes diagnosis/out/quadratic_bias_rule/v1.

suppressPackageStartupMessages(library(RoCE))
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

PILOT_ROOT <- "results/direct_tate_mc500_b5000/independent_inference_pilot_v19"
PROTOTYPE_ROWS <- file.path(PILOT_ROOT, "quadratic_weight_layer_n100_v1", "weight_layer_rows.csv")
PROTOTYPE_WEIGHTS <- file.path(PILOT_ROOT, "quadratic_bias_weight_n100_v1", "fold_source_weights.csv")
RHO_VALUES <- c(0, 0.5, 1, 1.5, 2, 2.5)

.package_quadratic_replay <- function(fit) {
  sizes <- fit$intermediates$sample_sizes
  n_t <- sizes$n_t
  n_source <- sizes$n_source
  K <- length(n_source)
  mass <- lapply(c(n_t, n_source), function(n) rep(1, n))
  inner <- fit$intermediates$inner_fold_info
  weights <- do.call(rbind, lapply(inner, function(info) {
    moments <- RoCE:::.inner_fold_moments(RoCE:::.inner_fold_records(info), mass)$moments
    RoCE:::.quadratic_bias_weights(moments)
  }))
  phi <- RoCE:::.compute_phase3_all_phi(
    n_folds = fit$n_folds, fold_weights = weights,
    fold_info = fit$intermediates$fold_info,
    n_t = n_t, n_source_full = n_source, N_all = n_t + sum(n_source),
    source_sites = paste0("s", seq_len(K)), K = K, verbose = FALSE
  )
  gradient <- RoCE:::.weight_layer_gradient(
    fit$intermediates$fold_info, inner, weights, fit$fold_lambdas,
    screening_rule = "quadratic_bias", rep(0, fit$n_folds), n_t, n_source
  )
  list(
    estimate = mean(phi), weights = weights,
    layer = RoCE:::.weight_layer_variance(phi, n_t, n_source, gradient)
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: quadratic_bias_rule.R OUTPUT_DIRECTORY")
  output <- file.path(args[1], "v1")
  if (file.exists(output)) stop("output already exists: ", output)
  prototype <- read.csv(PROTOTYPE_ROWS)
  prototype_weights <- read.csv(PROTOTYPE_WEIGHTS)
  rows <- list()
  for (sim_id in 10001:10100) {
    bundle <- readRDS(file.path(PILOT_ROOT, sprintf("seed_%06d", sim_id), "artifacts.rds"))
    for (rho in RHO_VALUES) {
      fit <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
      replay <- .package_quadratic_replay(fit)
      reference <- prototype[prototype$sim_id == sim_id & prototype$rho == rho, ]
      stopifnot(nrow(reference) == 1L)
      reference_weights <- prototype_weights[
        prototype_weights$sim_id == sim_id & prototype_weights$rho == rho,
      ]
      stopifnot(nrow(reference_weights) == length(replay$weights))
      reference_matrix <- matrix(NA_real_, nrow(replay$weights), ncol(replay$weights))
      reference_matrix[cbind(reference_weights$fold, match(reference_weights$source, paste0("s", seq_len(ncol(replay$weights)))))] <-
        reference_weights$weight
      rows[[length(rows) + 1L]] <- data.frame(
        sim_id = sim_id, rho = rho,
        estimate = replay$estimate, truth = reference$truth,
        weight_layer_se = sqrt(replay$layer$variance),
        fixed_weight_se = sqrt(replay$layer$fixed_variance),
        estimate_difference = abs(replay$estimate - reference$estimate),
        weight_difference = max(abs(replay$weights - reference_matrix)),
        se_difference = abs(sqrt(replay$layer$variance) - reference$weight_layer_se),
        fixed_se_difference = abs(sqrt(replay$layer$fixed_variance) - reference$fixed_weight_se),
        kink_cells = replay$layer$kink_cells
      )
    }
    rm(bundle)
  }
  rows <- do.call(rbind, rows)
  z <- stats::qnorm(0.975)
  rows$covered <- abs(rows$estimate - rows$truth) <= z * rows$weight_layer_se
  summary <- do.call(rbind, lapply(split(rows, rows$rho), function(g) data.frame(
    rho = g$rho[1L], n = nrow(g), bias = mean(g$estimate - g$truth),
    empirical_sd = stats::sd(g$estimate), mean_weight_layer_se = mean(g$weight_layer_se),
    coverage = mean(g$covered)
  )))
  print(summary, row.names = FALSE)
  cat(sprintf(
    "max |package - prototype|: estimate %.3e, weights %.3e, weight-layer SE %.3e, fixed SE %.3e\n",
    max(rows$estimate_difference), max(rows$weight_difference),
    max(rows$se_difference), max(rows$fixed_se_difference)
  ))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "package_vs_prototype.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "rho_summary.csv"), row.names = FALSE)
    writeLines(c("task=quadratic_bias_rule", "history_entry=0007",
                 paste0("max_estimate_difference=", max(rows$estimate_difference)),
                 paste0("max_weight_difference=", max(rows$weight_difference)),
                 paste0("max_weight_layer_se_difference=", max(rows$se_difference)),
                 paste0("max_fixed_se_difference=", max(rows$fixed_se_difference)),
                 paste0("package_library=", .libPaths()[1])),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "quadratic-bias rule package replay")
}

if (sys.nframe() == 0L) main()
