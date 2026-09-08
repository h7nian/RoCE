#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0005 Land the weight-layer variance in the package
# Task:    Replay the 100 saved v19 seeds with the in-package weight layer
#          (RoCE:::.weight_layer_gradient / .weight_layer_variance) and compare
#          with the Stage-1 prototype output diagnosis/out/weight_layer/v1
#          (HISTORY #0001). Writes diagnosis/out/weight_layer/v2. No refits.

suppressPackageStartupMessages(library(RoCE))
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

PILOT_ROOT <- "results/direct_tate_mc500_b5000/independent_inference_pilot_v19"
PROTOTYPE_ROWS <- "diagnosis/out/weight_layer/v1/replay_rows.csv"
RHO_VALUES <- c(0, 0.5, 1, 1.5, 2, 2.5)

.package_weight_layer <- function(fit) {
  sizes <- fit$intermediates$sample_sizes
  gradient <- RoCE:::.weight_layer_gradient(
    fit$intermediates$fold_info, fit$intermediates$inner_fold_info,
    fit$fold_weights, fit$fold_lambdas, penalized = TRUE,
    fit$fold_weight_psd_ridge, sizes$n_t, sizes$n_source
  )
  RoCE:::.weight_layer_variance(fit$all_phi_agg, sizes$n_t, sizes$n_source, gradient)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: weight_layer.R OUTPUT_DIRECTORY")
  output <- file.path(args[1], "v2")
  if (file.exists(output)) stop("output already exists: ", output)
  prototype <- read.csv(PROTOTYPE_ROWS)
  rows <- list()
  for (sim_id in 10001:10100) {
    bundle <- readRDS(file.path(PILOT_ROOT, sprintf("seed_%06d", sim_id), "artifacts.rds"))
    for (rho in RHO_VALUES) {
      fit <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
      layer <- .package_weight_layer(fit)
      reference <- prototype[prototype$sim_id == sim_id & prototype$rho == rho, ]
      stopifnot(nrow(reference) == 1L)
      rows[[length(rows) + 1L]] <- data.frame(
        sim_id = sim_id, rho = rho,
        package_weight_layer_se = sqrt(layer$variance),
        prototype_weight_layer_se = reference$weight_layer_se,
        se_difference = abs(sqrt(layer$variance) - reference$weight_layer_se),
        fixed_se_difference = abs(sqrt(layer$fixed_variance) - reference$fixed_weight_se),
        kink_cells = layer$kink_cells
      )
    }
    rm(bundle)
  }
  rows <- do.call(rbind, rows)
  cat(sprintf("max |package - prototype| weight-layer SE: %.3e; fixed SE: %.3e\n",
              max(rows$se_difference), max(rows$fixed_se_difference)))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "package_vs_prototype.csv"), row.names = FALSE)
    writeLines(c("task=weight_layer", "history_entry=0005",
                 paste0("max_weight_layer_se_difference=", max(rows$se_difference)),
                 paste0("max_fixed_se_difference=", max(rows$fixed_se_difference)),
                 paste0("package_library=", .libPaths()[1])),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "weight layer package replay")
}

if (sys.nframe() == 0L) main()
