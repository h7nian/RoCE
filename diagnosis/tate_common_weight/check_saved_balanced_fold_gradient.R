#!/usr/bin/env Rscript

# Apply the independently derived balanced-arm design map to the already
# checked analytic gradients. No choice is made using truth or coverage.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 3L) stop("usage: check_saved_balanced_fold_gradient.R GRADIENT_BUNDLE SEED_BUNDLE OUTPUT")
  gradients_root <- normalizePath(args[1], mustWork = TRUE)
  seed_root <- normalizePath(args[2], mustWork = TRUE)
  output <- args[3]
  if (file.exists(output)) stop("design-map output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/balanced_fold_gradient.R")
  checked_rds <- function(directory, name) {
    lines <- readLines(file.path(directory, "sha256.txt"))
    expected <- substr(lines[substring(lines, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, name)) == expected)
    readRDS(file.path(directory, name))
  }
  stopifnot(roce_sha256_file(file.path(seed_root, "sha256.txt")) ==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  gradients <- checked_rds(gradients_root, "analytic_gradients.rds")
  bundle <- checked_rds(seed_root, "artifacts.rds")
  rows <- records <- list()
  for (rho in names(gradients)) {
    artifact <- bundle$group_result$artifacts[[rho]]
    fit <- artifact$direct_tate_results$one_round_crossfit
    sites <- c("t", names(fit$intermediates$sample_sizes$n_source))
    mapped <- lapply(seq_along(sites), function(s) {
      treatment <- artifact$data_split[[sites[s]]]$A
      fold <- integer(length(treatment))
      for (k in seq_along(fit$intermediates$fold_info)) {
        info <- fit$intermediates$fold_info[[k]]
        ids <- if (s == 1L) info$target_idx else info$source_idx[[s-1L]]
        stopifnot(!any(fold[ids] != 0))
        fold[ids] <- k
      }
      answer <- .balanced_arm_fold_gradient(gradients[[rho]]$gradient[[s]], treatment, fold)
      answer$site <- sites[s]
      answer
    })
    records[[rho]] <- mapped
    rows[[rho]] <- data.frame(rho = as.numeric(rho),
      fixed_weight_se = sqrt(gradients[[rho]]$fixed_variance),
      empirical_mass_weight_layer_se = sqrt(gradients[[rho]]$weight_linearized_variance),
      balanced_arm_weight_layer_se = sqrt(sum(vapply(mapped, `[[`, numeric(1), "variance"))),
      nuisance_fits_differentiated = FALSE,
      integer_fold_remainders_differentiated = FALSE)
  }
  rows <- do.call(rbind, rows)
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "design_comparison.csv"), row.names = FALSE)
    saveRDS(records, file.path(stage, "design_gradients.rds"))
    writeLines(c("balanced_arm_design_mapping=complete", "sim_id=10013",
      "scope=first_order_design_component_not_full_estimator_inference",
      "actual_equal_within_arm_fold_counts_verified=TRUE", "resampling_draws=0",
      "nuisance_fits_differentiated=FALSE", "integer_remainders_differentiated=FALSE",
      "inference_validated=FALSE",
      paste0("input_gradient_manifest_sha256=", roce_sha256_file(file.path(gradients_root, "sha256.txt"))),
      paste0("design_mapper_sha256=", roce_sha256_file("diagnosis/tate_common_weight/balanced_fold_gradient.R"))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "saved balanced-arm design map")
  print(rows, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
