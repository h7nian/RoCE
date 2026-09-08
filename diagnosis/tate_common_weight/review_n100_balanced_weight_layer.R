#!/usr/bin/env Rscript

# Apply the actual equal-within-treatment fold-design map to every saved fit.
# No candidate choice, nuisance refit, resampling, or new point estimate.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: review_n100_balanced_weight_layer.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[1], mustWork = TRUE); output <- args[2]
  if (file.exists(output)) stop("balanced weight-layer review output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/quadratic_bias_weights.R")
  source("diagnosis/tate_common_weight/quadratic_weight_influence.R")
  source("diagnosis/tate_common_weight/balanced_fold_gradient.R")
  files <- c("diagnosis/tate_common_weight/quadratic_bias_weights.R",
    "diagnosis/tate_common_weight/quadratic_weight_influence.R",
    "diagnosis/tate_common_weight/balanced_fold_gradient.R",
    "diagnosis/tate_common_weight/review_n100_balanced_weight_layer.R")
  hashes <- vapply(files, roce_sha256_file, "")
  old_root <- file.path(root, "quadratic_weight_layer_n100_v1")
  lines <- readLines(file.path(old_root, "sha256.txt"))
  old_hash <- substr(lines[substring(lines, 67L) == "weight_layer_rows.csv"], 1L, 64L)
  stopifnot(length(old_hash) == 1L, roce_sha256_file(file.path(old_root, "weight_layer_rows.csv")) == old_hash)
  previous <- read.csv(file.path(old_root, "weight_layer_rows.csv"))
  stopifnot(nrow(previous) == 600L,
    roce_sha256_file(file.path(root, "summaries/n100/sha256.txt")) ==
      "bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe")
  refs <- read.csv(file.path(root, "summaries/n100/input_references.csv"))
  stopifnot(identical(refs$sim_id, 10001:10100))
  rows <- list()
  for (i in 1:100) {
    directory <- file.path(root, sprintf("seed_%06d", 10000L+i))
    stopifnot(roce_sha256_file(file.path(directory, "sha256.txt")) == refs$bundle_checksum_manifest_hash[i])
    lines <- readLines(file.path(directory, "sha256.txt"))
    expected <- substr(lines[substring(lines, 67L) == "artifacts.rds"], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, "artifacts.rds")) == expected)
    bundle <- readRDS(file.path(directory, "artifacts.rds"))
    for (rho in c(0, .5, 1, 1.5, 2, 2.5)) {
      artifact <- bundle$group_result$artifacts[[as.character(rho)]]
      fit <- artifact$direct_tate_results$one_round_crossfit
      answer <- .quadratic_weight_influence(fit)
      old <- previous[previous$sim_id == 10000L+i & previous$rho == rho, ]
      stopifnot(nrow(old) == 1L, abs(answer$estimate-old$estimate) < 1e-12,
                abs(answer$weight_linearized_variance-old$weight_layer_variance) < 1e-12)
      sites <- c("t", names(fit$intermediates$sample_sizes$n_source))
      mapped <- lapply(seq_along(sites), function(s) {
        A <- artifact$data_split[[sites[s]]]$A
        fold <- integer(length(A))
        for (k in seq_along(fit$intermediates$fold_info)) {
          info <- fit$intermediates$fold_info[[k]]
          index <- if (s == 1L) info$target_idx else info$source_idx[[s-1L]]
          stopifnot(!any(fold[index] != 0))
          fold[index] <- k
        }
        .balanced_arm_fold_gradient(answer$gradient[[s]], A, fold)
      })
      variance <- sum(vapply(mapped, `[[`, numeric(1L), "variance"))
      rows[[length(rows)+1L]] <- data.frame(sim_id = 10000L+i, rho,
        estimate = answer$estimate, truth = old$truth,
        fixed_weight_se = old$fixed_weight_se, mass_weight_layer_se = old$weight_layer_se,
        balanced_weight_layer_se = sqrt(variance),
        fixed_weight_coverage = old$fixed_weight_coverage,
        mass_weight_layer_coverage = old$weight_layer_coverage,
        balanced_weight_layer_coverage = abs(answer$estimate-old$truth) <= qnorm(.975)*sqrt(variance),
        within_cell_variance = sum(vapply(mapped, function(x) sum(x$within_cell^2), numeric(1))),
        treatment_total_variance = sum(vapply(mapped, function(x) sum(x$treatment_composition^2), numeric(1))))
    }
    rm(bundle, artifact, fit, answer, mapped)
    if (i %% 10L == 0L) cat("verified balanced-design seeds:", i, "\n")
  }
  rows <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(rows, rows$rho), function(x) {
    error <- x$estimate-x$truth
    data.frame(rho = x$rho[1], n = nrow(x), bias = mean(error), empirical_sd = sd(error),
      rmse = sqrt(mean(error^2)), mean_mass_se = mean(x$mass_weight_layer_se),
      mean_balanced_se = mean(x$balanced_weight_layer_se),
      fixed_coverage = mean(x$fixed_weight_coverage),
      mass_coverage = mean(x$mass_weight_layer_coverage),
      balanced_coverage = mean(x$balanced_weight_layer_coverage),
      mean_within_cell_variance = mean(x$within_cell_variance),
      mean_treatment_total_variance = mean(x$treatment_total_variance))
  }))
  stopifnot(nrow(rows) == 600L, identical(hashes, vapply(files, roce_sha256_file, "")),
            max(abs(rows$balanced_weight_layer_se^2-rows$within_cell_variance-rows$treatment_total_variance)) < 1e-12)
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "design_rows.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "rho_summary.csv"), row.names = FALSE)
    write.csv(data.frame(path = files, sha256 = hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c("balanced_weight_layer_reanalysis=complete", "scope=exploratory_saved_C1_K2_n100",
      "candidate_points_unchanged=TRUE", "power=0.75", "all_site_arm_fold_counts_verified=TRUE",
      "resampling_draws=0", "nuisance_refits=0", "nuisance_fits_differentiated=FALSE",
      "integer_fold_remainders_differentiated=FALSE", "independent_validation=FALSE",
      "inference_validated=FALSE", paste0("previous_rows_sha256=", old_hash)), file.path(stage, "metadata.txt"))
    payloads <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(payloads, roce_sha256_file, ""), basename(payloads), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "balanced weight-layer review")
  print(summary, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
