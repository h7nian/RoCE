#!/usr/bin/env Rscript

# Exploratory evaluation of the derived weight-layer gradient, not resampling.
# The estimator/power are identical to the already evaluated quadratic rule.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: review_n100_weight_derivative.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[1], mustWork = TRUE)
  output <- args[2]
  if (file.exists(output)) stop("weight-layer review output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/quadratic_bias_weights.R")
  source("diagnosis/tate_common_weight/quadratic_weight_influence.R")
  files <- c("diagnosis/tate_common_weight/quadratic_bias_weights.R",
    "diagnosis/tate_common_weight/quadratic_weight_influence.R",
    "diagnosis/tate_common_weight/review_n100_weight_derivative.R")
  hashes <- vapply(files, roce_sha256_file, "")
  candidate_root <- file.path(root, "quadratic_bias_weight_n100_v1")
  candidate_manifest <- readLines(file.path(candidate_root, "sha256.txt"))
  candidate_hash <- substr(candidate_manifest[substring(candidate_manifest, 67L) == "candidate_rows.csv"], 1L, 64L)
  stopifnot(length(candidate_hash) == 1L,
            roce_sha256_file(file.path(candidate_root, "candidate_rows.csv")) == candidate_hash)
  previous <- read.csv(file.path(candidate_root, "candidate_rows.csv"))
  stopifnot(nrow(previous) == 600L, !any(previous$failed), all(previous$power == .75))
  stopifnot(roce_sha256_file(file.path(root, "summaries/n100/sha256.txt")) ==
    "bfe8302a0b1a02405b6e4dab97f0332fce1d844751626dc20113052fda371efe")
  refs <- read.csv(file.path(root, "summaries/n100/input_references.csv"))
  stopifnot(identical(refs$sim_id, 10001:10100))
  rows <- list()
  for (i in 1:100) {
    directory <- file.path(root, sprintf("seed_%06d", 10000L+i))
    stopifnot(roce_sha256_file(file.path(directory, "sha256.txt")) == refs$bundle_checksum_manifest_hash[i])
    manifest <- readLines(file.path(directory, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == "artifacts.rds"], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, "artifacts.rds")) == expected)
    bundle <- readRDS(file.path(directory, "artifacts.rds"))
    for (rho in c(0, .5, 1, 1.5, 2, 2.5)) {
      fit <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
      answer <- .quadratic_weight_influence(fit)
      old <- previous[previous$sim_id == 10000L+i & previous$rho == rho, ]
      stopifnot(nrow(old) == 1L, abs(old$estimate-answer$estimate) < 1e-12,
                abs(old$se^2-answer$fixed_variance) < 1e-12)
      se <- sqrt(answer$weight_linearized_variance)
      rows[[length(rows)+1L]] <- data.frame(sim_id = 10000L+i, rho,
        estimate = answer$estimate, truth = old$truth,
        fixed_weight_se = old$se, weight_layer_se = se,
        fixed_weight_coverage = old$coverage,
        weight_layer_coverage = abs(answer$estimate-old$truth) <= qnorm(.975)*se,
        fixed_variance = answer$fixed_variance,
        indirect_variance = answer$indirect_variance,
        direct_indirect_cross_term = answer$direct_indirect_cross_term,
        weight_layer_variance = answer$weight_linearized_variance,
        max_site_gradient_sum = max(abs(vapply(answer$gradient, sum, numeric(1)))))
    }
    rm(bundle, fit, answer)
    if (i %% 10L == 0L) cat("analytically differentiated seeds:", i, "\n")
  }
  rows <- do.call(rbind, rows)
  summary <- do.call(rbind, lapply(split(rows, rows$rho), function(x) {
    error <- x$estimate-x$truth
    difference <- as.numeric(x$weight_layer_coverage)-as.numeric(x$fixed_weight_coverage)
    data.frame(rho = x$rho[1], n = nrow(x), bias = mean(error),
      empirical_sd = sd(error), rmse = sqrt(mean(error^2)),
      mean_fixed_se = mean(x$fixed_weight_se), mean_weight_layer_se = mean(x$weight_layer_se),
      mean_fixed_variance = mean(x$fixed_variance), mean_weight_layer_variance = mean(x$weight_layer_variance),
      fixed_weight_coverage = mean(x$fixed_weight_coverage),
      weight_layer_coverage = mean(x$weight_layer_coverage),
      coverage_difference = mean(difference), coverage_difference_mcse = sd(difference)/10,
      mean_indirect_variance = mean(x$indirect_variance),
      mean_direct_indirect_cross_term = mean(x$direct_indirect_cross_term),
      negative_cross_terms = sum(x$direct_indirect_cross_term < 0))
  }))
  stopifnot(nrow(rows) == 600L, identical(hashes, vapply(files, roce_sha256_file, "")))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(rows, file.path(stage, "weight_layer_rows.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "rho_summary.csv"), row.names = FALSE)
    write.csv(data.frame(path = files, sha256 = hashes), file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c("analytic_weight_layer_reanalysis=complete", "scope=exploratory_saved_C1_K2_n100",
      "power=0.75", "candidate_point_estimates_unchanged=TRUE",
      "new_mc_replications=0", "resampling_draws=0", "SE_multiplier_used=FALSE",
      "nuisance_fits_differentiated=FALSE", "fold_membership_differentiated=FALSE",
      "independent_validation=FALSE", "inference_validated=FALSE",
      paste0("maximum_site_gradient_sum=", max(rows$max_site_gradient_sum)),
      paste0("candidate_rows_sha256=", candidate_hash)), file.path(stage, "metadata.txt"))
    payloads <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(payloads, roce_sha256_file, ""), basename(payloads), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "analytic weight-layer reanalysis")
  print(summary, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
