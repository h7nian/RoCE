#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: check_quadratic_weight_influence.R SEED_BUNDLE OUTPUT")
  directory <- normalizePath(args[1], mustWork = TRUE)
  output <- args[2]
  if (file.exists(output)) stop("weight derivative check output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/quadratic_bias_weights.R")
  source("diagnosis/tate_common_weight/quadratic_weight_influence.R")
  stopifnot(roce_sha256_file(file.path(directory, "sha256.txt")) ==
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  lines <- readLines(file.path(directory, "sha256.txt"))
  expected <- substr(lines[substring(lines, 67L) == "artifacts.rds"], 1L, 64L)
  stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, "artifacts.rds")) == expected)
  bundle <- readRDS(file.path(directory, "artifacts.rds"))
  checks <- results <- list()
  for (rho in c(0, .5, 1, 1.5, 2, 2.5)) {
    fit <- bundle$group_result$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
    result <- .quadratic_weight_influence(fit)
    masses <- lapply(result$gradient, function(x) rep(1, length(x)))
    stopifnot(abs(.quadratic_weight_functional(fit, masses)-result$estimate) < 1e-12)
    for (j in 1:8) {
      direction <- lapply(seq_along(masses), function(s) sin(seq_along(masses[[s]])*j+s))
      epsilon <- 1e-5
      plus <- Map(function(x, d) x+epsilon*d, masses, direction)
      minus <- Map(function(x, d) x-epsilon*d, masses, direction)
      numerical <- (.quadratic_weight_functional(fit, plus)-.quadratic_weight_functional(fit, minus))/(2*epsilon)
      analytical <- sum(unlist(result$gradient)*unlist(direction))
      direct_only <- sum(unlist(result$direct)*unlist(direction))
      checks[[length(checks)+1L]] <- data.frame(rho, direction = j, numerical,
        analytical, absolute_error = abs(numerical-analytical),
        direct_only_error = abs(numerical-direct_only))
      stopifnot(abs(numerical-analytical) < 1e-8)
    }
    results[[as.character(rho)]] <- result
  }
  checks <- do.call(rbind, checks)
  stopifnot(max(checks$direct_only_error) > 1e-5)
  summaries <- do.call(rbind, lapply(names(results), function(rho) {
    x <- results[[rho]]
    data.frame(rho = as.numeric(rho), estimate = x$estimate,
      fixed_weight_se = sqrt(x$fixed_variance),
      weight_linearized_se = sqrt(x$weight_linearized_variance),
      indirect_variance = x$indirect_variance,
      direct_indirect_cross_term = x$direct_indirect_cross_term)
  }))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(checks, file.path(stage, "directional_checks.csv"), row.names = FALSE)
    write.csv(summaries, file.path(stage, "variance_components.csv"), row.names = FALSE)
    saveRDS(results, file.path(stage, "analytic_gradients.rds"))
    writeLines(c("analytic_weight_derivative_check=passed", "sim_id=10013",
      paste0("maximum_directional_error=", max(checks$absolute_error)),
      paste0("maximum_direct_only_error=", max(checks$direct_only_error)),
      "empirical_mass_perturbations=deterministic_derivative_checks_only",
      "resampling_draws=0", "nuisance_fits_differentiated=FALSE",
      "new_mc_replications=0", "inference_validated=FALSE",
      paste0("gradient_code_sha256=", roce_sha256_file("diagnosis/tate_common_weight/quadratic_weight_influence.R"))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "analytic weight derivative check")
  print(summaries, row.names = FALSE)
  cat("maximum_directional_error=", max(checks$absolute_error), "\n")
}

if (sys.nframe() == 0L) main()
