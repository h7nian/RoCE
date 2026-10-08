#!/usr/bin/env Rscript
# Rscript run.R CONFIG.R OUTPUT_DIR [all|roce|sample_size|inverse_variance|federated_dr|pooled_dr|preflight]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) stop("Usage: Rscript run.R CONFIG.R OUTPUT_DIR [METHOD|preflight]")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
source(file.path(dirname(normalizePath(script)), "helpers.R"), local = TRUE)
suppressPackageStartupMessages(library(RoCE))
config_path <- normalizePath(args[1L], mustWork = TRUE)
config <- source(config_path, local = new.env(parent = globalenv()))$value
stopifnot(is.list(config), is.function(config$build_data), is.list(config$fit),
          is.numeric(config$seed), length(config$seed) == 1L)
method <- if (length(args) == 3L) args[3L] else "all"
allowed <- c("roce", "sample_size", "inverse_variance", "federated_dr", "pooled_dr")
if (!method %in% c("all", "preflight", allowed)) stop("Unknown method: ", method)
cores <- as.integer(Sys.getenv("ROCE_CORES", "1"))
if (is.na(cores) || cores < 1L) stop("ROCE_CORES must be a positive integer.")
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
config_hash <- digest::digest(config_path, algo = "sha256", file = TRUE)
package_files <- list.files(find.package("RoCE"), recursive = TRUE, full.names = TRUE)
package_hash <- digest::digest(vapply(package_files, digest::digest, character(1L),
  algo = "sha256", file = TRUE), algo = "sha256")
manifest <- list(config_sha256 = config_hash, package_sha256 = package_hash,
                 package_version = as.character(packageVersion("RoCE")), cores = cores,
                 runner_sha256 = digest::digest(script, algo = "sha256", file = TRUE),
                 helper_sha256 = digest::digest(file.path(dirname(script), "helpers.R"),
                   algo = "sha256", file = TRUE))
manifest_path <- file.path(output, "configuration.rds")
if (file.exists(manifest_path) && !identical(readRDS(manifest_path), manifest)) {
  stop("Output belongs to a different configuration/build/CPU count; choose a new directory.")
}
saveRDS(manifest, manifest_path)
stopifnot(file.copy(config_path, file.path(output, "configuration.R"), overwrite = TRUE))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))

prepare <- function() {
  set.seed(config$seed)
  prepared <- config$build_data()
  data <- prepared$data_split
  RoCE:::validate_algorithm_inputs(data, family = config$fit$family)
  if (any(vapply(data, function(site) !identical(site$W_outcome, site$Z_site), logical(1L)))) {
    stop("This calibrated profile requires the same outcome and weighting feature map.")
  }
  columns <- colnames(data$t$W_outcome)
  if (is.null(columns) || any(vapply(data, function(site)
      !identical(colnames(site$W_outcome), columns), logical(1L)))) {
    stop("All sites must have the same named feature columns in the same order.")
  }
  for (arm in 0:1) RoCE:::validate_crossfit_sample_sizes(data, config$fit$n_folds, A_val = arm)
  data_hash <- digest::digest(data, algo = "sha256")
  path <- file.path(output, "input_sha256.txt")
  if (file.exists(path) && !identical(readLines(path), data_hash)) {
    stop("Input changed since the previous run; choose a new output directory.")
  }
  writeLines(data_hash, path)
  write.csv(data.frame(site = names(data), n = vapply(data, `[[`, integer(1L), "n"),
    treated = vapply(data, function(site) sum(site$A == 1), numeric(1L)),
    control = vapply(data, function(site) sum(site$A == 0), numeric(1L)),
    features = ncol(data$t$W_outcome)), file.path(output, "site_summary.csv"), row.names = FALSE)
  saveRDS(prepared$metadata, file.path(output, "data_metadata.rds"))
  prepared
}
# Validate the complete input before fitting; no patient-level records are exported.
invisible(prepare())
if (method == "preflight") {
  cat("Preflight passed. No nuisance models fitted.\n")
  quit(save = "no", status = 0L)
}
methods <- if (method == "all") c("roce", config$baselines) else method
if (any(!methods %in% allowed)) stop("Invalid configured baseline.")
for (selected in methods) {
  directory <- file.path(output, selected)
  dir.create(directory, showWarnings = FALSE)
  if (file.exists(file.path(directory, "COMPLETE"))) next
  warnings <- character()
  writeLines(c("status: running", paste0("method: ", selected)), file.path(directory, "status.txt"))
  tryCatch({
    prepared <- prepare()  # Reset the same input/RNG stream for each method.
    fit_args <- config$fit
    fit_args$data_split <- prepared$data_split
    fit_args$precomputed_folds <- prepared$folds
    fit_args$n_cores <- cores
    fit_args$parallel_arms <- FALSE
    fit_args$checkpoint_dir <- file.path(directory, "checkpoints")
    capture_warnings <- function(expression) withCallingHandlers(expression, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    if (selected == "roce") {
      fit <- capture_warnings(do.call(RoCE::run_tate_crossfit, fit_args))
      rows <- rbind(real_data_method_row("Target-only", fit$target_only), real_data_method_row("RoCE", fit))
      weights <- fit$weights_by_arm
      if (is.null(weights)) weights <- list(mu1 = fit$weights, mu0 = fit$weights)
      write.csv(data.frame(source = names(weights$mu1), weight_mu1 = weights$mu1,
        weight_mu0 = weights$mu0), file.path(directory, "weights.csv"), row.names = FALSE)
      write.csv(fit$fold_weights, file.path(directory, "fold_weights.csv"), row.names = FALSE)
      write.csv(fit$nuisance_fit_diagnostics, file.path(directory, "nuisance_diagnostics.csv"), row.names = FALSE)
    } else {
      # Match the reported baseline streams; these methods retain their own fitting protocols.
      invisible(RoCE:::.set_nuisance_solver(config$fit$nuisance_solver))
      invisible(RoCE:::.set_nuisance_cv_certificate(config$fit$nuisance_cv_certificate))
      reference <- capture_warnings(RoCE:::.target_tate_reference(prepared$folds$target_folds,
        config$fit$family, config$fit$nuisance_lambda_rule))
      write.csv(real_data_method_row("Ordinary target reference", reference),
        file.path(directory, "ordinary_target_reference.csv"), row.names = FALSE)
      dr_weights <- if (selected %in% c("federated_dr", "pooled_dr")) {
        capture_warnings(RoCE:::.resolve_dr_weights_by_site(prepared$data_split,
          n_cores = cores, caller = "Real-data baseline"))
      } else NULL
      baseline_args <- list(data_split = prepared$comparison_data, family = config$fit$family,
        n_folds = config$fit$n_folds, methods = selected, include_tilted = FALSE,
        variance_method = "bootstrap", n_bootstrap = config$n_bootstrap,
        n_cores = cores, dr_weights_by_site = dr_weights)
      arm1 <- capture_warnings(do.call(RoCE::run_all_comparisons, c(baseline_args, list(A_val = 1L))))
      fit <- capture_warnings(do.call(RoCE::run_all_comparisons_tate,
        c(baseline_args, list(mu1_results = arm1)))[[selected]])
      rows <- real_data_method_row(selected, fit)
      curvature_rows <- real_data_curvature_rows(fit)
      if (nrow(curvature_rows)) write.csv(curvature_rows,
        file.path(directory, "curvature_diagnostics.csv"), row.names = FALSE)
    }
    write.csv(rows, file.path(directory, "methods.csv"), row.names = FALSE)
    # Full fits may contain patient-level scores; keep them only on explicit request.
    if (identical(Sys.getenv("ROCE_SAVE_FIT"), "1")) saveRDS(fit, file.path(directory, "fit.rds"))
    writeLines(c("status: completed", paste0("method: ", selected)), file.path(directory, "status.txt"))
    writeLines("COMPLETE", file.path(directory, "COMPLETE"))
    print(rows)
  }, error = function(error) {
    writeLines(c("status: failed", paste0("method: ", selected),
      paste0("error: ", conditionMessage(error))), file.path(directory, "status.txt"))
    stop(error)
  }, finally = writeLines(warnings, file.path(directory, "warnings.txt")))
}
cat("Finished. Outputs: ", output, "\n", sep = "")
