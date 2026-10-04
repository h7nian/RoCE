#!/usr/bin/env Rscript
# Baselines on the same generated records and target folds as the RoCE run.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/"))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
.libPaths(c(configuration$library, .libPaths()))
library(RoCE)
stopifnot(identical(normalizePath(find.package("RoCE")),
                    normalizePath(file.path(configuration$library, "RoCE"))))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$task_id == as.integer(args[2L]), , drop = FALSE]
stopifnot(nrow(task) == 1L, length(configuration$baseline_methods) > 0L)
output <- file.path(root, "tasks", args[2L])
stopifnot(dir.exists(output))
options(error = function() {
  writeLines("FAILED", file.path(output, "status.txt"))
  traceback(20L)
  quit(status = 1L)
})
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(configuration$cpus)))
stopifnot(is.finite(cores), cores >= 1L)
invisible(RoCE:::set_nuisance_solver_cpp("proximal_newton"))
if (isTRUE(configuration$nuisance_cv_certificate)) {
  if (!".set_nuisance_cv_certificate" %in% ls(asNamespace("RoCE"), all.names = TRUE)) {
    stop("The frozen library does not support the CV certificate control.")
  }
  invisible(RoCE:::.set_nuisance_cv_certificate(TRUE))
}
writeLines("RUNNING", file.path(output, "status.txt"))
set.seed(task$sim_id)
n_deviated_sites <- if (is.null(configuration$n_deviated_sites)) 1L else configuration$n_deviated_sites
stopifnot(n_deviated_sites >= 0L, n_deviated_sites <= task$K,
          n_deviated_sites == as.integer(n_deviated_sites))
data <- generate_simulation_data(
  n_total = configuration$n_per_site * (task$K + 1L), K = task$K,
  p = task$p, config = task$config, dgp_type = configuration$dgp_type,
  dgp_control = configuration$dgp_control, outcome_type = "binary",
  estimand_type = "superpopulation", ate_deviation = task$rho,
  n_deviated_sites = as.integer(n_deviated_sites), deviation_mechanism = task$deviation_mechanism,
  n_target = configuration$n_per_site,
  n_source_sizes = rep(configuration$n_per_site, task$K), warn_ignored = FALSE)
data_split <- split_data_by_site(data)
stopifnot(all(vapply(data_split, function(site) {
  site$n == configuration$n_per_site && ncol(site$W_outcome) == task$p
}, logical(1L))))
data_hash <- digest::digest(data_split, algo = "sha256")
truth <- data$mu1_true - data$mu0_true
folds <- RoCE:::build_crossfit_folds(data_split, configuration$n_folds)
checkpoint_dir <- file.path(output, "checkpoints")
dir.create(checkpoint_dir, showWarnings = FALSE)
warning_rows <- list()

run_checkpointed_stage <- function(stage_name, compute) {
  path <- file.path(checkpoint_dir, paste0(stage_name, ".rds"))
  input_hash <- digest::digest(list(data_hash, configuration, task, stage_name), algo = "sha256")
  if (file.exists(path)) {
    entry <- readRDS(path)
    stopifnot(identical(entry$input_hash, input_hash),
              identical(entry$value_hash, digest::digest(entry$value, algo = "sha256")))
    assign(".Random.seed", entry$rng_after, envir = .GlobalEnv)
  } else {
    cat("Starting", stage_name, "for task", task$task_id, "at", format(Sys.time()), "\n")
    messages <- character(0)
    value <- withCallingHandlers(compute(), warning = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
    })
    entry <- list(input_hash = input_hash, value = value,
      value_hash = digest::digest(value, algo = "sha256"),
      rng_after = .Random.seed, warnings = messages)
    temporary <- paste0(path, ".tmp.", Sys.getpid())
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(entry, temporary)
    stopifnot(file.rename(temporary, path))
  }
  if (length(entry$warnings)) {
    warning_rows[[stage_name]] <<- data.frame(stage = stage_name, message = entry$warnings)
  }
  entry$value
}

reference <- run_checkpointed_stage("target_reference", function() {
  RoCE:::.target_tate_reference(folds$target_folds, "binomial", "min")
})
dr_weights <- run_checkpointed_stage("density_weights", function() {
  if (!any(configuration$baseline_methods %in% c("federated_dr", "pooled_dr"))) return(NULL)
  RoCE:::.resolve_dr_weights_by_site(data_split, n_cores = cores)
})
arm1 <- run_checkpointed_stage("arm1", function() {
  run_all_comparisons(data_split, n_folds = configuration$n_folds,
    methods = configuration$baseline_methods, include_tilted = FALSE,
    variance_method = configuration$baseline_variance_method,
    n_bootstrap = configuration$n_bootstrap, n_cores = cores,
    dr_weights_by_site = dr_weights, A_val = 1L)
})
fits <- run_checkpointed_stage("tate", function() {
  run_all_comparisons_tate(data_split, n_folds = configuration$n_folds,
    methods = configuration$baseline_methods, include_tilted = FALSE,
    variance_method = configuration$baseline_variance_method,
    n_bootstrap = configuration$n_bootstrap, n_cores = cores,
    dr_weights_by_site = dr_weights, mu1_results = arm1)
})
stopifnot(setequal(names(fits), configuration$baseline_methods))
for (method in names(fits)) {
  fit <- fits[[method]]
  components <- fit$components
  independent_variance <- sum(vapply(fit$influence_blocks, function(block) {
    values <- as.numeric(block$influence)
    block$weight^2 * sum((values - mean(values))^2) / length(values)^2
  }, numeric(1L)))
  covariance_variance <- components$mu1_influence_se^2 + components$mu0_influence_se^2 -
    2 * components$cross_arm_covariance
  stopifnot(abs(fit$estimate - (components$mu1_estimate - components$mu0_estimate)) < 1e-12,
            abs(independent_variance - components$variance_analytic) < 1e-10,
            abs(covariance_variance - components$variance_analytic) < 1e-10)
}

pairing <- list()
identity_fields <- c("config", "K", "p", "rho", "sim_id", "protocol", "deviation_mechanism")
for (reference_root in configuration$baseline_reference_roots) {
  reference_manifest <- read.csv(file.path(reference_root, "manifest.csv"), stringsAsFactors = FALSE)
  matches <- rep(TRUE, nrow(reference_manifest))
  for (field in identity_fields) matches <- matches & reference_manifest[[field]] == task[[field]]
  reference_task <- reference_manifest[matches, , drop = FALSE]
  if (!nrow(reference_task)) next
  stopifnot(nrow(reference_task) == 1L)
  reference_configuration <- jsonlite::fromJSON(file.path(reference_root, "configuration.json"))
  reference_deviated <- if (is.null(reference_configuration$n_deviated_sites)) 1L else reference_configuration$n_deviated_sites
  if (task$rho != 0 && reference_deviated != n_deviated_sites) {
    stop("The paired method uses a different number of deviated sources.")
  }
  reference_output <- file.path(reference_root, "tasks", reference_task$task_id)
  paired <- file.exists(file.path(reference_output, "COMPLETE"))
  if (paired) {
    stopifnot(identical(readLines(file.path(reference_output, "data_sha256.txt")), data_hash))
    reference_results <- read.csv(file.path(reference_output, "results.csv"), stringsAsFactors = FALSE)
    target_row <- reference_results[reference_results$method == "target_only_ate", , drop = FALSE]
    stopifnot(nrow(target_row) == 1L, abs(reference$estimate - target_row$estimate) < 1e-12,
              abs(reference$se - target_row$se) < 1e-12, abs(truth - target_row$truth) < 1e-12)
  }
  pairing[[length(pairing) + 1L]] <- data.frame(reference_root = reference_root,
    reference_task_id = reference_task$task_id, verified = paired)
}
fits$target_only <- reference
rows <- lapply(names(fits), function(method) {
  fit <- fits[[method]]
  stopifnot(is.finite(fit$estimate), is.finite(fit$se), fit$se > 0)
  diagnostics <- RoCE:::.comparison_variance_diagnostics(fit)
  row <- data.frame(sim_id = task$sim_id, method = paste0(method, "_ate"),
    estimand_scope = "tate", estimate = fit$estimate, se = fit$se, truth = truth,
    bias = fit$estimate - truth, coverage = abs(fit$estimate - truth) <= 1.959963984540054 * fit$se,
    n_total = configuration$n_per_site * (task$K + 1L), K = task$K, p = task$p,
    n_deviated_sites = n_deviated_sites,
    config = task$config, n_folds = configuration$n_folds, dgp_type = configuration$dgp_type,
    pilot_recipe = "baselines", pilot_task_id = task$task_id, data_sha256 = data_hash)
  for (name in names(diagnostics)) row[[name]] <- diagnostics[[name]]
  if (isTRUE(configuration$nuisance_cv_certificate)) row$nuisance_cv_certificate <- TRUE
  row
})
results <- do.call(rbind, rows)
saveRDS(list(results = results, fits = fits, data_sha256 = data_hash), file.path(output, "baseline_results.rds"))
write.csv(results, file.path(output, "results.csv"), row.names = FALSE)
if (length(warning_rows)) write.csv(do.call(rbind, warning_rows), file.path(output, "warnings.csv"), row.names = FALSE)
if (length(pairing)) write.csv(do.call(rbind, pairing), file.path(output, "pairing.csv"), row.names = FALSE)
writeLines(data_hash, file.path(output, "data_sha256.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("COMPLETE", file.path(output, "status.txt"))
writeLines("All requested baselines completed and variance identities checked.", file.path(output, "COMPLETE.tmp"))
stopifnot(file.rename(file.path(output, "COMPLETE.tmp"), file.path(output, "COMPLETE")))
unlink(checkpoint_dir, recursive = TRUE)
