#!/usr/bin/env Rscript
# One frozen RHC analysis per job; nuisance and completed-stage checkpoints survive requeue.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE)
stopifnot(startsWith(root, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
if (!is.null(configuration$collation)) {
  stopifnot(identical(Sys.setlocale("LC_COLLATE", configuration$collation), configuration$collation))
}
.libPaths(c(configuration$library, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = configuration$library))
stopifnot(identical(normalizePath(find.package("RoCE")),
                    normalizePath(file.path(configuration$library, "RoCE"))))
interface <- new.env(parent = asNamespace("RoCE"))
sys.source(file.path(root, "workflow/real_data_rhc.R"), envir = interface)
source(file.path(root, "workflow/layer_comparison_summary.R"))
source(file.path(root, "workflow/rhc_stage_signatures.R"))
manifest <- read.csv(file.path(root, "manifest.csv"), stringsAsFactors = FALSE)
task <- manifest[manifest$task_id == as.integer(arguments[2L]), , drop = FALSE]
stopifnot(nrow(task) == 1L)
output <- file.path(root, "tasks", arguments[2L])
stopifnot(dir.exists(output))

main <- function() {
  stopifnot(identical(digest::digest(rhc_csv_path(), algo = "sha256", file = TRUE),
                      configuration$data_sha256))
  invisible(RoCE:::.set_nuisance_solver("proximal_newton"))
  invisible(RoCE:::.set_nuisance_cv_certificate(TRUE))
  writeLines("RUNNING", file.path(output, "status.txt"))
  set.seed(42L)
  excluded <- c("Medicare & Medicaid", "No insurance")
  recode <- stats::setNames(rep(NA_character_, length(excluded)), excluded)
  covariate_profile <- configuration$covariate_profile
  if (is.null(covariate_profile)) covariate_profile <- "historical"
  preprocessing <- if ("preprocessing" %in% names(task)) task$preprocessing else "cohort"
  preprocessing <- match.arg(preprocessing, c("cohort", "outer_fold"))
  expected_features <- c(historical = 61L, log_missing = 62L, grouped_log_missing = 60L)[[covariate_profile]]
  target_site <- configuration$target_site
  stopifnot(is.character(target_site), length(target_site) == 1L)
  cohort <- interface$build_rhc_cohort(outcome = "death30", site_var = "ninsclas", site_recode = recode,
    covariate_profile = covariate_profile)
  data_split <- interface$build_rhc_data_split(cohort, K = 4L, target_site = target_site, seed = 42L)
  site_mapping <- attr(data_split, "site_mapping")
  expected_sizes <- as.integer(table(cohort$site_var)[unname(site_mapping)])
  stopifnot(nrow(cohort) == 5039L,
    identical(unname(site_mapping["target"]), target_site),
    identical(as.integer(vapply(data_split, `[[`, integer(1L), "n")), expected_sizes),
    all(vapply(data_split, function(site) ncol(site$W_outcome) == expected_features &&
      identical(site$W_outcome, site$Z_site), logical(1L))))
  data_hash <- digest::digest(data_split, algo = "sha256")
  writeLines(data_hash, file.path(output, "data_sha256.txt"))
  checkpoints <- file.path(output, "checkpoints")
  dir.create(checkpoints, showWarnings = FALSE)
  warning_rows <- list()
  stage <- function(name, compute) {
    path <- file.path(checkpoints, paste0(name, ".rds"))
    input_hash <- digest::digest(list(configuration, task, data_hash, name), algo = "sha256")
    if (file.exists(path)) {
      entry <- readRDS(path)
      stopifnot(identical(entry$input_hash, input_hash))
      assign(".Random.seed", entry$rng_after, envir = .GlobalEnv)
    } else {
      messages <- character()
      value <- withCallingHandlers(compute(), warning = function(warning) {
        messages <<- c(messages, conditionMessage(warning))
        invokeRestart("muffleWarning")
      })
      entry <- list(input_hash = input_hash, value = value, rng_after = .Random.seed, warnings = messages)
      temporary <- paste0(path, ".tmp")
      saveRDS(entry, temporary)
      stopifnot(file.rename(temporary, path))
    }
    if (length(entry$warnings)) warning_rows[[name]] <<- data.frame(stage = name, warning = entry$warnings)
    entry$value
  }
  write_table <- function(value, name) {
    destination <- file.path(output, name)
    temporary <- paste0(destination, ".tmp")
    write.csv(value, temporary, row.names = FALSE)
    stopifnot(file.rename(temporary, destination))
  }
  cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", as.character(configuration$cpus)))
  stopifnot(cores >= 2L)
  if (task$role == "method") {
    source_rules <- NULL
    rule_columns <- c(initial_weight = "source_initial_weight_rule", weight = "source_weight_rule",
                      outcome = "source_outcome_rule")
    if (all(rule_columns %in% names(task))) {
      values <- vapply(rule_columns, function(column) as.character(task[[column]]), character(1L))
      selected <- !is.na(values) & nzchar(values)
      if (any(selected)) source_rules <- as.list(values[selected])
    }
    target_radius <- if ("target_radius" %in% names(task)) task$target_radius else task$radius
    stopifnot(is.finite(target_radius), target_radius > 0)
    control <- interface$.rhc_calibration_control(task$target_program, task$source_program, target_radius,
      source_lambda_rules = source_rules)
    fit <- stage("analysis", function() interface$run_rhc_tate_experiment(
      K = 4L, target_site = target_site, outcome = "death30", n_folds = 10L,
      site_var = "ninsclas", site_recode = recode, seed = 42L,
      covariate_profile = covariate_profile,
      preprocessing = preprocessing,
      fold_seed = if (task$fold_seed == 0L) NULL else as.integer(task$fold_seed),
      aggregation_mode = "joint_tate", aggregation_lambda = .5,
      target_nuisance_method = task$target_program, source_validation_method = "outer_fit",
      crossfit_layers = 2L, calibration_control = control, calibration_layout = "compact",
      nuisance_solver = "proximal_newton", nuisance_tol = 1e-10, nuisance_cv_certificate = TRUE,
      M_tau = task$radius, M_tau_inference = task$radius, nlambda_init = 100L,
      nuisance_lambda_rule = task$nuisance_rule, comparison_methods = character(),
      checkpoint_dir = file.path(checkpoints, "nuisance"),
      n_cores = as.integer(cores %/% 2L), parallel_arms = TRUE, verbose = TRUE))
    stopifnot(identical(fit$metadata$data_sha256, data_hash), fit$metadata$K == 3L,
      identical(unname(fit$metadata$target_site), target_site),
      fit$metadata$n_sites == 4L, fit$metadata$crossfit_layers == 2L,
      all(is.finite(fit$methods$estimate)), all(is.finite(fit$methods$se)), all(fit$methods$se > 0))
    write_table(fit$methods, "methods.csv")
    write_table(fit$pairwise, "sources.csv")
    write_table(fit$tate_fit$nuisance_fit_diagnostics, "nuisance_fits.csv")
    jsonlite::write_json(fit$metadata, file.path(output, "metadata.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    arm <- summarize_layer_arms(fit$tate_fit, data_split, keep_influence = TRUE)
    write_table(data.frame(arm = names(arm$estimates), estimate = as.numeric(arm$estimates),
      variance = diag(arm$covariance), variance_fixed_weights = diag(arm$fixed_covariance),
      covariance_mu1_mu0 = arm$covariance[1L, 2L],
      covariance_mu1_mu0_fixed_weights = arm$fixed_covariance[1L, 2L]), "arm_summary.csv")
    saveRDS(arm, file.path(output, "arm_influence.rds"))
    saveRDS(rhc_stage_signatures(fit$tate_fit), file.path(output, "stage_signatures.rds"))
    saveRDS(fit$tate_fit$clip_diagnostics, file.path(output, "clip_diagnostics.rds"))
    if (task$config == "RHC_primary") {
      sensitivity <- stage("reaggregation", function() {
        do.call(rbind, lapply(c("joint_tate", "separate_arms", "common_tate"), function(mode) {
          do.call(rbind, lapply(c(1, 2, 3), function(cutoff) {
            selected <- reaggregate_tate_crossfit(data_split, fit$tate_fit,
              M_tau_inference = task$radius, lambda_selection = 1 / cutoff, aggregation_mode = mode)
            data.frame(aggregation_mode = mode, cutoff = cutoff, estimate = selected$estimate,
              se = selected$se, se_fixed_weights = selected$se_fixed_weights,
              ci_lower = selected$ci_lower, ci_upper = selected$ci_upper)
          }))
        }))
      })
      write_table(sensitivity, "aggregation_sensitivity.csv")
    }
  } else {
    stopifnot(task$role == "baseline", task$baseline %in% configuration$baseline_methods)
    folds <- RoCE:::build_crossfit_folds(data_split, 10L)
    comparison_data <- data_split
    if (preprocessing == "outer_fold") {
      prepared <- RoCE:::.prepare_rhc_outer_preprocessing(cohort, data_split, folds,
        cohort_arguments = list(outcome = "death30", site_var = "ninsclas",
          site_recode = recode, covariate_profile = covariate_profile),
        split_arguments = list(K = 4L, target_site = target_site, phi = base::identity, seed = 42L))
      folds <- prepared$folds
      comparison_data <- prepared$baseline_data
      saveRDS(attr(folds, "preprocessing"), file.path(output, "preprocessing.rds"))
    }
    reference <- stage("target_reference", function() RoCE:::.target_tate_reference(folds$target_folds, "binomial", "min"))
    dr_weights <- stage("density_weights", function() {
      if (task$baseline %in% c("federated_dr", "pooled_dr")) {
        RoCE:::.resolve_dr_weights_by_site(data_split, n_cores = cores, caller = "RHC baseline")
      } else NULL
    })
    arm1 <- stage("arm1", function() run_all_comparisons(comparison_data, n_folds = 10L,
      methods = task$baseline, include_tilted = FALSE, variance_method = "bootstrap",
      n_bootstrap = 5000L, n_cores = cores, dr_weights_by_site = dr_weights, A_val = 1L))
    fit <- stage("analysis", function() run_all_comparisons_tate(comparison_data, n_folds = 10L,
      methods = task$baseline, include_tilted = FALSE, variance_method = "bootstrap",
      n_bootstrap = 5000L, n_cores = cores, dr_weights_by_site = dr_weights, mu1_results = arm1)[[task$baseline]])
    components <- fit$components
    stopifnot(is.finite(fit$estimate), is.finite(fit$se), fit$se > 0,
      abs(fit$estimate - (components$mu1_estimate - components$mu0_estimate)) < 1e-12,
      abs(components$mu1_influence_se^2 + components$mu0_influence_se^2 -
        2 * components$cross_arm_covariance - components$variance_analytic) < 1e-10)
    write_table(rbind(interface$.rhc_method_row("Target-only", reference$estimate, reference$se),
      interface$.rhc_method_row(task$baseline, fit$estimate, fit$se)), "methods.csv")
    write_table(as.data.frame(RoCE:::.comparison_variance_diagnostics(fit)), "comparison_diagnostics.csv")
  }
  if (length(warning_rows)) write_table(do.call(rbind, warning_rows), "warnings.csv")
  writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
  writeLines("COMPLETE", file.path(output, "status.txt"))
  writeLines("RHC analysis completed; no simulation truth or coverage is defined.", file.path(output, "COMPLETE.tmp"))
  stopifnot(file.rename(file.path(output, "COMPLETE.tmp"), file.path(output, "COMPLETE")))
}

tryCatch(main(), error = function(error) {
  writeLines(conditionMessage(error), file.path(output, "analysis_FAILED.txt"))
  writeLines("FAILED", file.path(output, "status.txt"))
  stop(error)
})
