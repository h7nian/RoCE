#!/usr/bin/env Rscript

main <- function() {
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "resource_topology.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
package_provenance <- roce_runtime_package_provenance(project_library)
loaded_library <- package_provenance$library
package_fingerprint <- package_provenance$fingerprint
rhc_workflow_fingerprint <- roce_sha256_environment(
  "ROCE_RHC_WORKFLOW_FINGERPRINT"
)
rhc_data_fingerprint <- roce_sha256_environment(
  "ROCE_RHC_DATA_FINGERPRINT"
)
rhc_data_path <- normalizePath(rhc_csv_path(), mustWork = TRUE)
rhc_sha_output <- system2("sha256sum", rhc_data_path, stdout = TRUE,
                          stderr = TRUE)
rhc_sha_status <- attr(rhc_sha_output, "status")
if (is.null(rhc_sha_status)) rhc_sha_status <- 0L
observed_rhc_data_fingerprint <- if (
    rhc_sha_status == 0L && length(rhc_sha_output) == 1L) {
  strsplit(rhc_sha_output[[1L]], "[[:space:]]+")[[1L]][[1L]]
} else {
  NA_character_
}
if (!identical(
    tolower(observed_rhc_data_fingerprint), rhc_data_fingerprint
)) {
  stop("loaded RHC data fingerprint does not match the locked input.")
}

output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT",
  "results/direct_tate_mc500_b5000/rhc_supported_primary"
)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
artifact_names <- c(
  "rhc_direct_tate.rds",
  "rhc_direct_tate_methods.csv",
  "rhc_direct_tate_sources.csv",
  "rhc_direct_tate_metadata.csv",
  "rhc_direct_tate_diagnostics.csv",
  "rhc_direct_tate_comparison_diagnostics.csv",
  "rhc_direct_tate_nuisance_fits.csv",
  "rhc_direct_tate_forest.pdf"
)
output_files <- file.path(output_root, artifact_names)
existing_outputs <- output_files[file.exists(output_files)]
if (length(existing_outputs) > 0L) {
  stop(
    "RHC output directory already contains result artifacts; use a fresh ",
    "ROCE_OUTPUT_ROOT to preserve existing results: ",
    paste(existing_outputs, collapse = ", ")
  )
}
lock_path <- file.path(output_root, ".rhc_direct_tate.lock")
release_lock <- roce_claim_task_lock(lock_path, c(
  paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", "local")),
  paste0("claimed_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
))
on.exit(release_lock(), add = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required to produce the complete RHC artifact set.")
}
staging_directory <- tempfile(pattern = ".rhc_staging_", tmpdir = output_root)
if (!dir.create(staging_directory, showWarnings = FALSE)) {
  stop("failed to create the RHC staging directory: ", staging_directory)
}
committed_files <- character(0)
commit_complete <- FALSE
on.exit({
  if (!commit_complete && length(committed_files) > 0L) {
    unlink(committed_files, force = TRUE)
  }
  unlink(staging_directory, recursive = TRUE, force = TRUE)
}, add = TRUE)
artifact_path <- function(name) file.path(staging_directory, name)

allocated_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(allocated_cores) || allocated_cores < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.")
}
n_folds <- 10L
total_sites <- suppressWarnings(as.integer(
  Sys.getenv("ROCE_RHC_TOTAL_SITES", "4")
))
if (length(total_sites) != 1L || is.na(total_sites) || total_sites < 2L) {
  stop("ROCE_RHC_TOTAL_SITES must be one integer >= 2.")
}
nuisance_cv_threads <- roce_read_positive_integer_env(
  "ROCE_NUISANCE_CV_THREADS", 1L, n_folds
)
resource_plan <- roce_slurm_resource_plan(
  source_count = total_sites - 1L,
  n_folds = n_folds,
  allocated_cores = allocated_cores,
  nuisance_cv_threads = nuisance_cv_threads
)
parallel_arms <- resource_plan$parallel_treatment_arms
source_cores_per_arm <- resource_plan$source_workers
nuisance_lambda_rule <- Sys.getenv("ROCE_NUISANCE_LAMBDA_RULE", "min")
if (!nuisance_lambda_rule %in% c("min", "1se")) {
  stop("ROCE_NUISANCE_LAMBDA_RULE must be 'min' or '1se'.")
}
n_bootstrap <- suppressWarnings(as.integer(
  Sys.getenv("ROCE_N_BOOTSTRAP", "5000")
))
if (length(n_bootstrap) != 1L || is.na(n_bootstrap) || n_bootstrap < 2L) {
  stop("ROCE_N_BOOTSTRAP must be one integer >= 2.")
}
nlambda_init <- suppressWarnings(as.integer(
  Sys.getenv("ROCE_NLAMBDA_INIT", "100")
))
if (length(nlambda_init) != 1L || is.na(nlambda_init) ||
    nlambda_init < 2L) {
  stop("ROCE_NLAMBDA_INIT must be one integer >= 2.")
}
excluded_site_text <- trimws(Sys.getenv(
  "ROCE_RHC_EXCLUDE_SITES",
  "Medicare & Medicaid,No insurance"
))
excluded_sites <- if (nzchar(excluded_site_text)) {
  unique(trimws(strsplit(excluded_site_text, ",", fixed = TRUE)[[1L]]))
} else {
  character(0)
}
if (any(!nzchar(excluded_sites))) {
  stop("ROCE_RHC_EXCLUDE_SITES contains an empty comma-separated level.")
}
site_recode <- if (length(excluded_sites) > 0L) {
  stats::setNames(rep(NA_character_, length(excluded_sites)), excluded_sites)
} else {
  NULL
}
message(sprintf(
  paste0(
    "[resource] allocated=%d cv_threads=%d source_workers_per_arm=%d ",
    "parallel_arms=%s"
  ),
  resource_plan$allocated_cores,
  resource_plan$nuisance_cv_threads,
  resource_plan$source_workers,
  resource_plan$parallel_treatment_arms
))

result <- run_rhc_tate_experiment(
  # The production driver defaults to the pre-specified empirical support
  # screen: four retained strata (one target plus three sources) after dropping
  # both raw components of the unsupported collapsed source. Environment
  # overrides remain available for separately named sensitivity analyses.
  K = total_sites,
  target_site = "Private",
  outcome = "death30",
  site_var = "ninsclas",
  site_recode = site_recode,
  n_folds = n_folds,
  nlambda_init = nlambda_init,
  nuisance_lambda_rule = nuisance_lambda_rule,
  seed = 42L,
  n_cores = source_cores_per_arm,
  parallel_arms = parallel_arms,
  M_tau = 5,
  M_tau_inference = 5,
  aggregation_lambda = RoCE:::AGG_WALD_LAMBDA,
  comparison_methods = c(
    "sample_size", "inverse_variance", "federated_dr", "pooled_dr"
  ),
  variance_method = "bootstrap",
  n_bootstrap = n_bootstrap,
  verbose = TRUE
)

result$metadata$package_library <- loaded_library
result$metadata$package_fingerprint <- package_fingerprint
result$metadata$rhc_workflow_fingerprint <- rhc_workflow_fingerprint
result$metadata$rhc_data_path <- rhc_data_path
result$metadata$rhc_data_fingerprint <- rhc_data_fingerprint
result$metadata$target_anchor_nlambda <- RoCE:::LAMBDA_GRID_SIZE_STANDARD
result$metadata$allocated_cores <- resource_plan$allocated_cores
result$metadata$nuisance_cv_threads <- resource_plan$nuisance_cv_threads
result$metadata$source_workers_per_arm <- resource_plan$source_workers
result$metadata$parallel_treatment_arms <-
  resource_plan$parallel_treatment_arms
scheduler_provenance <- roce_scheduler_provenance()
for (field in names(scheduler_provenance)) {
  result$metadata[[field]] <- scheduler_provenance[[field]]
}

saveRDS(result, artifact_path("rhc_direct_tate.rds"))
write.csv(
  result$methods,
  artifact_path("rhc_direct_tate_methods.csv"),
  row.names = FALSE
)
write.csv(
  result$pairwise,
  artifact_path("rhc_direct_tate_sources.csv"),
  row.names = FALSE
)

metadata <- data.frame(
  name = names(unlist(result$metadata, recursive = FALSE)),
  value = vapply(
    unlist(result$metadata, recursive = FALSE),
    function(value) paste(value, collapse = ","),
    character(1L)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  artifact_path("rhc_direct_tate_metadata.csv"),
  row.names = FALSE
)

diagnostics <- c(
  RoCE:::.summarize_tate_crossfit_timing(
    result$tate_fit$arm_results$mu1,
    result$tate_fit$arm_results$mu0
  ),
  RoCE:::.summarize_tate_nuisance_diagnostics(
    result$tate_fit$arm_results$mu1,
    result$tate_fit$arm_results$mu0
  ),
  RoCE:::.summarize_tate_aggregation_diagnostics(result$tate_fit)
)
write.csv(
  as.data.frame(as.list(diagnostics), check.names = FALSE),
  artifact_path("rhc_direct_tate_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  result$tate_fit$nuisance_fit_diagnostics,
  artifact_path("rhc_direct_tate_nuisance_fits.csv"),
  row.names = FALSE
)

comparison_diagnostics <- do.call(rbind, lapply(
  names(result$comparison_results),
  function(method_name) {
    fit <- result$comparison_results[[method_name]]
    component <- fit$components
    variance_diagnostics <-
      RoCE:::.comparison_variance_diagnostics(fit)
    data.frame(
      method = method_name,
      estimate = fit$estimate,
      se = fit$se,
      variance_method = component$variance_method,
      se_analytic = component$se_analytic,
      se_bootstrap = component$se_bootstrap,
      n_bootstrap = component$n_bootstrap,
      mu1_reported_se = component$mu1_reported_se,
      mu0_reported_se = component$mu0_reported_se,
      mu1_influence_se = component$mu1_influence_se,
      mu0_influence_se = component$mu0_influence_se,
      cross_arm_covariance = component$cross_arm_covariance,
      cross_arm_correlation = component$cross_arm_correlation,
      dr_weight_n = variance_diagnostics$dr_weight_n,
      dr_weight_n_clipped = variance_diagnostics$dr_weight_n_clipped,
      dr_weight_fraction_clipped =
        variance_diagnostics$dr_weight_fraction_clipped,
      dr_weight_max_site_fraction_clipped =
        variance_diagnostics$dr_weight_max_site_fraction_clipped,
      dr_weight_min_before_clipping =
        variance_diagnostics$dr_weight_min_before_clipping,
      dr_weight_max_before_clipping =
        variance_diagnostics$dr_weight_max_before_clipping,
      stringsAsFactors = FALSE
    )
  }
))
write.csv(
  comparison_diagnostics,
  artifact_path("rhc_direct_tate_comparison_diagnostics.csv"),
  row.names = FALSE
)

forest <- plot_forest_methods(
  result$methods,
  highlight = "RoCE",
  xlab = "Target average treatment effect (risk difference)"
)
forest <- forest + ggplot2::theme(
  axis.title = ggplot2::element_text(size = 21),
  axis.text = ggplot2::element_text(size = 20)
)
save_plot(
  forest,
  artifact_path("rhc_direct_tate_forest.pdf"),
  width = 8.2,
  height = 4.2
)

staged_files <- file.path(staging_directory, artifact_names)
staged_info <- file.info(staged_files)
valid_artifacts <- file.exists(staged_files) &
  !is.na(staged_info$size) & staged_info$size > 0
invalid_artifacts <- staged_files[!valid_artifacts]
if (length(invalid_artifacts) > 0L) {
  stop(
    "RHC staging did not produce every non-empty artifact: ",
    paste(basename(invalid_artifacts), collapse = ", ")
  )
}
if (any(file.exists(output_files))) {
  stop("a final RHC artifact appeared while the run was staging outputs.")
}
for (index in seq_along(staged_files)) {
  if (!file.rename(staged_files[[index]], output_files[[index]])) {
    stop("failed to commit RHC artifact: ", basename(staged_files[[index]]))
  }
  committed_files <- c(committed_files, output_files[[index]])
}
commit_complete <- TRUE

message("[done] RHC TATE outputs written to ", output_root)
}

main()
