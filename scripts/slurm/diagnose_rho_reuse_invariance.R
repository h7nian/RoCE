#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/rho_reuse_invariance_v14"
}
project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
provenance <- roce_runtime_package_provenance(project_library)

simulation_seed <- 1L
source_count <- 4L
n_site <- 1000L
dimension <- 100L
n_folds <- 5L
configuration <- "C3"
rho_grid <- c(0, 0.5, 1, 1.5, 2, 2.5)

generate_design <- function(rho) {
  set.seed(simulation_seed)
  data <- RoCE:::generate_simulation_data(
    n_total = n_site * (source_count + 1L),
    K = source_count,
    p = dimension,
    config = configuration,
    estimand_type = "superpopulation",
    outcome_type = "binary",
    dgp_type = "face",
    ate_deviation = rho,
    n_deviated_sites = if (rho > 0) 1L else 0L,
    warn_ignored = FALSE
  )
  split <- RoCE:::split_data_by_site(data)
  folds <- RoCE:::build_crossfit_folds(split, n_folds)
  list(data = data, split = split, folds = folds, rng_state = .Random.seed)
}

designs <- lapply(rho_grid, generate_design)
names(designs) <- format(rho_grid, scientific = FALSE, trim = TRUE)
reference <- designs[[1L]]
fold_indices <- function(folds) {
  list(
    target = lapply(folds$target_folds, `[[`, "original_idx"),
    sources = lapply(folds$source_folds, function(source_folds) {
      lapply(source_folds, `[[`, "original_idx")
    })
  )
}
reference_fold_indices <- fold_indices(reference$folds)

stable_data_fields <- c(
  "X", "X_dagger", "R", "A", "Y_0", "p_treat_true", "Z_site",
  "W_outcome", "mu1_true", "mu0_true"
)
nondeviated_indices <- reference$data$R != "s1"
target_indices <- reference$data$R == "t"
control_indices <- reference$data$A == 0L

rows <- do.call(rbind, lapply(seq_along(designs), function(index) {
  candidate <- designs[[index]]
  rho <- rho_grid[[index]]
  stable_fields_identical <- vapply(
    stable_data_fields,
    function(field) identical(reference$data[[field]], candidate$data[[field]]),
    logical(1L)
  )
  data.frame(
    package_library = provenance$library,
    package_fingerprint = provenance$fingerprint,
    sim_id = simulation_seed,
    config = configuration,
    p = dimension,
    K = source_count,
    rho = rho,
    n_site = n_site,
    n_folds = n_folds,
    stable_fields_exact = all(stable_fields_identical),
    stable_field_failures = paste(
      names(stable_fields_identical)[!stable_fields_identical],
      collapse = ";"
    ),
    rng_state_exact = identical(reference$rng_state, candidate$rng_state),
    folds_exact = identical(
      reference_fold_indices, fold_indices(candidate$folds)
    ),
    target_y1_exact = identical(
      reference$data$Y_1[target_indices], candidate$data$Y_1[target_indices]
    ),
    nondeviated_y1_exact = identical(
      reference$data$Y_1[nondeviated_indices],
      candidate$data$Y_1[nondeviated_indices]
    ),
    control_observed_y_exact = identical(
      reference$data$Y[control_indices], candidate$data$Y[control_indices]
    ),
    deviated_source_y1_changed = if (rho == 0) {
      FALSE
    } else {
      any(
        reference$data$Y_1[reference$data$R == "s1"] !=
          candidate$data$Y_1[candidate$data$R == "s1"]
      )
    },
    stringsAsFactors = FALSE
  )
}))

required_invariance <- with(rows,
  stable_fields_exact & rng_state_exact & folds_exact & target_y1_exact &
    nondeviated_y1_exact & control_observed_y_exact &
    (rho == 0 | deviated_source_y1_changed)
)
rows$required_invariance_passed <- required_invariance

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_directory, "rho_reuse_invariance.csv")
sentinel_path <- file.path(output_directory, "rho_reuse_invariance_passed.txt")
if (file.exists(output_path) || file.exists(sentinel_path)) {
  stop("refusing to overwrite an existing rho-reuse audit.", call. = FALSE)
}
temporary_output <- tempfile(
  pattern = ".rho_reuse_invariance_", tmpdir = output_directory,
  fileext = ".csv"
)
write.csv(rows, temporary_output, row.names = FALSE)
if (!file.rename(temporary_output, output_path)) {
  unlink(temporary_output)
  stop("failed to atomically write rho-reuse diagnostics.", call. = FALSE)
}
if (!all(required_invariance)) {
  stop(
    "rho-reuse invariance failed; do not reuse fitted components across rho.",
    call. = FALSE
  )
}
temporary_sentinel <- tempfile(
  pattern = ".rho_reuse_invariance_passed_", tmpdir = output_directory,
  fileext = ".tmp"
)
writeLines(
  c(
    "rho_reuse_invariance=passed",
    paste0("package_fingerprint=", provenance$fingerprint),
    "scope=C3,p100,K4,seed1,rho0_to_2p5"
  ),
  temporary_sentinel
)
if (!file.rename(temporary_sentinel, sentinel_path)) {
  unlink(temporary_sentinel)
  stop("failed to atomically write rho-reuse sentinel.", call. = FALSE)
}
print(rows, row.names = FALSE)
message("[pass] DGP, RNG, folds, control arm, target, and nondeviated sources are rho-invariant.")
