#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) {
  stop(
    paste(
      "usage: probe_nuisance_cv_tail.R OUTPUT.csv CONFIG K SIM_ID",
      "A_VAL SOURCE_INDEX OUTER_FOLD INNER_FOLD"
    ),
    call. = FALSE
  )
}

output_path <- args[[1L]]
configuration <- args[[2L]]
source_count <- as.integer(args[[3L]])
simulation_id <- as.integer(args[[4L]])
treatment_value <- as.integer(args[[5L]])
source_index <- as.integer(args[[6L]])
outer_fold <- as.integer(args[[7L]])
inner_fold <- as.integer(args[[8L]])

if (!configuration %in% c("C1", "C2", "C3") ||
    is.na(source_count) || !source_count %in% c(2L, 4L, 8L) ||
    is.na(simulation_id) || simulation_id < 1L ||
    is.na(treatment_value) || !treatment_value %in% c(0L, 1L) ||
    is.na(source_index) || source_index < 1L ||
    source_index > source_count ||
    is.na(outer_fold) || !outer_fold %in% 1:5 ||
    is.na(inner_fold) || !inner_fold %in% 1:5 ||
    outer_fold == inner_fold) {
  stop("invalid production-path probe arguments.", call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (!nzchar(project_library)) {
  stop("ROCE_PROJECT_LIB must identify the tested package library.",
       call. = FALSE)
}
.libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))

package_provenance <- roce_runtime_package_provenance(project_library)
set.seed(simulation_id)
generated <- RoCE:::generate_simulation_data(
  n_total = 1000L * (source_count + 1L),
  K = source_count,
  p = 100L,
  config = configuration,
  estimand_type = "superpopulation",
  outcome_type = "binary",
  dgp_type = "face",
  ate_deviation = 0,
  n_deviated_sites = 0L,
  warn_ignored = FALSE
)
data_split <- RoCE:::split_data_by_site(generated)
folds <- RoCE:::build_crossfit_folds(data_split, n_folds = 5L)
training_folds <- setdiff(1:5, c(outer_fold, inner_fold))
target_train <- RoCE:::combine_folds(
  folds$target_folds, training_folds
)
source_name <- paste0("s", source_index)
source_train <- RoCE:::combine_folds(
  folds$source_folds[[source_name]], training_folds
)
mean_phi <- c(1, colMeans(target_train$Z_site))

lambda_max <- RoCE:::compute_lambda_max_initial_dr(
  source_train$Z_site, source_train$A, mean_phi, treatment_value
)
lambda_min_ratio <- if (
    nrow(source_train$Z_site) > ncol(source_train$Z_site)) {
  get("LAMBDA_MIN_RATIO_LOW_DIM", envir = asNamespace("RoCE"))
} else {
  get("LAMBDA_MIN_RATIO_HIGH_DIM", envir = asNamespace("RoCE"))
}
lambda_grid <- RoCE:::build_lambda_grid(
  lambda_max = lambda_max,
  lambda_min_ratio = lambda_min_ratio,
  nlambda = 100L
)
support <- RoCE:::.density_ratio_support_penalty_floor(
  source_train$Z_site, source_train$A, mean_phi, treatment_value,
  "probe_nuisance_cv_tail"
)
lambda_grid <- RoCE:::.constrain_density_ratio_lambda_grid(
  lambda_grid, support
)

started_at <- proc.time()[["elapsed"]]
cv <- RoCE:::select_lambda_cv_initial_density_ratio_cpp(
  source_train$Z_site,
  source_train$A,
  mean_phi,
  lambda_grid,
  5L,
  get("MAX_ITER_DEFAULT", envir = asNamespace("RoCE")),
  get("TOL_DEFAULT", envir = asNamespace("RoCE")),
  treatment_value
)
elapsed_seconds <- proc.time()[["elapsed"]] - started_at

descending_order <- order(lambda_grid, decreasing = TRUE)
valid_folds <- as.integer(cv$n_valid_folds)
path_valid_folds <- valid_folds[descending_order]
fully_valid <- path_valid_folds == 5L
first_incomplete <- match(FALSE, fully_valid)
valid_after_first_incomplete <- if (is.na(first_incomplete)) {
  FALSE
} else {
  any(fully_valid[seq.int(first_incomplete, length(fully_valid))])
}
increase_after_first_incomplete <- if (
    is.na(first_incomplete) || first_incomplete == length(path_valid_folds)) {
  FALSE
} else {
  any(diff(path_valid_folds[first_incomplete:length(path_valid_folds)]) > 0L)
}

result <- data.frame(
  path_position = seq_along(descending_order),
  original_lambda_index = descending_order,
  lambda = as.numeric(lambda_grid[descending_order]),
  n_valid_folds = path_valid_folds,
  fully_valid = fully_valid,
  selected_lambda_min =
    descending_order == as.integer(cv$idx_min),
  selected_lambda_1se =
    descending_order == as.integer(cv$idx_1se),
  first_incomplete_path_position = first_incomplete,
  any_fully_valid_after_first_incomplete = valid_after_first_incomplete,
  any_valid_fold_count_increase_after_first_incomplete =
    increase_after_first_incomplete,
  elapsed_seconds = as.numeric(elapsed_seconds),
  configuration = configuration,
  K = source_count,
  p = 100L,
  sim_id = simulation_id,
  A_val = treatment_value,
  source = source_name,
  outer_fold = outer_fold,
  inner_fold = inner_fold,
  support_penalty_floor = support$floor,
  support_constant_feature_count = support$constant_feature_count,
  package_library = package_provenance$library,
  package_fingerprint = package_provenance$fingerprint,
  stringsAsFactors = FALSE
)

output_directory <- dirname(output_path)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile(
  pattern = ".cv_tail_", tmpdir = output_directory, fileext = ".csv"
)
utils::write.csv(result, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path)
  stop("failed to atomically write CV-tail probe output.", call. = FALSE)
}
message(sprintf(
  paste0(
    "[done] %s K=%d sim=%d A=%d %s k1=%d k2=%d: ",
    "first incomplete=%s, valid recovery=%s, count increase=%s, elapsed=%.1fs"
  ),
  configuration, source_count, simulation_id, treatment_value,
  source_name, outer_fold, inner_fold,
  if (is.na(first_incomplete)) "none" else first_incomplete,
  valid_after_first_incomplete, increase_after_first_incomplete,
  elapsed_seconds
))
