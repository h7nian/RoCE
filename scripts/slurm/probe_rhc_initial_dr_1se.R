#!/usr/bin/env Rscript

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))

output_root <- Sys.getenv(
  "ROCE_OUTPUT_ROOT", "results/direct_tate_v1/rhc_initial_dr_1se_probe"
)
lambda_rule <- Sys.getenv("ROCE_NUISANCE_LAMBDA_RULE", "1se")
if (!lambda_rule %in% c("min", "1se")) {
  stop("ROCE_NUISANCE_LAMBDA_RULE must be 'min' or '1se'.")
}
output_basename <- Sys.getenv(
  "ROCE_RHC_PROBE_FILENAME", "rhc_initial_dr_1se_probe.csv"
)
if (!nzchar(output_basename) || basename(output_basename) != output_basename) {
  stop("ROCE_RHC_PROBE_FILENAME must be one non-empty file basename.")
}
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(output_root, output_basename)
if (file.exists(output_path)) {
  stop("Probe output already exists; use a fresh ROCE_OUTPUT_ROOT: ", output_path)
}

cohort <- build_rhc_cohort(outcome = "death30", site_var = "ninsclas")
data_split <- build_rhc_data_split(
  cohort = cohort, K = 5L, target_site = "Private", seed = 42L
)
n_folds <- 10L
folds <- RoCE:::build_crossfit_folds(data_split, n_folds)
target_folds <- folds$target_folds
source_folds <- folds$source_folds

# These are the only outer-fold/source combinations that failed under
# lambda.min in the archived primary fit.  All other arm/source paths
# converged and are intentionally excluded from this focused gate.
probe_cases <- data.frame(
  outer_fold = c(5L, 9L, 10L, 10L),
  source = c("s4", "s4", "s1", "s4"),
  stringsAsFactors = FALSE
)

rows <- list()
for (case_index in seq_len(nrow(probe_cases))) {
  outer_fold <- probe_cases$outer_fold[[case_index]]
  source <- probe_cases$source[[case_index]]
  secondary_folds <- setdiff(seq_len(n_folds), outer_fold)
  cached_lambda <- NULL
  previous_gamma <- NULL

  for (inner_fold in secondary_folds) {
    training_folds <- setdiff(seq_len(n_folds), c(outer_fold, inner_fold))
    target_train <- RoCE:::combine_folds(target_folds, training_folds)
    source_train <- RoCE:::combine_folds(
      source_folds[[source]], training_folds
    )
    mean_phi <- c(1, colMeans(target_train$Z_site))

    gamma <- fit_initial_density_ratio(
      source_train$Z_site,
      source_train$A,
      mean_phi,
      lambda = cached_lambda,
      max_iter = 10000L,
      A_val = 0L,
      warm_start = previous_gamma,
      nlambda = 100L,
      lambda_rule = lambda_rule
    )
    if (is.null(cached_lambda)) {
      cached_lambda <- attr(gamma, "lambda_used")
    }
    previous_gamma <- RoCE:::.safe_density_ratio_warm_start(gamma)

    rows[[length(rows) + 1L]] <- data.frame(
      outer_fold = outer_fold,
      source = source,
      inner_fold = inner_fold,
      n_source_train = source_train$n,
      n_source_control = sum(source_train$A == 0L),
      lambda_used = as.numeric(attr(gamma, "lambda_used")),
      lambda_min = as.numeric(attr(gamma, "lambda_min")),
      lambda_1se = as.numeric(attr(gamma, "lambda_1se")),
      lambda_selected_on_path = as.numeric(
        attr(gamma, "lambda_selected_on_path")
      ),
      lambda_selected_before_support_floor = as.numeric(
        attr(gamma, "lambda_selected_before_support_floor")
      ),
      support_penalty_floor = as.numeric(
        attr(gamma, "support_penalty_floor")
      ),
      support_penalty_floor_applied =
        isTRUE(attr(gamma, "support_penalty_floor_applied")),
      support_constant_feature_count = as.integer(
        attr(gamma, "support_constant_feature_count")
      ),
      converged = isTRUE(attr(gamma, "converged")),
      iterations = as.integer(attr(gamma, "iterations")),
      update_to_threshold_ratio = as.numeric(
        attr(gamma, "update_to_threshold_ratio")
      ),
      line_search_failures = as.integer(attr(gamma, "line_search_failures")),
      max_abs_coefficient = max(abs(as.numeric(gamma))),
      stringsAsFactors = FALSE
    )
  }
}

result <- do.call(rbind, rows)
temporary_path <- tempfile(
  pattern = ".rhc_initial_dr_probe_", tmpdir = output_root, fileext = ".csv"
)
write.csv(result, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path)
  stop("Failed to atomically install probe output: ", output_path)
}

cat("Probe fits:", nrow(result), "\n")
cat("Lambda rule:", lambda_rule, "\n")
cat("Nonconverged:", sum(!result$converged), "\n")
cat("Line-search failures:", sum(result$line_search_failures), "\n")
cat("Support-floor activations:",
    sum(result$support_penalty_floor_applied), "\n")
cat("Maximum |gamma|:", max(result$max_abs_coefficient), "\n")
print(result, row.names = FALSE)
if (any(!result$converged) || any(result$line_search_failures > 0L) ||
    any(!is.finite(result$max_abs_coefficient)) ||
    any(result$max_abs_coefficient >= 0.9 * RoCE:::PARAM_MAX)) {
  stop("The RHC initial-density-ratio probe failed stability checks.")
}
