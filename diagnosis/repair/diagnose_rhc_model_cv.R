#!/usr/bin/env Rscript
# Reproduce a saved smoke-test failure with its original inputs and RNG stream.
# Expanded lambda ranges here are diagnostics only, not deployed fitting changes.
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
root <- normalizePath(arguments[1L], mustWork = TRUE); output <- arguments[2L]
stopifnot(startsWith(output, paste0(root, "/diagnostics/")), !dir.exists(output))
configuration <- jsonlite::fromJSON(file.path(root, "configuration.json"))
.libPaths(c(configuration$library, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = configuration$library))
source(file.path(root, "workflow/rhc_model_validation.R"))
data <- rhc_generate_population_sample(readRDS(configuration$template), "W_overlap", 929001L)$data
folds <- RoCE:::build_crossfit_folds(data, 10L)
source_data <- RoCE:::combine_folds(folds$source_folds$s2, 1:8)
target_data <- RoCE:::combine_folds(folds$target_folds, 1:8)
moment <- c(1, colMeans(target_data$Z_site)); arm <- 1L
seed <- RoCE:::.nuisance_training_seed(RoCE:::.nuisance_training_key(
  "fit_initial_density_ratio", "s2", 1:8, arm))
n_cv <- RoCE:::.nuisance_cv_fold_count(source_data$A, arm, "diagnosis")
arm_design <- source_data$Z_site[source_data$A == arm, , drop = FALSE]
# Literal R counterpart of src/utils.h create_folds_cpp, verified below against
# the selector's internal RNG path by comparing successful selector outputs.
set.seed(seed)
random <- runif(nrow(arm_design)); indices <- seq_len(nrow(arm_design))
for (i in length(indices):2L) {
  j <- min(floor(random[i] * i) + 1L, i)
  temporary <- indices[i]; indices[i] <- indices[j]; indices[j] <- temporary
}
cv_fold <- ((indices - 1L) %% n_cv) + 1L
full_max <- RoCE:::compute_lambda_max_initial_dr(source_data$Z_site, source_data$A, moment, arm)
fold_max <- vapply(seq_len(n_cv), function(k) max(abs(moment[-1L] -
  colMeans(arm_design[cv_fold != k, , drop = FALSE]))), numeric(1L))
support <- RoCE:::.density_ratio_support_penalty_floor(source_data$Z_site,
  source_data$A, moment, arm, "diagnosis")
dir.create(output, recursive = TRUE)
saveRDS(list(source = source_data, target = target_data, arm = arm, seed = seed,
  cv_fold = cv_fold, full_lambda_max = full_max, cv_lambda_max = fold_max, support = support),
  file.path(output, "inputs.rds"))
print(list(seed = seed, full_lambda_max = full_max, cv_lambda_max = fold_max, support = support))
records <- list()
for (solver in c("proximal_newton", "coordinate_descent")) for (certificate in c(TRUE, FALSE)) {
  set.seed(seed)
  fit <- tryCatch(fit_initial_density_ratio(source_data$Z_site, source_data$A, moment,
    A_val = arm, M_tau = 2, nlambda = 100L, lambda_rule = "min", tol = 1e-10,
    nuisance_solver = solver, nuisance_cv_certificate = certificate), error = function(error) conditionMessage(error))
  saveRDS(fit, file.path(output, paste0(solver, "_", certificate, ".rds")))
  records[[length(records) + 1L]] <- data.frame(solver = solver, certificate = certificate,
    grid = "original", success = !is.character(fit),
    lambda = if (is.character(fit)) NA_real_ else attr(fit, "lambda_used"),
    error = if (is.character(fit)) fit else "")
}
RoCE:::set_nuisance_solver_cpp("proximal_newton")
safe_max <- max(full_max, fold_max) * (1 + 1e-6)
grid <- RoCE:::build_lambda_grid(safe_max, lambda_min_ratio = 1e-4, nlambda = 100L)
grid <- RoCE:::.constrain_density_ratio_lambda_grid(grid, support)
selector <- function(explicit) {
  set.seed(seed)
  RoCE:::select_lambda_cv_initial_density_ratio_cpp(source_data$Z_site, source_data$A,
    moment, grid, n_cv, 10000L, 1e-10, arm, 2,
    cv_fold_id = if (explicit) cv_fold else NULL, use_kkt_certificate = TRUE)
}
expanded <- tryCatch(selector(FALSE), error = function(error) conditionMessage(error))
if (!is.character(expanded)) {
  explicit <- selector(TRUE)
  for (name in intersect(names(expanded), names(explicit)))
    stopifnot(isTRUE(all.equal(expanded[[name]], explicit[[name]], tolerance = 1e-12)))
}
saveRDS(expanded, file.path(output, "expanded_cv_grid.rds"))
records[[length(records) + 1L]] <- data.frame(solver = "proximal_newton", certificate = TRUE,
  grid = "CV_intercept_bound", success = !is.character(expanded),
  lambda = if (is.character(expanded)) NA_real_ else expanded$lambda_min,
  error = if (is.character(expanded)) expanded else "")
write.csv(do.call(rbind, records), file.path(output, "comparison.csv"), row.names = FALSE)
print(do.call(rbind, records))
