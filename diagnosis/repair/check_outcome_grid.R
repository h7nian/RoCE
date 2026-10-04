#!/usr/bin/env Rscript
# Profile exactly the cached task or validate instrumentation on its path prefix.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("usage: check_outcome_grid.R CHECK_ROOT INPUT_RDS NEW_OUTPUT [PREFIX_LENGTH]")
}
check <- normalizePath(args[1L])
.libPaths(c(file.path(check, "Rlib"), .libPaths()))
library(RoCE)
saved <- readRDS(args[2L])
arguments <- saved$arguments
output <- args[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
writeLines("RUNNING", file.path(output, "status.txt"))
options(error = function() {
  writeLines("FAILED", file.path(output, "status.txt"))
  traceback(15L)
  quit(status = 1L)
})
driver <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
cache <- Sys.getenv("ROCE_DIAGNOSTIC_CPP_CACHE", file.path(output, "compiled"))
dir.create(cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(file.path(check, "source", "src"))),
           PKG_CXXFLAGS = "-g0 -O2", ROCE_NUISANCE_SOLVER = "proximal_newton")
Rcpp::sourceCpp(file.path(dirname(normalizePath(driver)), "profile_outcome_grid.cpp"), cacheDir = cache)
RoCE:::set_nuisance_solver_cpp("proximal_newton")
lambda_max <- RoCE:::compute_lambda_max_outcome(arguments$W_outcome, arguments$Y,
  arguments$A, arguments$gamma_s, A_val = arguments$A_val,
  family_int = arguments$family_int, link_int = arguments$link_int,
  Z_site = arguments$Z_site, calibrated = TRUE, M_tau = arguments$M_tau,
  use_weight_derivative = arguments$use_weight_derivative)
ratio <- if (nrow(arguments$W_outcome) > ncol(arguments$W_outcome)) 1e-4 else 0.01
grid <- RoCE:::build_lambda_grid(lambda_max, ratio, arguments$nlambda)
if (length(args) == 4L) {
  prefix <- as.integer(args[4L])
  stopifnot(prefix >= 2L, prefix <= length(grid))
  grid <- grid[seq_len(prefix)]
}
saveRDS(list(input = normalizePath(args[2L]), input_hash = saved$input_hash,
  grid = grid, cv_seed = saved$cv_seed), file.path(output, "configuration.rds"))
n_folds <- RoCE:::.nuisance_cv_fold_count(arguments$A, arguments$A_val, "outcome profile")
set.seed(saved$cv_seed)
started <- proc.time()[["elapsed"]]
profile <- profile_outcome_grid_cpp(arguments$W_outcome, arguments$Y, arguments$A,
  arguments$gamma_s, arguments$Z_site, grid, arguments$A_val, n_folds,
  arguments$M_tau, arguments$tol, arguments$max_iter,
  arguments$family_int, arguments$link_int, arguments$use_weight_derivative,
  file.path(output, "fold_lambda_profile.csv"))
elapsed <- proc.time()[["elapsed"]] - started
saveRDS(profile, file.path(output, "profile.rds"))
if (length(args) == 4L) {
  native <- RoCE:::select_lambda_cv_calibrated_outcome_cpp(arguments$W_outcome,
    arguments$Y, arguments$A, arguments$gamma_s, grid, n_folds, arguments$max_iter,
    arguments$tol, arguments$A_val, arguments$M_tau, arguments$Z_site,
    arguments$family_int, arguments$link_int, cv_fold_id = profile$cv_fold_id,
    use_weight_derivative = arguments$use_weight_derivative)
  stopifnot(isTRUE(all.equal(profile$summary$cv_scores, native$cv_scores, tolerance = 1e-12)),
            identical(profile$summary$lambda_min, native$lambda_min))
} else {
  stopifnot(identical(profile$summary$lambda_min, attr(saved$reference_fit, "lambda_used")),
            identical(profile$summary$lambda_1se, attr(saved$reference_fit, "lambda_1se")),
            profile$summary$invalid_lambdas == attr(saved$reference_fit, "cv_invalid_lambdas"),
            profile$summary$path_tail_skipped_fold_fits == attr(saved$reference_fit, "cv_path_tail_skipped_fold_fits"))
}
weights <- profile$arm_weights
write.csv(data.frame(elapsed_seconds = elapsed, arm_observations = length(weights),
  positive_weight_observations = sum(weights > 0), effective_weighted_n = sum(weights)^2 / sum(weights^2),
  selected_lambda = profile$summary$lambda_min), file.path(output, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("PASSED", file.path(output, "status.txt"))
