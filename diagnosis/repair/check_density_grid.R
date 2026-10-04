#!/usr/bin/env Rscript
# Isolate grid length/range from upstream fitting on saved p100 data.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(3L, 4L)) {
  stop("usage: check_density_grid.R CHECK_ROOT INPUT_RDS NEW_OUTPUT [comma-separated variants]")
}
check_root <- normalizePath(args[1L])
.libPaths(c(file.path(check_root, "Rlib"), .libPaths()))
library(RoCE)
output <- args[3L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
compiled_cache <- Sys.getenv("ROCE_DIAGNOSTIC_CPP_CACHE", file.path(output, "compiled"))
dir.create(compiled_cache, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PKG_CPPFLAGS = paste("-I", shQuote(file.path(check_root, "source", "src"))),
           PKG_CXXFLAGS = "-g0 -O2", ROCE_NUISANCE_SOLVER = "proximal_newton")
Rcpp::sourceCpp("diagnosis/repair/profile_density_grid.cpp", cacheDir = compiled_cache)
saved <- readRDS(args[2L])
arguments <- saved$arguments
key <- RoCE:::.nuisance_training_key("fit_initial_density_ratio", "s1", 4:10, 1L)
seed <- RoCE:::.nuisance_training_seed(key)
fold_count <- RoCE:::.nuisance_cv_fold_count(arguments$A, 1L, "grid diagnosis")
maximum <- RoCE:::compute_lambda_max_initial_dr(arguments$Z_site, arguments$A, arguments$mean_phi, 1L)
grids <- list(grid100 = RoCE:::build_lambda_grid(maximum, 1e-4, 100L),
              grid50 = RoCE:::build_lambda_grid(maximum, 1e-4, 50L),
              grid25 = RoCE:::build_lambda_grid(maximum, 1e-4, 25L))
grids$prefix26 <- grids$grid100[1:26]
saveRDS(list(arguments = arguments, grids = grids, seed = seed), file.path(output, "inputs.rds"))
writeLines(c("One initial source-weight task: C3, 1000/site, raw p100 / 200 features, ten original folds",
  "Training folds 4:10: outer/validation/calibration folds 1/2/3 excluded",
  "Same data, objective and CV folds across grids; no upstream nuisance refits",
  "This component study is not TATE coverage or approval to reduce the production grid"),
  file.path(output, "scope.txt"))
invisible(RoCE:::set_nuisance_solver_cpp("proximal_newton"))
fold_ids <- NULL
summaries <- fits <- list()
variants <- if (length(args) == 4L) strsplit(args[4L], ",", fixed = TRUE)[[1L]] else
  c(names(grids), "exhaustive100", "centered100")
stopifnot(variants[1L] == "grid100", !anyDuplicated(variants),
          all(variants %in% c(names(grids), "exhaustive100", "centered100")))
for (name in variants) {
  grid <- if (name %in% c("exhaustive100", "centered100")) grids$grid100 else grids[[name]]
  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  result <- profile_density_grid_cpp(arguments$Z_site, arguments$A, arguments$mean_phi,
    rep(1, length(arguments$A)), grid, 1L, fold_count, 5, 1e-10, 10000L,
    exhaustive = name == "exhaustive100", center_design = name == "centered100",
    cv_fold_id = fold_ids)
  elapsed <- proc.time()[["elapsed"]] - started
  if (is.null(fold_ids)) fold_ids <- result$cv_fold_id
  saveRDS(result, file.path(output, paste0(name, ".rds")))
  write.csv(result$profile, file.path(output, paste0(name, "_profile.csv")), row.names = FALSE)
  if (nzchar(result$error)) stop(name, ": ", result$error)
  native <- RoCE:::select_lambda_cv_initial_density_ratio_cpp(arguments$Z_site, arguments$A,
    arguments$mean_phi, grid, fold_count, 10000L, 1e-10, 1L, 5, cv_fold_id = fold_ids)
  if (!name %in% c("exhaustive100", "centered100")) {
    stopifnot(isTRUE(all.equal(result$summary$cv_scores, native$cv_scores, tolerance = 1e-12)),
              identical(result$summary$lambda_min, native$lambda_min))
  }
  fit <- tryCatch(RoCE::fit_initial_density_ratio(arguments$Z_site, arguments$A, arguments$mean_phi,
    lambda = result$summary$lambda_min, M_tau = 5, tol = 1e-10), error = identity)
  refit_ok <- !inherits(fit, "error") && isTRUE(attr(fit, "converged"))
  if (name != "centered100") stopifnot(refit_ok)
  saveRDS(fit, file.path(output, paste0(name, "_refit.rds")))
  fits[[name]] <- fit
  profile <- result$profile
  summaries[[name]] <- data.frame(grid = name, points = length(grid),
    lambda_max = max(grid), lambda_min = min(grid), selected_lambda = result$summary$lambda_min,
    eligible_lambdas = sum(is.finite(result$summary$cv_scores)),
    attempted_fold_fits = sum(profile$attempted), skipped_fold_fits = sum(!profile$attempted),
    failed_attempts = sum(profile$attempted & !profile$converged),
    converged_seconds = sum(profile$seconds[profile$converged]),
    failed_seconds = sum(profile$seconds[profile$attempted & !profile$converged]),
    elapsed_seconds = elapsed,
    raw_refit_converged = refit_ok,
    coefficient_difference = if (refit_ok)
      max(abs(as.numeric(fit) - as.numeric(fits$grid100))) else NA_real_)
  write.csv(do.call(rbind, summaries), file.path(output, "comparison.csv"), row.names = FALSE)
  cat(name, "finished:", elapsed, "seconds; selected lambda", result$summary$lambda_min, "\n")
}
saveRDS(fits, file.path(output, "selected_fits.rds"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("COMPLETED; compare grids before accepting any change", file.path(output, "status.txt"))
