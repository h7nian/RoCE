#!/usr/bin/env Rscript
# Bounded p100 check of one deepest-layer initial density-ratio fit. This does
# not estimate TATE or coverage; it verifies the numerical path before a full
# three-level study is scheduled.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: check_large_basis.R INSTALLED_LIBRARY NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
set.seed(1081)
data <- split_data_by_site(generate_simulation_data(
  n_total = 3000L, n_target = 1000L, n_source_sizes = c(1000L, 1000L),
  K = 2L, p = 100L, config = "C3", dgp_type = "face", outcome_type = "binary",
  warn_ignored = FALSE))
folds <- RoCE:::build_crossfit_folds(data, 10L)
source_train <- RoCE:::combine_folds(folds$source_folds$s1, 4:10)
target_train <- RoCE:::combine_folds(folds$target_folds, 4:10)
arguments <- list(Z_site = source_train$Z_site, A = source_train$A,
  mean_phi = c(1, colMeans(target_train$Z_site)), A_val = 1L,
  nlambda = 100L, M_tau = 5, tol = 1e-10)
saveRDS(list(data = data, folds = folds, arguments = arguments), file.path(output, "inputs.rds"))
writeLines(c("n=1000/site; raw p=100; 200-column quadratic working basis; ten original folds",
  "initial fit excludes outer 1, validation 2, calibration 3",
  "initial transport nuisance only; no TATE/coverage claim",
  "prespecified checks: lambda difference <=1e-10; prediction difference <=1e-6; KKT <=1e-10"),
  file.path(output, "scope.txt"))
fits <- list()
summary <- list()
for (solver in c("coordinate_descent", "proximal_newton")) {
  cat("Starting", solver, "\n")
  started <- proc.time()[["elapsed"]]
  fit <- tryCatch(RoCE:::.fit_nuisance_training_subset("fit_initial_density_ratio",
    c(arguments, list(nuisance_solver = solver)), "s1", 4:10), error = identity)
  saveRDS(fit, file.path(output, paste0(solver, ".rds")))
  if (inherits(fit, "error")) {
    writeLines(conditionMessage(fit), file.path(output, paste0(solver, "_FAILED.txt")))
    stop(conditionMessage(fit))
  }
  fits[[solver]] <- fit
  summary[[solver]] <- data.frame(solver = solver, n_train = nrow(source_train$Z_site),
    n_arm = sum(source_train$A == 1), basis_columns = ncol(source_train$Z_site),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    lambda = attr(fit, "lambda_used"), kkt = attr(fit, "kkt_residual"),
    invalid_cv_lambdas = attr(fit, "cv_invalid_lambdas"))
  write.csv(do.call(rbind, summary), file.path(output, "fits.csv"), row.names = FALSE)
}
delta <- as.numeric(fits[[1L]]) - as.numeric(fits[[2L]])
predictor_difference <- max(vapply(data, function(site) {
  max(abs(drop(cbind(1, site$Z_site) %*% delta)))
}, numeric(1L)))
lambda_difference <- abs(attr(fits[[1L]], "lambda_used") - attr(fits[[2L]], "lambda_used"))
comparison <- data.frame(max_coefficient_difference = max(abs(delta)),
  max_predictor_difference = predictor_difference, lambda_difference = lambda_difference)
write.csv(comparison, file.path(output, "comparison.csv"), row.names = FALSE)
passed <- predictor_difference <= 1e-6 && lambda_difference <= 1e-10 &&
  all(vapply(fits, function(fit) attr(fit, "kkt_residual") <= 1e-10, logical(1L)))
writeLines(if (passed) "PASSED" else "FAILED", file.path(output, "status.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
if (!passed) stop("Large-basis equivalence criteria failed; results retained.")
