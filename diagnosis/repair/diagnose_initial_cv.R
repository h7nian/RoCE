#!/usr/bin/env Rscript
# Capture a failing nuisance-CV task from the unchanged 120/site regression
# fixture, then inspect every CV fold at the original and null-model penalties.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: diagnose_initial_cv.R INSTALLED_LIBRARY NEW_OUTPUT")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
ns <- asNamespace("RoCE")
.captured_initial_cv <- NULL
trace(".call_nuisance_cv_with_groups", where = ns, print = FALSE, tracer = quote({
  if (identical(cv_function, get("select_lambda_cv_initial_density_ratio_cpp", asNamespace("RoCE")))) {
    assign(".captured_initial_cv", list(arguments = args, seed = .Random.seed), globalenv())
  }
}))
set.seed(7102)
data <- split_data_by_site(generate_simulation_data(
  n_total = 360, K = 2, p = 3, config = "C1", dgp_type = "face",
  outcome_type = "continuous", n_target = 120, n_source_sizes = c(120, 120), warn_ignored = FALSE))
failure <- tryCatch(run_tate_crossfit(data, n_folds = 3, communication_mode = "one_round",
  family = "gaussian", nlambda_init = 5, n_cores = 1, verbose = FALSE), error = identity)
untrace(".call_nuisance_cv_with_groups", where = ns)
stopifnot(inherits(failure, "error"), !is.null(.captured_initial_cv))
writeLines(conditionMessage(failure), file.path(output, "original_failure.txt"))
saveRDS(.captured_initial_cv, file.path(output, "cv_inputs.rds"))
inputs <- .captured_initial_cv$arguments
design <- as.matrix(inputs[[1L]])
treatment <- inputs[[2L]]
moment <- inputs[[3L]]
penalties <- inputs[[4L]]
n_folds <- inputs[[5L]]
arm <- inputs[[8L]]
radius <- inputs[[9L]]
selected <- which(treatment == arm)
source_scale <- length(selected) / length(treatment)
assign(".Random.seed", .captured_initial_cv$seed, globalenv())
uniforms <- runif(length(selected))
permutation <- seq_along(selected)
for (index in seq.int(length(selected), 2L)) {
  other <- min(floor(uniforms[index] * index) + 1L, index)
  permutation[c(index, other)] <- permutation[c(other, index)]
}
fold_ids <- ((seq_along(selected) - 1L) %% n_folds + 1L)[permutation]
summary <- list()
for (fold in seq_len(n_folds)) {
  x <- design[selected[fold_ids != fold], , drop = FALSE]
  null_lambda <- max(abs(moment[-1L] - moment[1L] * colMeans(x)))
  write.csv(cbind(intercept = 1, x), file.path(output, paste0("fold_", fold, "_design.csv")), row.names = FALSE)
  write.csv(data.frame(moment = moment), file.path(output, "moment.csv"), row.names = FALSE)
  for (solver in c("coordinate_descent", "proximal_newton")) {
    for (lambda in unique(c(penalties, null_lambda * (1 + 1e-6)))) {
      fit <- RoCE::fit_initial_density_ratio(x, rep(1, nrow(x)), moment / source_scale,
        lambda = lambda / source_scale, M_tau = radius, max_iter = 10000L,
        tol = 1e-9, nuisance_solver = solver)
      linear <- drop(cbind(1, x) %*% fit)
      gradient <- moment - source_scale * colMeans(cbind(1, x) * exp(-pmax(-radius, pmin(radius, linear))))
      slopes <- as.numeric(fit)[-1L]
      residual <- c(gradient[1L], ifelse(abs(slopes) > 1e-8,
        gradient[-1L] + lambda * sign(slopes), pmax(abs(gradient[-1L]) - lambda, 0)))
      summary[[length(summary) + 1L]] <- data.frame(
        fold = fold, solver = solver, lambda = lambda, original_max = max(penalties),
        null_lambda = null_lambda, converged = isTRUE(attr(fit, "converged")),
        iterations = attr(fit, "iterations"), max_abs_coefficient = max(abs(fit)),
        kkt = max(abs(residual)), n_arm_train = nrow(x), source_scale = source_scale)
    }
  }
}
write.csv(do.call(rbind, summary), file.path(output, "fold_fit_summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
