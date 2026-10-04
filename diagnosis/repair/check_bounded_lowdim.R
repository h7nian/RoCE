#!/usr/bin/env Rscript
# One low-dimensional development repeat, with oracle benchmarks on identical
# data. This is a staged diagnosis, not evidence of nominal coverage.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 2:6) stop("usage: check_bounded_lowdim.R LIBRARY NEW_OUTPUT [CONFIG] [COVARIATE_SHIFT] [SEED] [DIMENSION]")
.libPaths(c(args[1L], .libPaths()))
library(RoCE)
output <- args[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/"), !file.exists(output))
dir.create(output, recursive = TRUE)
options(error = function() {
  writeLines("FAILED", file.path(output, "status.txt"))
  traceback(15L)
  quit(status = 1L)
})
config <- if (length(args) >= 3L) match.arg(args[3L], c("C1", "C2", "C3")) else "C1"
covariate_shift <- if (length(args) >= 4L) as.numeric(args[4L]) else 1
seed <- if (length(args) >= 5L) as.integer(args[5L]) else 1L
dimension <- if (length(args) >= 6L) as.integer(args[6L]) else 10L
stopifnot(is.finite(seed), seed >= 1L, is.finite(dimension), dimension >= 4L)
configuration <- list(sim_id = seed, K = 2L, p = dimension, config = config, dgp_type = "bounded",
  n_target = 1000L, n_source_sizes = c(1000L, 1000L), dgp_control = list(covariate_shift = covariate_shift),
  n_folds = 10L, nlambda_init = 100L, estimate_ate = TRUE,
  methods = c("one_round_crossfit", "target_only"), return_fitted_tate = TRUE,
  include_quadratic_bias_rule = FALSE, additional_aggregation_modes = c("separate_arms", "joint_tate"),
  target_nuisance_method = "hou_calibrated", source_validation_method = "calibrated",
  calibration_layout = "compact", nuisance_solver = "proximal_newton", nuisance_tol = 1e-10,
  M_tau = 12, M_tau_inference = 12,
  calibration_control = list(recipe = "score_derivative",
    target_propensity_initialization = "calibrated", target_radius = 12),
  n_cores_internal = 1L, parallel_treatment_arms = TRUE, verbose = TRUE)
saveRDS(configuration, file.path(output, "configuration.rds"))
writeLines("RUNNING", file.path(output, "status.txt"))
started <- proc.time()[["elapsed"]]
results <- do.call(run_single_simulation, configuration)
saveRDS(results, file.path(output, "fitted.rds"))
artifacts <- attr(results, "roce_simulation_artifacts")
set.seed(configuration$sim_id)
data <- generate_simulation_data(K = 2L, p = dimension, config = config, dgp_type = "bounded",
  n_target = 1000L, n_source_sizes = c(1000L, 1000L),
  dgp_control = configuration$dgp_control, warn_ignored = FALSE)
stopifnot(identical(split_data_by_site(data), artifacts$data_split))
target <- data$R == "t"
mean_contrast <- data$oracle$mu1[target] - data$oracle$mu0[target]
target_score <- mean_contrast + data$A[target] / data$p_treat_true[target] *
  (data$Y[target] - data$oracle$mu1[target]) - (1 - data$A[target]) /
  (1 - data$p_treat_true[target]) * (data$Y[target] - data$oracle$mu0[target])
oracle <- list(data.frame(method = "target_oracle", estimate = mean(target_score),
  se = sd(target_score) / sqrt(sum(target))))
for (site in c("s1", "s2")) {
  selected <- data$R == site
  residual_contrast <- data$A[selected] * data$oracle$source_weights[selected, "mu1"] *
    (data$Y[selected] - data$oracle$mu1[selected]) - (1 - data$A[selected]) *
    data$oracle$source_weights[selected, "mu0"] * (data$Y[selected] - data$oracle$mu0[selected])
  oracle[[length(oracle) + 1L]] <- data.frame(method = paste0(site, "_oracle_assisted"),
    estimate = mean(mean_contrast) + mean(residual_contrast),
    se = sqrt(var(mean_contrast) / sum(target) + var(residual_contrast) / sum(selected)))
}
oracle <- do.call(rbind, oracle)
oracle$truth <- data$mu1_true - data$mu0_true
oracle$error <- oracle$estimate - oracle$truth
write.csv(oracle, file.path(output, "oracle.csv"), row.names = FALSE)
attr(results, "roce_simulation_artifacts") <- NULL
write.csv(results, file.path(output, "results.csv"), row.names = FALSE)
write.csv(data.frame(elapsed_seconds = proc.time()[["elapsed"]] - started),
  file.path(output, "timing.csv"), row.names = FALSE)
writeLines(c(sprintf("One %s/K2/seed%d repeat, %d direct working features, four active slopes, 1000/site.", config, seed, dimension),
  sprintf("Ten folds, three-level fitting, 100 lambda candidates, common source tilt=%g.", covariate_shift),
  "Source treatment scale=1: even zero common tilt would retain the mixture shape difference.",
  "Radius12 corresponds to population-margin constant6; the structural audit checks its scale.",
  "Aggregation retains the current fixed-cutoff implementation; this does not validate the latest growing-cutoff theorem.",
  "Oracle estimates remain random on a finite dataset. They are compared with fixed integrated truth.",
  "No empirical coverage or RMSE conclusion follows from this single repeat."), file.path(output, "scope.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session_info.txt"))
writeLines("COMPLETE: low-dimensional development repeat only", file.path(output, "status.txt"))
