#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0006 Land the common-basis DGP with X-dagger misspecification
# Task:    Reproduce the #0002 pilot cells with the in-package DGP: for C1 and for
#          C2/C3/C4 at the pre-registered strength, generate seed 20001 with
#          generate_simulation_data(), compare standardization, calibration and truth
#          with the prototype's stored fit bundle, refit with run_tate_crossfit() and
#          compare estimate, fixed-weight SE and fold weights. Array task = one
#          configuration; writes diagnosis/out/dgp_common_basis/package_<config>/.

suppressPackageStartupMessages(library(RoCE))
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

P <- 100L
N_TARGET <- 1000L
N_SOURCE <- c(1000L, 1000L)
SIM_ID <- 20001L
CONFIGS <- c("C1", "C2", "C3", "C4")
CELL_ROOT <- "diagnosis/out/dgp_common_basis"

.prototype_cell <- function(config) {
  strength <- if (config == "C1") 0 else RoCE:::FACE_MISSPECIFICATION_STRENGTH
  directory <- file.path(CELL_ROOT, sprintf("%s_omega%.2f", config, strength))
  list(
    bundle = readRDS(file.path(directory, "fit.rds")),
    summary = read.csv(file.path(directory, "summary.csv"))
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: dgp_common_basis.R OUTPUT_DIRECTORY CONFIG_INDEX")
  config <- CONFIGS[as.integer(args[2])]
  output <- file.path(args[1], paste0("package_", config))
  if (file.exists(output)) stop("output already exists: ", output)
  prototype <- .prototype_cell(config)

  set.seed(SIM_ID)
  data <- generate_simulation_data(
    K = length(N_SOURCE), p = P, config = config, dgp_type = "face",
    outcome_type = "binary", estimand_type = "superpopulation",
    n_target = N_TARGET, n_source_sizes = N_SOURCE, warn_ignored = FALSE
  )
  reference <- RoCE:::.face_reference_population(
    P, RoCE:::FACE_KAPPA, config, data$misspecification_strength, "binary"
  )
  difference <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
  dgp_checks <- data.frame(
    config = config,
    misspecification_strength = data$misspecification_strength,
    standardization_difference = max(
      difference(reference$standardization$mean, prototype$bundle$standardization$mean),
      difference(reference$standardization$sd, prototype$bundle$standardization$sd)
    ),
    calibration_difference = max(
      difference(reference$calibration$eta_mean, prototype$bundle$calibration$eta_mean),
      difference(reference$calibration$eta_sd, prototype$bundle$calibration$eta_sd)
    ),
    truth_difference = difference(data$mu1_true - data$mu0_true, prototype$summary$truth)
  )
  print(t(dgp_checks))

  fit <- run_tate_crossfit(
    split_data_by_site(data), n_folds = 5L, communication_mode = "one_round",
    lambda_selection = RoCE:::AGG_WALD_LAMBDA, verbose = FALSE,
    M_tau = 5, M_tau_inference = 5, n_cores = length(N_SOURCE),
    nlambda_init = 100L, family = "binomial", nuisance_lambda_rule = "min"
  )
  stored <- prototype$bundle$fit
  fit_checks <- data.frame(
    config = config,
    estimate_difference = difference(fit$estimate, stored$estimate),
    fixed_se_difference = difference(fit$se_fixed_weights, stored$se),
    fold_weight_difference = difference(fit$fold_weights, stored$fold_weights),
    target_only_difference = difference(fit$target_only$estimate, stored$target_only$estimate),
    weight_layer_se = fit$se,
    fixed_weight_se = fit$se_fixed_weights
  )
  print(t(fit_checks))

  roce_write_atomic_directory(output, function(stage) {
    write.csv(dgp_checks, file.path(stage, "dgp_checks.csv"), row.names = FALSE)
    write.csv(fit_checks, file.path(stage, "fit_checks.csv"), row.names = FALSE)
    writeLines(c("task=dgp_common_basis", "history_entry=0006",
                 paste0("config=", config), paste0("sim_id=", SIM_ID),
                 paste0("package_library=", .libPaths()[1])),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "common basis DGP package reproduction")
}

if (sys.nframe() == 0L) main()
