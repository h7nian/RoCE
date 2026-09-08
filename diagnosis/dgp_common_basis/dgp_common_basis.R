#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0002 Common working basis and X-dagger misspecification: strength pilot
# Task:    One seed (20001), p = 100, K = 2, 1000 per site, rho = 0. Both working
#          bases equal [X - kappa, X^2]; C2/C3/C4 misspecify the outcome / propensity
#          / both mechanisms through eta_mix = (1 - omega) eta(X) + omega eta(X_dagger).
#          Array task = one (config, omega) cell; writes
#          diagnosis/out/dgp_common_basis/<config>_omega<omega>/.

suppressPackageStartupMessages(library(RoCE))
`%||%` <- RoCE:::`%||%`
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

P <- 100L
N_TARGET <- 1000L
N_SOURCE <- c(1000L, 1000L)
SIM_ID <- 20001L
N_REF <- 100000L
REF_SEED <- 99999L
OMEGA_GRID <- c(0.25, 0.5, 0.75, 1)
CELLS <- rbind(
  data.frame(config = "C1", omega = 0),
  expand.grid(config = c("C2", "C3", "C4"), omega = OMEGA_GRID, stringsAsFactors = FALSE)
)
CURRENT_C1_TATE_SEED20001 <- 0.194329335788676

# Kang-Schafer-style transforms of the four signal coordinates, standardized on
# the target reference population to mean kappa and unit sd so that the same
# coefficient structure applies to X and X_dagger.
.transform_signal_coordinates <- function(X) {
  cbind(exp(X[, 1] / 2),
        X[, 2] / (1 + exp(X[, 1])) + 10,
        (X[, 1] * X[, 3] / 25 + 0.6)^3,
        (X[, 2] + X[, 4] + 20)^2)
}

.make_x_dagger <- function(X, standardization, kappa) {
  transformed <- .transform_signal_coordinates(X)
  X_dagger <- X
  for (j in 1:4) {
    X_dagger[, j] <- kappa + (transformed[, j] - standardization$mean[j]) / standardization$sd[j]
  }
  X_dagger
}

.linear_quadratic_predictor <- function(X, linear, squared, kappa) {
  as.numeric(sweep(X, 2L, kappa, "-") %*% linear + (X^2) %*% squared)
}

.mixed_predictor <- function(X, X_dagger, linear, squared, kappa, omega) {
  (1 - omega) * .linear_quadratic_predictor(X, linear, squared, kappa) +
    omega * .linear_quadratic_predictor(X_dagger, linear, squared, kappa)
}

.working_basis <- function(X, kappa) {
  cbind(sweep(X, 2L, kappa, "-"), X^2)
}

# Reference-population quantities: standardization constants, binary outcome
# calibration of the (possibly mixed) predictor, truth, and the R^2 of each
# true predictor on the working basis.
.reference_population <- function(omega_outcome, omega_ps, kappa) {
  out_params <- get_face_outcome_parameters(P)
  ps_params <- get_face_ps_parameters(P)
  RoCE:::with_seed(REF_SEED, {
    X_ref <- RoCE:::generate_face_covariates(N_REF, integer(0), P, kappa = kappa)$X
    transformed <- .transform_signal_coordinates(X_ref)
    standardization <- list(mean = colMeans(transformed), sd = apply(transformed, 2L, stats::sd))
    X_dagger_ref <- .make_x_dagger(X_ref, standardization, kappa)
    eta_outcome <- .mixed_predictor(X_ref, X_dagger_ref, out_params$beta_linear,
                                    out_params$beta_squared, kappa, omega_outcome)
    eta_ps <- .mixed_predictor(X_ref, X_dagger_ref, ps_params$alpha1, ps_params$alpha2, 0, omega_ps)
    calibration <- list(eta_mean = mean(eta_outcome),
                        eta_sd = max(stats::sd(eta_outcome), RoCE:::EPSILON_DEFAULT),
                        signal_sd = RoCE:::FACE_BINARY_SIGNAL_SD)
    g_ref <- RoCE:::face_binary_logit(eta_outcome, calibration)
    basis <- .working_basis(X_ref, kappa)
    r_squared <- function(y) summary(stats::lm(y ~ basis))$r.squared
    p_ref <- RoCE:::logistic(eta_ps)
    list(
      standardization = standardization,
      calibration = calibration,
      mu0_true = mean(RoCE:::logistic(g_ref)),
      mu1_true = mean(RoCE:::logistic(g_ref + RoCE:::FACE_BINARY_ATE_TARGET)),
      r_squared_outcome = r_squared(eta_outcome),
      r_squared_propensity = r_squared(eta_ps),
      propensity_clip_fraction = mean(p_ref < RoCE:::POSITIVITY_LOWER | p_ref > RoCE:::POSITIVITY_UPPER),
      propensity_quantiles = stats::quantile(p_ref, c(0.001, 0.01, 0.5, 0.99, 0.999))
    )
  })
}

# Mirrors generate_face_data()'s random-number order (covariates, treatment,
# Y(1), Y(0)) so that C1 reproduces the current DGP up to basis centering.
.generate_common_basis_data <- function(config, omega, reference, kappa) {
  omega_outcome <- if (config %in% c("C2", "C4")) omega else 0
  omega_ps <- if (config %in% c("C3", "C4")) omega else 0
  out_params <- get_face_outcome_parameters(P)
  ps_params <- get_face_ps_parameters(P)
  set.seed(SIM_ID)
  covariates <- RoCE:::generate_face_covariates(N_TARGET, N_SOURCE, P, kappa = kappa)
  X <- covariates$X
  R <- covariates$R
  X_dagger <- .make_x_dagger(X, reference$standardization, kappa)
  eta_ps <- .mixed_predictor(X, X_dagger, ps_params$alpha1, ps_params$alpha2, 0, omega_ps)
  p_treat <- pmax(pmin(RoCE:::logistic(eta_ps), RoCE:::POSITIVITY_UPPER), RoCE:::POSITIVITY_LOWER)
  A <- stats::rbinom(length(R), 1L, p_treat)
  ate_map <- RoCE:::build_face_ate_map(length(N_SOURCE), 0, 0L, base_ate = RoCE:::FACE_BINARY_ATE_TARGET)
  eta_outcome <- .mixed_predictor(X, X_dagger, out_params$beta_linear,
                                  out_params$beta_squared, kappa, omega_outcome)
  g <- RoCE:::face_binary_logit(eta_outcome, reference$calibration)
  draw_outcome <- function(a) stats::rbinom(length(R), 1L, RoCE:::logistic(g + ate_map[R] * a))
  Y_1 <- draw_outcome(1L)
  Y_0 <- draw_outcome(0L)
  basis <- .working_basis(X, kappa)
  list(
    n = length(R), K = length(N_SOURCE), p = P, X = X, X_dagger = X_dagger, R = R, A = A,
    Y = ifelse(A == 1L, Y_1, Y_0), Y_1 = Y_1, Y_0 = Y_0, p_treat_true = p_treat,
    mu1_true = reference$mu1_true, mu0_true = reference$mu0_true,
    mu1_realized = mean(Y_1[R == "t"]), mu0_realized = mean(Y_0[R == "t"]),
    gamma_params = NULL, alpha1_true = NULL, alpha0_true = NULL, out_params = out_params,
    ate_map = ate_map, treatment_shift_map = ate_map, config = config,
    Z_site = basis, W_outcome = basis, estimand_type = "superpopulation",
    outcome_type = "binary", dgp_type = "face"
  )
}

.fit_summary <- function(fit, config, omega, data, reference) {
  clip <- fit$clip_diagnostics
  diagnostics <- fit$nuisance_fit_diagnostics
  nonconverged <- if (is.list(diagnostics)) {
    flags <- unlist(diagnostics[grepl("converged", names(diagnostics))])
    if (length(flags)) sum(!as.logical(flags)) else NA_integer_
  } else NA_integer_
  data.frame(
    config = config, omega = omega, sim_id = SIM_ID,
    tate_estimate = fit$estimate, tate_se = fit$se,
    target_only_estimate = fit$target_only$estimate, target_only_se = fit$target_only$se,
    truth = data$mu1_true - data$mu0_true,
    weight_s1 = fit$weights[1], weight_s2 = fit$weights[2],
    wald_s1_mean = mean(fit$fold_wald_statistics[, 1]),
    wald_s2_mean = mean(fit$fold_wald_statistics[, 2]),
    inference_logit_truncation_fraction = clip$logit_truncation_fraction %||% NA_real_,
    max_raw_weight = clip$max_raw_weight %||% NA_real_,
    max_abs_logit = clip$max_abs_logit %||% NA_real_,
    nonconverged_fits = nonconverged,
    r_squared_outcome = reference$r_squared_outcome,
    r_squared_propensity = reference$r_squared_propensity,
    propensity_clip_fraction = reference$propensity_clip_fraction,
    c1_reference_difference = if (config == "C1") fit$estimate - CURRENT_C1_TATE_SEED20001 else NA_real_
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: dgp_common_basis.R OUTPUT_DIRECTORY CELL_INDEX")
  cell <- CELLS[as.integer(args[2]), ]
  config <- as.character(cell$config)
  omega <- as.numeric(cell$omega)
  output <- file.path(args[1], sprintf("%s_omega%.2f", config, omega))
  if (file.exists(output)) stop("output already exists: ", output)
  kappa <- RoCE:::FACE_KAPPA
  reference <- .reference_population(
    omega_outcome = if (config %in% c("C2", "C4")) omega else 0,
    omega_ps = if (config %in% c("C3", "C4")) omega else 0,
    kappa = kappa
  )
  if (config == "C1") {
    package_calibration <- RoCE:::get_face_binary_calibration(P, kappa = kappa)
    stopifnot(abs(reference$mu1_true - package_calibration$mu1_superpop) < 1e-12,
              abs(reference$mu0_true - package_calibration$mu0_superpop) < 1e-12)
  }
  data <- .generate_common_basis_data(config, omega, reference, kappa)
  data_split <- split_data_by_site(data)
  fit <- run_tate_crossfit(
    data_split, n_folds = 5L, communication_mode = "one_round",
    lambda_selection = RoCE:::AGG_WALD_LAMBDA, verbose = FALSE,
    M_tau = 5, M_tau_inference = 5, n_cores = length(N_SOURCE),
    nlambda_init = 100L, family = "binomial", nuisance_lambda_rule = "min"
  )
  summary <- .fit_summary(fit, config, omega, data, reference)
  print(t(summary))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(summary, file.path(stage, "summary.csv"), row.names = FALSE)
    write.csv(data.frame(quantile = names(reference$propensity_quantiles),
                         value = as.numeric(reference$propensity_quantiles)),
              file.path(stage, "propensity_quantiles.csv"), row.names = FALSE)
    saveRDS(list(fit = fit, standardization = reference$standardization,
                 calibration = reference$calibration),
            file.path(stage, "fit.rds"))
    writeLines(c("task=dgp_common_basis", "history_entry=0002",
                 paste0("config=", config), paste0("omega=", omega),
                 paste0("sim_id=", SIM_ID), "rho=0",
                 paste0("code_sha256=", roce_sha256_file("diagnosis/dgp_common_basis/dgp_common_basis.R"))),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "common basis DGP pilot")
}

if (sys.nframe() == 0L) main()
