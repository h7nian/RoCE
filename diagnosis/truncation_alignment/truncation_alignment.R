#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0003 Truncation alignment of the tilting calibrated loss
# Task:    (A) On a synthetic two-site example with an active truncation radius,
#          check that the calibrated tilting fit solves the TRUNCATED score
#          E_t[g' Z] = E_s[I exp(-T_M(Z'gamma)) g' Z] (candidate library) versus the
#          untruncated score (baseline library). (B, candidate only) Refit the #0002
#          C1 cell and compare with its stored fit: no unit exceeds the radius, so
#          the common-basis C1 results must be unchanged.
#          Usage: truncation_alignment.R OUTPUT_DIRECTORY LABEL (candidate|baseline)

suppressPackageStartupMessages(library(RoCE))
`%||%` <- RoCE:::`%||%`
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

TRUNCATION_RADIUS <- 2
C1_REFERENCE_FIT <- "diagnosis/out/dgp_common_basis/C1_omega0.00/fit.rds"

.logit_derivative <- function(eta) {
  mu <- RoCE:::logistic(eta)
  mu * (1 - mu)
}

# Two sites, p = 3, a strong mean shift on the first coordinate so that the
# fitted tilt exceeds the radius for a non-trivial share of source units.
.synthetic_sites <- function(n = 800L, p = 3L) {
  RoCE:::with_seed(20260907L, {
    Z_target <- matrix(stats::rnorm(n * p), n, p)
    Z_target[, 1] <- Z_target[, 1] + 1.2
    Z_source <- matrix(stats::rnorm(n * p), n, p)
    A_source <- stats::rbinom(n, 1L, 0.5)
    list(Z_target = Z_target, Z_source = Z_source, A_source = A_source)
  })
}

.score_check <- function() {
  sites <- .synthetic_sites()
  alpha_init <- c(0, 0.3, 0, 0)
  family_int <- RoCE:::FAMILY_BINOMIAL
  link_int <- RoCE:::LINK_LOGIT
  mean_grad_psi <- RoCE:::.mean_glm_gradient_site_basis(
    sites$Z_target, sites$Z_target, alpha_init, family_int, link_int
  )
  fit <- RoCE:::fit_unified_density_ratio_cpp(
    sites$Z_source, sites$A_source, mean_grad_psi, alpha_init,
    0, 20000L, 1e-10, TRUE, TRUNCATION_RADIUS, sites$Z_source, 1L,
    family_int, link_int, numeric(0)
  )
  gamma <- as.numeric(fit$gamma)
  treated <- sites$A_source == 1L
  Z_int <- cbind(1, sites$Z_source[treated, , drop = FALSE])
  predictor <- drop(Z_int %*% gamma)
  truncated <- pmax(pmin(predictor, TRUNCATION_RADIUS), -TRUNCATION_RADIUS)
  outcome_predictor <- drop(Z_int %*% alpha_init)
  psi_prime <- .logit_derivative(pmax(pmin(outcome_predictor, TRUNCATION_RADIUS), -TRUNCATION_RADIUS))
  arm_fraction <- mean(treated)
  source_moment <- function(weight) {
    arm_fraction * colMeans(Z_int * (weight * psi_prime))
  }
  data.frame(
    converged = isTRUE(fit$converged),
    iterations = as.integer(fit$iterations),
    truncation_fraction = mean(abs(predictor) > TRUNCATION_RADIUS),
    max_abs_predictor = max(abs(predictor)),
    truncated_score_max_abs = max(abs(mean_grad_psi - source_moment(exp(-truncated)))),
    untruncated_score_max_abs = max(abs(mean_grad_psi - source_moment(exp(-predictor)))),
    gamma_norm = sqrt(sum(gamma^2))
  )
}

.c1_identity_check <- function() {
  source("diagnosis/dgp_common_basis/dgp_common_basis.R", local = TRUE)
  kappa <- RoCE:::FACE_KAPPA
  reference <- .reference_population(0, 0, kappa)
  data <- .generate_common_basis_data("C1", 0, reference, kappa)
  fit <- run_tate_crossfit(
    split_data_by_site(data), n_folds = 5L, communication_mode = "one_round",
    lambda_selection = RoCE:::AGG_WALD_LAMBDA, verbose = FALSE,
    M_tau = 5, M_tau_inference = 5, n_cores = length(N_SOURCE),
    nlambda_init = 100L, family = "binomial", nuisance_lambda_rule = "min"
  )
  stored <- readRDS(C1_REFERENCE_FIT)$fit
  difference <- function(a, b) max(abs(as.numeric(a) - as.numeric(b)))
  # The stored fit predates the weight-layer SE (HISTORY #0005); compare the
  # fixed-weight SE, which the truncation change must leave untouched.
  data.frame(
    estimate_difference = difference(fit$estimate, stored$estimate),
    se_difference = difference(fit$se_fixed_weights %||% fit$se, stored$se_fixed_weights %||% stored$se),
    fold_weight_difference = difference(fit$fold_weights, stored$fold_weights),
    target_only_difference = difference(
      c(fit$target_only$estimate, fit$target_only$se),
      c(stored$target_only$estimate, stored$target_only$se)
    ),
    source_estimate_difference = difference(fit$source_estimates, stored$source_estimates),
    truncation_fraction = fit$clip_diagnostics$total$logit_truncation_fraction %||% NA_real_
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L || !args[2] %in% c("candidate", "baseline")) {
    stop("usage: truncation_alignment.R OUTPUT_DIRECTORY candidate|baseline")
  }
  output <- file.path(args[1], args[2])
  if (file.exists(output)) stop("output already exists: ", output)
  score <- .score_check()
  print(t(score))
  identity <- if (args[2] == "candidate") .c1_identity_check() else NULL
  if (!is.null(identity)) print(t(identity))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(score, file.path(stage, "score_check.csv"), row.names = FALSE)
    if (!is.null(identity)) {
      write.csv(identity, file.path(stage, "c1_identity.csv"), row.names = FALSE)
    }
    writeLines(c("task=truncation_alignment", "history_entry=0003",
                 paste0("library=", args[2]),
                 paste0("package_library=", .libPaths()[1]),
                 paste0("truncation_radius=", TRUNCATION_RADIUS)),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "truncation alignment check")
}

if (sys.nframe() == 0L) main()
