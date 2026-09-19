#!/usr/bin/env Rscript
# Oracle check: is the target-side estimand/truth/influence-function chain sound?
#
# Builds the TRUE per-observation outcome probabilities and uses the TRUE
# propensity that the generator actually drew treatment with, then forms the
# ordinary doubly robust estimator on the target site only. With correct
# nuisances this must be unbiased for the target ATE, so any systematic
# departure from mu1_true - mu0_true indicts the estimand definition, the truth
# constant, or the harness rather than nuisance estimation.
#
# Nothing is fitted and nothing is modified; every quantity is reconstructed
# from what generate_face_data() already returns plus the generator's own
# helpers, so this measures the shipped DGP, not a re-implementation of it.
#
# Usage: oracle_check.R CONFIG K N_SEEDS
args <- commandArgs(trailingOnly = TRUE)
configuration <- if (length(args) >= 1L) args[[1L]] else "C1"
source_count <- if (length(args) >= 2L) as.integer(args[[2L]]) else 2L
n_seeds <- if (length(args) >= 3L) as.integer(args[[3L]]) else 300L

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))

n_site <- 1000L
dimension <- 100L
kappa <- RoCE:::FACE_KAPPA

strengths <- RoCE:::.face_misspecification_strengths(
  configuration, RoCE:::FACE_MISSPECIFICATION_STRENGTH
)
outcome_strength <- strengths$outcome
calibration <- RoCE:::get_face_binary_calibration(
  p = dimension, config = configuration, kappa = kappa,
  misspecification_strength = outcome_strength
)

one_seed <- function(seed) {
  set.seed(seed)
  d <- RoCE:::generate_simulation_data(
    n_total = n_site * (source_count + 1L), K = source_count, p = dimension,
    config = configuration, estimand_type = "superpopulation",
    outcome_type = "binary", dgp_type = "face",
    ate_deviation = 0, n_deviated_sites = 0L,
    deviation_mechanism = "treated_arm", warn_ignored = FALSE
  )
  target <- d$R == "t"
  # True conditional outcome probabilities, exactly as generate_face_outcomes()
  # forms them: logistic(face_binary_logit(eta) + shift + tau * a).
  eta <- RoCE:::.face_mixed_predictor(
    d$X, d$X_dagger, d$out_params$beta_linear, d$out_params$beta_squared,
    kappa, outcome_strength
  )
  g <- RoCE:::face_binary_logit(eta, calibration)
  delta <- unname(d$ate_map[d$R])
  shift <- if (is.null(d$outcome_shift_map)) rep(0, length(d$R)) else unname(d$outcome_shift_map[d$R])
  mu1 <- RoCE:::logistic(g + shift + delta)
  mu0 <- RoCE:::logistic(g + shift)

  y <- d$Y[target]; a <- d$A[target]; e <- d$p_treat_true[target]
  m1 <- mu1[target]; m0 <- mu0[target]
  psi <- (m1 - m0) + a * (y - m1) / e - (1 - a) * (y - m0) / (1 - e)
  c(oracle = mean(psi),
    plugin = mean(m1 - m0),
    truth = as.numeric(d$mu1_true - d$mu0_true),
    min_e = min(e), max_e = max(e))
}

res <- as.data.frame(do.call(rbind, lapply(seq_len(n_seeds), one_seed)))
cat(sprintf("\nORACLE CHECK  config=%s  K=%d  seeds=%d  (no fitting)\n",
            configuration, source_count, n_seeds))
cat(sprintf("  superpopulation truth            = %.6f\n", res$truth[1]))
cat(sprintf("  oracle DR estimate   mean        = %.6f\n", mean(res$oracle)))
cat(sprintf("  oracle DR bias                   = %+.6f   (MC se %.6f)\n",
            mean(res$oracle) - res$truth[1], sd(res$oracle) / sqrt(n_seeds)))
cat(sprintf("  oracle DR empirical sd           = %.6f\n", sd(res$oracle)))
cat(sprintf("  plug-in only (no DR correction)  = %+.6f bias\n",
            mean(res$plugin) - res$truth[1]))
cat(sprintf("  true propensity range            = [%.4f, %.4f]\n",
            min(res$min_e), max(res$max_e)))
cat(sprintf("\n  for comparison, the fitted estimator's bias here is about %s\n",
            if (configuration == "C3") "-0.0140" else "-0.0104"))
