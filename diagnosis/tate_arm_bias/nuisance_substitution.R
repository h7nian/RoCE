#!/usr/bin/env Rscript
# Which nuisance carries the residual confounding?
#
# Target site only, cross-fitted AIPW for each arm, in a 2x2 design: the outcome
# model and the propensity are each either the TRUE function from the generator
# or a cross-fitted lasso on the same working basis the estimator uses. The
# (true, true) cell is the oracle and must be unbiased; whichever substitution
# brings the bias back is the nuisance responsible.
#
# Usage: nuisance_substitution.R CONFIG FIRST_SEED LAST_SEED [OUT_CSV]
suppressWarnings(suppressPackageStartupMessages(library(glmnet)))
args <- commandArgs(trailingOnly = TRUE)
configuration <- if (length(args) >= 1L) args[[1L]] else "C1"
first_seed <- if (length(args) >= 2L) as.integer(args[[2L]]) else 1L
last_seed <- if (length(args) >= 3L) as.integer(args[[3L]]) else 50L
out_csv <- if (length(args) >= 4L) args[[4L]] else ""
seeds <- seq(first_seed, last_seed)

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) .libPaths(c(project_library, .libPaths()))
suppressPackageStartupMessages(library(RoCE))

kappa <- RoCE:::FACE_KAPPA
strengths <- RoCE:::.face_misspecification_strengths(
  configuration, RoCE:::FACE_MISSPECIFICATION_STRENGTH
)
calibration <- RoCE:::get_face_binary_calibration(
  p = 100L, config = configuration, kappa = kappa,
  misspecification_strength = strengths$outcome
)
n_folds <- 10L

# Cross-fitted lasso prediction on the working basis, fitted on `subset` rows.
crossfit_predict <- function(basis, response, folds, subset_mask) {
  prediction <- rep(NA_real_, length(response))
  for (k in seq_len(n_folds)) {
    train <- which(folds != k & subset_mask)
    test <- which(folds == k)
    if (length(train) < 30L || length(unique(response[train])) < 2L) next
    fit <- glmnet::cv.glmnet(basis[train, , drop = FALSE], response[train],
                             family = "binomial", standardize = TRUE, nfolds = 5L)
    prediction[test] <- as.numeric(stats::predict(
      fit, newx = basis[test, , drop = FALSE], s = "lambda.min", type = "response"
    ))
  }
  prediction
}

one_seed <- function(seed) {
  set.seed(seed)
  d <- RoCE:::generate_simulation_data(
    n_total = 3000L, K = 2L, p = 100L, config = configuration,
    estimand_type = "superpopulation", outcome_type = "binary",
    dgp_type = "face", ate_deviation = 0, n_deviated_sites = 0L,
    deviation_mechanism = "treated_arm", warn_ignored = FALSE
  )
  target <- d$R == "t"
  eta <- RoCE:::.face_mixed_predictor(
    d$X, d$X_dagger, d$out_params$beta_linear, d$out_params$beta_squared,
    kappa, strengths$outcome
  )
  g <- RoCE:::face_binary_logit(eta, calibration)
  m1_true <- RoCE:::logistic(g + 1)[target]
  m0_true <- RoCE:::logistic(g)[target]
  e_true <- d$p_treat_true[target]

  basis <- as.matrix(d$W_outcome)[target, , drop = FALSE]
  y <- d$Y[target]; a <- d$A[target]
  set.seed(seed + 1e6)
  folds <- sample(rep_len(seq_len(n_folds), length(y)))

  m1_hat <- crossfit_predict(basis, y, folds, a == 1L)
  m0_hat <- crossfit_predict(basis, y, folds, a == 0L)
  e_hat <- crossfit_predict(basis, a, folds, rep(TRUE, length(a)))
  e_hat <- pmin(pmax(e_hat, 0.01), 0.99)

  arm_means <- function(m1, m0, e) {
    c(mu1 = mean(m1 + a * (y - m1) / e),
      mu0 = mean(m0 + (1 - a) * (y - m0) / (1 - e)))
  }
  # Double cross-fitting (Zivich & Breskin 2021): three disjoint splits, with the
  # outcome model and the propensity fitted on DIFFERENT splits and both applied
  # to the third. Because the two nuisances never share training data their
  # errors are independent, so the second-order product term that survives when
  # both are estimated on the same split has expectation near zero.
  split3 <- ((folds - 1L) %% 3L) + 1L
  m1_dcf <- rep(NA_real_, length(y)); m0_dcf <- m1_dcf; e_dcf <- m1_dcf
  for (s3 in 1:3) {
    outcome_split <- (s3 %% 3L) + 1L
    propensity_split <- ((s3 + 1L) %% 3L) + 1L
    predict_on <- which(split3 == s3)
    fit_arm <- function(arm) {
      train <- which(split3 == outcome_split & a == arm)
      if (length(train) < 30L || length(unique(y[train])) < 2L) return(rep(NA_real_, length(predict_on)))
      f <- glmnet::cv.glmnet(basis[train, , drop = FALSE], y[train],
                             family = "binomial", standardize = TRUE, nfolds = 5L)
      as.numeric(stats::predict(f, newx = basis[predict_on, , drop = FALSE],
                                s = "lambda.min", type = "response"))
    }
    m1_dcf[predict_on] <- fit_arm(1L)
    m0_dcf[predict_on] <- fit_arm(0L)
    train_ps <- which(split3 == propensity_split)
    fp <- glmnet::cv.glmnet(basis[train_ps, , drop = FALSE], a[train_ps],
                            family = "binomial", standardize = TRUE, nfolds = 5L)
    e_dcf[predict_on] <- as.numeric(stats::predict(
      fp, newx = basis[predict_on, , drop = FALSE], s = "lambda.min", type = "response"))
  }
  e_dcf <- pmin(pmax(e_dcf, 0.01), 0.99)

  cells <- list(
    oracle        = arm_means(m1_true, m0_true, e_true),
    est_outcome   = arm_means(m1_hat,  m0_hat,  e_true),
    est_propensity= arm_means(m1_true, m0_true, e_hat),
    both_est      = arm_means(m1_hat,  m0_hat,  e_hat),
    double_cf     = arm_means(m1_dcf,  m0_dcf,  e_dcf)
  )
  out <- do.call(rbind, lapply(names(cells), function(nm) {
    v <- cells[[nm]]
    data.frame(cell = nm, mu1 = v[["mu1"]], mu0 = v[["mu0"]],
               tate = v[["mu1"]] - v[["mu0"]], stringsAsFactors = FALSE)
  }))
  out$mu1_true <- mean(m1_true); out$mu0_true <- mean(m0_true)
  out$seed <- seed
  out
}

res <- do.call(rbind, lapply(seeds, one_seed))
res$config <- configuration
if (nzchar(out_csv)) utils::write.csv(res, out_csv, row.names = FALSE)
n_seeds <- length(seeds)
res$tate_true <- res$mu1_true - res$mu0_true
cat(sprintf("\nNUISANCE SUBSTITUTION  config=%s  seeds=%d  (target site, cross-fitted AIPW)\n",
            configuration, n_seeds))
cat(sprintf("%-15s %12s %12s %12s\n", "nuisances", "mu1 bias", "mu0 bias", "TATE bias"))
for (nm in c("oracle", "est_outcome", "est_propensity", "both_est", "double_cf")) {
  s <- res[res$cell == nm, ]
  cat(sprintf("%-15s %+12.5f %+12.5f %+12.5f   (mcse %.5f)\n", nm,
              mean(s$mu1 - s$mu1_true), mean(s$mu0 - s$mu0_true),
              mean(s$tate - s$tate_true),
              sd(s$tate - s$tate_true) / sqrt(nrow(s))))
}
