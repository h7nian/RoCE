.c2_nuisance_lambda_enabled <- function() {
  env_enabled <- Sys.getenv("FACEHD_RUN_C2_NUISANCE_LAMBDA", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-nuisance-lambda", filter, fixed = TRUE)
}

.c2_nuisance_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_logistic <- function(eta) {
  eta <- pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP)
  1 / (1 + exp(-eta))
}

.c2_predict <- function(W, alpha) {
  .c2_logistic(as.numeric(cbind(1, as.matrix(W)) %*% as.numeric(alpha)))
}

.c2_dr_weights <- function(Z, gamma) {
  exp(-as.numeric(cbind(1, as.matrix(Z)) %*% as.numeric(gamma)))
}

source({
  helper_path <- file.path("diagnosis", "c2", "c2_true_gamma_utils.R")
  if (file.exists(helper_path)) helper_path else file.path("..", "..", helper_path)
})

.c2_ess <- function(w) {
  if (length(w) == 0L) return(NA_real_)
  s1 <- sum(w)
  s2 <- sum(w^2)
  if (s2 <= 0) NA_real_ else s1^2 / s2
}

.c2_balance_gap <- function(Z_target, Z_source, A_source, gamma, A_val = 1L) {
  treated <- A_source == A_val
  w <- .c2_dr_weights(Z_source, gamma)[treated]
  Zs <- as.matrix(Z_source)[treated, , drop = FALSE]
  if (length(w) == 0L || sum(w) <= 0) return(NA_real_)
  target_mean <- colMeans(as.matrix(Z_target))
  source_mean <- colSums(Zs * w) / sum(w)
  sqrt(mean((target_mean - source_mean)^2))
}

.c2_dr_bias <- function(W_target, W_source, Z_source, Y_source, A_source,
                        alpha, gamma, truth, A_val = 1L) {
  m_target <- .c2_predict(W_target, alpha)
  m_source <- .c2_predict(W_source, alpha)
  w_source <- .c2_dr_weights(Z_source, gamma)
  correction <- mean(as.numeric(A_source == A_val) * w_source *
                       (Y_source - m_source))
  mean(m_target) + correction - truth
}

.c2_cv_density_lambda <- function(stage, Z_source, A_source, Z_target,
                                  W_source, W_target, alpha_init,
                                  nlambda, n_cv_folds, A_val = 1L) {
  if (identical(stage, "initial")) {
    mean_phi <- c(1, colMeans(as.matrix(Z_target)))
    lambda_max <- compute_lambda_max_initial_dr(
      Z_source, A_source, mean_phi, A_val = A_val
    )
    lambda_grid <- build_lambda_grid(
      lambda_max = lambda_max,
      lambda_min_ratio = if (nrow(Z_source) > ncol(Z_source)) {
        LAMBDA_MIN_RATIO_LOW_DIM
      } else {
        LAMBDA_MIN_RATIO_HIGH_DIM
      },
      nlambda = nlambda
    )
    cv <- select_lambda_cv_initial_density_ratio_cpp(
      Z_source, A_source, mean_phi, lambda_grid, n_cv_folds,
      MAX_ITER_DEFAULT, TOL_DEFAULT, A_val
    )
    list(cv = cv, lambda_grid = lambda_grid)
  } else {
    mean_grad <- .mean_glm_gradient_site_basis(
      W_target, Z_target, alpha_init, FAMILY_BINOMIAL, LINK_LOGIT
    )
    lambda_max <- compute_lambda_max_refined_dr(
      Z_source, A_source, mean_grad, alpha_init,
      A_val = A_val,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT,
      W_outcome = W_source,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT
    )
    lambda_grid <- build_lambda_grid(
      lambda_max = lambda_max,
      lambda_min_ratio = if (nrow(Z_source) > ncol(Z_source)) {
        LAMBDA_MIN_RATIO_LOW_DIM
      } else {
        LAMBDA_MIN_RATIO_HIGH_DIM
      },
      nlambda = nlambda
    )
    cv <- select_lambda_cv_calibrated_density_ratio_cpp(
      Z_source, A_source, mean_grad, alpha_init, lambda_grid,
      n_cv_folds, MAX_ITER_DEFAULT, TOL_DEFAULT, M_TAU_DEFAULT,
      W_source, A_val, FAMILY_BINOMIAL, LINK_LOGIT
    )
    list(cv = cv, lambda_grid = lambda_grid, mean_grad = mean_grad)
  }
}

.c2_fit_gamma_at_lambda <- function(stage, lambda, Z_source, A_source, Z_target,
                                    W_source, alpha_init, mean_grad = NULL,
                                    A_val = 1L) {
  if (identical(stage, "initial")) {
    fit_initial_density_ratio(
      Z_source, A_source, c(1, colMeans(Z_target)),
      lambda = lambda, A_val = A_val
    )
  } else {
    fit_unified_density_ratio(
      Z_site = Z_source, A = A_source,
      mean_grad_psi = mean_grad,
      alpha_init = alpha_init,
      lambda = lambda,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT,
      W_outcome = W_source, A_val = A_val,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    )
  }
}

test_that("c2-nuisance-lambda probe compares gamma lambda choices directly", {
  skip_if_not(.c2_nuisance_lambda_enabled(),
              message = "set FACEHD_RUN_C2_NUISANCE_LAMBDA=1 or run with --filter c2-nuisance-lambda")

  n_total <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_N", 900L)
  K_sites <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_K", 2L)
  p <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_P", 50L)
  nlambda <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_NLAMBDA", 8L)
  n_cv_folds <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_CV_FOLDS", 3L)
  base_seed <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_SEED", 42050L)
  reps <- .c2_nuisance_int_env("FACEHD_C2_NUISANCE_REPS", 1L)
  A_val <- 1L

  stages <- c("initial", "calibrated")
  rows <- list()

  for (rep_idx in seq_len(reps)) {
    seed <- base_seed + rep_idx - 1L
    set.seed(seed)
    data <- generate_simulation_data(
      n_total = n_total, K = K_sites, p = p,
      config = "C2",
      estimand_type = "superpopulation",
      site_allocation = "model",
      transform_type = "mild",
      outcome_type = "binary",
      heterogeneity_type = "none",
      shift_strength = 0.5,
      dgp_type = "facehd",
      warn_ignored = FALSE
    )
    split <- split_data_by_site(data)
    source_name <- "s1"
    target <- split[["t"]]
    source <- split[[source_name]]

    W_target <- as.matrix(target$W_outcome)
    Z_target <- as.matrix(target$Z_site)
    W_source <- as.matrix(source$W_outcome)
    Z_source <- as.matrix(source$Z_site)
    A_source <- source$A
    Y_source <- source$Y
    truth <- as.numeric(data$mu1_true)

    alpha_init <- fit_initial_outcome(
      W_source, Y_source, A_source, A_val = A_val,
      lambda = 0.01, nlambda = nlambda, family = "binomial"
    )

    gamma_true <- c2_true_source_calibration_gamma(data, source_name, A_val)
    gamma_true_bias <- .c2_dr_bias(
      W_target, W_source, Z_source, Y_source, A_source,
      alpha_init, gamma_true, truth, A_val
    )

    for (stage in stages) {
      cv_info <- .c2_cv_density_lambda(
        stage = stage,
        Z_source = Z_source, A_source = A_source, Z_target = Z_target,
        W_source = W_source, W_target = W_target,
        alpha_init = alpha_init,
        nlambda = nlambda, n_cv_folds = n_cv_folds, A_val = A_val
      )
      cv <- cv_info$cv
      lambda_grid <- as.numeric(cv_info$lambda_grid)
      for (li in seq_along(lambda_grid)) {
        lambda <- lambda_grid[[li]]
        selected_min <- isTRUE(all.equal(lambda, cv$lambda_min, tolerance = 1e-12))
        selected_1se <- isTRUE(all.equal(lambda, cv$lambda_1se, tolerance = 1e-12))
        gamma_hat <- .c2_fit_gamma_at_lambda(
          stage = stage, lambda = lambda,
          Z_source = Z_source, A_source = A_source, Z_target = Z_target,
          W_source = W_source, alpha_init = alpha_init,
          mean_grad = cv_info$mean_grad, A_val = A_val
        )
        treated_w <- .c2_dr_weights(Z_source, gamma_hat)[A_source == A_val]
        rows[[length(rows) + 1L]] <- data.frame(
          rep = rep_idx,
          seed = seed,
          source = source_name,
          stage = stage,
          lambda_index = li,
          lambda = lambda,
          selected_min = selected_min,
          selected_1se = selected_1se,
          lambda_min = cv$lambda_min,
          lambda_1se = cv$lambda_1se,
          lambda_1se_over_min = cv$lambda_1se / cv$lambda_min,
          support = sum(abs(as.numeric(gamma_hat)[-1L]) > 1e-8),
          gamma_l2 = sqrt(mean(as.numeric(gamma_hat)^2)),
          gamma_diff_true_l2 = sqrt(mean((as.numeric(gamma_hat) - gamma_true)^2)),
          balance_gap = .c2_balance_gap(Z_target, Z_source, A_source, gamma_hat, A_val),
          ess_treated = .c2_ess(treated_w),
          max_weight_treated = max(treated_w),
          dr_bias = .c2_dr_bias(
            W_target, W_source, Z_source, Y_source, A_source,
            alpha_init, gamma_hat, truth, A_val
          ),
          gamma_true_dr_bias = gamma_true_bias,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  df <- do.call(rbind, rows)
  df$n_total <- n_total
  df$K <- K_sites
  df$p <- p

  out_dir <- "c2_nuisance_lambda_output"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(
    out_dir,
    paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_c2_nuisance_lambda.csv")
  )
  write.csv(df, out_path, row.names = FALSE)

  cat(sprintf("\n[c2-nuisance-lambda] CSV: %s\n", out_path))
  cat(sprintf("[c2-nuisance-lambda] reps=%d, seed range=%d:%d\n",
              reps, base_seed, base_seed + reps - 1L))
  print(df[, c("rep", "seed", "stage", "lambda_index", "lambda",
               "selected_min", "selected_1se", "support",
               "gamma_diff_true_l2", "balance_gap", "ess_treated",
               "max_weight_treated", "dr_bias", "gamma_true_dr_bias")],
        row.names = FALSE, digits = 4)

  selected <- df[df$selected_min | df$selected_1se, , drop = FALSE]
  best_abs_bias <- df[order(abs(df$dr_bias), df$stage, df$lambda_index),
                      c("rep", "seed", "stage", "lambda_index", "lambda",
                        "selected_min", "selected_1se", "support",
                        "balance_gap", "ess_treated", "dr_bias",
                        "gamma_true_dr_bias"),
                      drop = FALSE]
  cat("\n[c2-nuisance-lambda] CV-selected rows:\n")
  print(selected[, c("rep", "seed", "stage", "lambda_index", "lambda",
                     "selected_min", "selected_1se", "support", "balance_gap",
                     "ess_treated", "dr_bias", "gamma_true_dr_bias")],
        row.names = FALSE, digits = 4)
  cat("\n[c2-nuisance-lambda] Rows sorted by |DR bias|:\n")
  print(best_abs_bias, row.names = FALSE, digits = 4)

  expect_equal(nrow(df), reps * length(stages) * nlambda)
  expect_true(any(df$selected_min))
  expect_true(any(df$selected_1se))
  expect_true(all(is.finite(df$lambda)))
  expect_true(all(is.finite(df$dr_bias)))
  expect_true(all(is.finite(df$balance_gap)))
  expect_true(all(df$lambda_1se >= df$lambda_min))
})
