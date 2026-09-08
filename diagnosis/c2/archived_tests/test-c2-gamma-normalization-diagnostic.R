# test-c2-gamma-normalization-diagnostic.R
#
# Read-only diagnostic for the C2 coverage problem.
#
# The DGP stores multinomial-logit gamma coefficients for
# P(R=s_j,A=a | X) / P(R=t | X).  The source-assisted estimating equations,
# however, average over the conditional source-site sample R=s_j and use
# I(A=a) weights.  This file checks whether the stored raw gamma coefficients
# need an intercept normalization before they are used as "true" density
# ratios in misspecified-outcome DR diagnostics.

library(testthat)

.c2_gamma_norm_enabled <- function() {
  env_enabled <- Sys.getenv("ROCE_RUN_C2_GAMMA_NORMALIZATION", "0") %in%
    c("1", "TRUE", "true", "True")
  filter <- Sys.getenv("TEST_FILTER", "")
  env_enabled || grepl("c2-gamma-normalization", filter, fixed = TRUE)
}

.c2_gamma_norm_int_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0L) default else parsed
}

.c2_gamma_norm_int_vector_env <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(trimws(strsplit(value, ",", fixed = TRUE)[[1L]])))
  parsed <- parsed[is.finite(parsed) & parsed > 0L]
  if (length(parsed) == 0L) default else parsed
}

.c2_gamma_norm_source_env <- function(default) {
  value <- Sys.getenv("ROCE_C2_GAMMA_NORM_SOURCES", unset = "")
  if (!nzchar(value)) return(default)
  parsed <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  parsed <- parsed[nzchar(parsed)]
  if (length(parsed) == 0L) default else parsed
}

.c2_gamma_norm_logistic <- function(eta) {
  eta <- pmax(pmin(eta, LOGISTIC_CLIP), -LOGISTIC_CLIP)
  1 / (1 + exp(-eta))
}

.c2_gamma_norm_predict <- function(W, alpha) {
  .c2_gamma_norm_logistic(as.numeric(cbind(1, as.matrix(W)) %*% as.numeric(alpha)))
}

.c2_gamma_norm_weights <- function(Z, gamma) {
  exp(-as.numeric(cbind(1, as.matrix(Z)) %*% as.numeric(gamma)))
}

.c2_gamma_norm_shift <- function(gamma, log_target_over_source) {
  gamma <- as.numeric(gamma)
  gamma[1L] <- gamma[1L] + log_target_over_source
  gamma
}

.c2_gamma_norm_moments <- function(Z_target, Z_source, A_source, gamma,
                                   A_val = 1L) {
  phi_target <- cbind(1, as.matrix(Z_target))
  phi_source <- cbind(1, as.matrix(Z_source))
  w <- .c2_gamma_norm_weights(Z_source, gamma)
  source_moment <- colMeans(phi_source * (as.numeric(A_source == A_val) * w))
  target_moment <- colMeans(phi_target)
  diff <- source_moment - target_moment
  list(
    moment_l2 = sqrt(mean(diff^2)),
    intercept_error = diff[1L],
    source_intercept = source_moment[1L],
    target_intercept = target_moment[1L],
    balance_gap = {
      treated_w <- w[A_source == A_val]
      if (length(treated_w) == 0L || sum(treated_w) <= 0) {
        NA_real_
      } else {
        zbar <- colSums(as.matrix(Z_source)[A_source == A_val, , drop = FALSE] *
                          treated_w) / sum(treated_w)
        sqrt(mean((colMeans(as.matrix(Z_target)) - zbar)^2))
      }
    }
  )
}

.c2_gamma_norm_dr <- function(W_target, W_source, Z_source, Y_source, A_source,
                              alpha, gamma, mu_superpop, mu_realized,
                              A_val = 1L) {
  m_target <- .c2_gamma_norm_predict(W_target, alpha)
  m_source <- .c2_gamma_norm_predict(W_source, alpha)
  w_source <- .c2_gamma_norm_weights(Z_source, gamma)
  correction <- mean(as.numeric(A_source == A_val) * w_source *
                       (Y_source - m_source))
  estimate <- mean(m_target) + correction
  list(
    estimate = estimate,
    target_mean = mean(m_target),
    correction = correction,
    bias_superpop = estimate - mu_superpop,
    bias_realized = estimate - mu_realized,
    mean_weight_treated = mean(w_source[A_source == A_val]),
    max_weight = max(w_source)
  )
}

.c2_gamma_norm_pop_ratio <- function(data, source_name, A_val = 1L) {
  site_num <- as.integer(gsub("s", "", source_name))
  probs <- calculate_site_probabilities(
    as.matrix(data$Z_site_true %||% data$Z_site),
    data$gamma_params,
    data$K
  )
  p_target <- mean(probs$p_target)
  p_source <- mean(probs[[paste0("p_s", site_num, "_0")]] +
                     probs[[paste0("p_s", site_num, "_1")]])
  list(
    p_target = p_target,
    p_source = p_source,
    log_target_over_source = log(p_target / p_source)
  )
}

.c2_gamma_norm_rows <- function(data, source_names, A_val = 1L,
                                nlambda = 6L) {
  split <- split_data_by_site(data)
  target <- split[["t"]]
  W_target <- as.matrix(target$W_outcome)
  Z_target <- as.matrix(target$Z_site)
  W_target_true <- as.matrix(target$W_outcome_true %||% target$W_outcome)
  Z_target_true <- as.matrix(target$Z_site_true %||% target$Z_site)

  rows <- list()
  for (source_name in source_names) {
    source <- split[[source_name]]
    W_source <- as.matrix(source$W_outcome)
    Z_source <- as.matrix(source$Z_site)
    W_source_true <- as.matrix(source$W_outcome_true %||% source$W_outcome)
    Z_source_true <- as.matrix(source$Z_site_true %||% source$Z_site)
    A_source <- source$A
    Y_source <- source$Y

    site_num <- as.integer(gsub("s", "", source_name))
    gamma_key <- paste0("s", site_num, "_", A_val)
    gamma_raw <- as.numeric(data$gamma_params[[gamma_key]])
    if (length(gamma_raw) == 0L) {
      stop(sprintf("Missing gamma key '%s'.", gamma_key), call. = FALSE)
    }
    gamma_opp_key <- paste0("s", site_num, "_", 1L - A_val)
    gamma_raw_opposite <- as.numeric(data$gamma_params[[gamma_opp_key]])
    if (length(gamma_raw_opposite) == 0L) {
      stop(sprintf("Missing gamma key '%s'.", gamma_opp_key), call. = FALSE)
    }

    pop_ratio <- .c2_gamma_norm_pop_ratio(data, source_name, A_val)
    gamma_shift_pop <- .c2_gamma_norm_shift(gamma_raw, pop_ratio$log_target_over_source)
    gamma_shift_sample <- .c2_gamma_norm_shift(
      gamma_raw,
      log(nrow(W_target) / nrow(W_source))
    )
    gamma_shift_pop_opposite <- .c2_gamma_norm_shift(
      gamma_raw_opposite, pop_ratio$log_target_over_source
    )
    gamma_shift_sample_opposite <- .c2_gamma_norm_shift(
      gamma_raw_opposite,
      log(nrow(W_target) / nrow(W_source))
    )

    mean_phi_target <- c(1, colMeans(Z_target))
    alpha_init <- as.numeric(fit_initial_outcome(
      W_source, Y_source, A_source,
      A_val = A_val, lambda = 0.01, nlambda = nlambda,
      family = "binomial"
    ))
    gamma_fit_l0 <- as.numeric(fit_initial_density_ratio(
      Z_source, A_source, mean_phi_target,
      lambda = 0, A_val = A_val
    ))
    gamma_fit_l001 <- as.numeric(fit_initial_density_ratio(
      Z_source, A_source, mean_phi_target,
      lambda = 0.01, A_val = A_val
    ))
    alpha_cal <- as.numeric(fit_unified_outcome(
      W_outcome = W_source, Y = Y_source, A = A_source,
      A_val = A_val, gamma_s = gamma_fit_l001,
      lambda = 0.01,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT, Z_site = Z_source
    ))
    mean_grad <- .mean_glm_gradient_site_basis(
      W_outcome = W_target, Z_site = Z_target, alpha = alpha_init,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    )
    gamma_cal_l001 <- as.numeric(fit_unified_density_ratio(
      Z_site = Z_source, A = A_source,
      mean_grad_psi = mean_grad, alpha_init = alpha_init,
      lambda = 0.01,
      calibrated = TRUE, M_tau = M_TAU_DEFAULT,
      W_outcome = W_source, A_val = A_val,
      family_int = FAMILY_BINOMIAL, link_int = LINK_LOGIT
    ))

    gamma_cases <- list(
      raw = list(gamma = gamma_raw, Z_t = Z_target_true, Z_s = Z_source_true),
      shift_pop = list(gamma = gamma_shift_pop, Z_t = Z_target_true, Z_s = Z_source_true),
      shift_sample = list(gamma = gamma_shift_sample, Z_t = Z_target_true, Z_s = Z_source_true),
      raw_opposite = list(gamma = gamma_raw_opposite, Z_t = Z_target_true, Z_s = Z_source_true),
      shift_pop_opposite = list(gamma = gamma_shift_pop_opposite, Z_t = Z_target_true, Z_s = Z_source_true),
      shift_sample_opposite = list(gamma = gamma_shift_sample_opposite, Z_t = Z_target_true, Z_s = Z_source_true),
      fit_l0 = list(gamma = gamma_fit_l0, Z_t = Z_target, Z_s = Z_source),
      fit_l001 = list(gamma = gamma_fit_l001, Z_t = Z_target, Z_s = Z_source),
      cal_l001 = list(gamma = gamma_cal_l001, Z_t = Z_target, Z_s = Z_source)
    )
    alpha_cases <- list(
      true = list(alpha = as.numeric(data$alpha1_true),
                  W_t = W_target_true, W_s = W_source_true),
      init = list(alpha = alpha_init, W_t = W_target, W_s = W_source),
      cal = list(alpha = alpha_cal, W_t = W_target, W_s = W_source)
    )

    for (gamma_method in names(gamma_cases)) {
      gc <- gamma_cases[[gamma_method]]
      mom <- .c2_gamma_norm_moments(
        gc$Z_t, gc$Z_s, A_source, gc$gamma, A_val
      )
      for (alpha_method in names(alpha_cases)) {
        ac <- alpha_cases[[alpha_method]]
        dr <- .c2_gamma_norm_dr(
          ac$W_t, ac$W_s, gc$Z_s, Y_source, A_source,
          ac$alpha, gc$gamma,
          mu_superpop = as.numeric(data$mu1_true),
          mu_realized = as.numeric(data$mu1_realized),
          A_val = A_val
        )
        rows[[length(rows) + 1L]] <- data.frame(
          source = source_name,
          alpha_method = alpha_method,
          gamma_method = gamma_method,
          gamma_intercept = gc$gamma[1L],
          raw_intercept = gamma_raw[1L],
          raw_opposite_intercept = gamma_raw_opposite[1L],
          sample_log_target_over_source = log(nrow(W_target) / nrow(W_source)),
          pop_log_target_over_source = pop_ratio$log_target_over_source,
          pop_p_target = pop_ratio$p_target,
          pop_p_source = pop_ratio$p_source,
          moment_l2 = mom$moment_l2,
          intercept_error = mom$intercept_error,
          source_intercept = mom$source_intercept,
          balance_gap = mom$balance_gap,
          estimate = dr$estimate,
          target_mean = dr$target_mean,
          correction = dr$correction,
          bias_superpop = dr$bias_superpop,
          bias_realized = dr$bias_realized,
          mean_weight_treated = dr$mean_weight_treated,
          max_weight = dr$max_weight,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

test_that("c2-gamma-normalization diagnostic compares true-gamma intercept conventions", {
  skip_if_not(.c2_gamma_norm_enabled(),
              message = "run with --filter c2-gamma-normalization or set ROCE_RUN_C2_GAMMA_NORMALIZATION=1")

  n_total <- .c2_gamma_norm_int_env("ROCE_C2_GAMMA_NORM_N", 5000L)
  K_sites <- .c2_gamma_norm_int_env("ROCE_C2_GAMMA_NORM_K", 3L)
  p_values <- .c2_gamma_norm_int_vector_env("ROCE_C2_GAMMA_NORM_P", c(10L, 50L))
  reps <- .c2_gamma_norm_int_env("ROCE_C2_GAMMA_NORM_REPS", 3L)
  nlambda <- .c2_gamma_norm_int_env("ROCE_C2_GAMMA_NORM_NLAMBDA", 6L)
  base_seed <- .c2_gamma_norm_int_env("ROCE_C2_GAMMA_NORM_SEED", 52050L)
  source_names <- .c2_gamma_norm_source_env(c("s1", "s2"))

  rows <- list()
  for (p in p_values) {
    for (rep_idx in seq_len(reps)) {
      seed <- base_seed + 1000L * which(p_values == p)[1L] + rep_idx - 1L
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
        dgp_type = "roce",
        warn_ignored = FALSE
      )
      available <- setdiff(names(split_data_by_site(data)), "t")
      chosen <- intersect(source_names, available)
      if (length(chosen) != length(source_names)) {
        stop(sprintf("Requested sources not available. requested={%s}, available={%s}",
                     paste(source_names, collapse = ","),
                     paste(available, collapse = ",")), call. = FALSE)
      }
      out <- .c2_gamma_norm_rows(data, chosen, A_val = 1L, nlambda = nlambda)
      out$p <- p
      out$rep <- rep_idx
      out$seed <- seed
      out$n_total <- n_total
      out$K <- K_sites
      rows[[length(rows) + 1L]] <- out
    }
  }
  df <- do.call(rbind, rows)

  out_dir <- "c2_gamma_normalization_output"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(
    out_dir,
    paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_c2_gamma_normalization.csv")
  )
  write.csv(df, out_path, row.names = FALSE)

  cat(sprintf("\n[c2-gamma-normalization] CSV: %s\n", out_path))
  cat(sprintf("[c2-gamma-normalization] n=%d K=%d p={%s} reps=%d sources={%s}\n",
              n_total, K_sites, paste(p_values, collapse = ","),
              reps, paste(source_names, collapse = ",")))

  key <- with(df, paste(p, alpha_method, gamma_method, sep = "|"))
  summary <- do.call(rbind, lapply(split(df, key), function(sub) {
    data.frame(
      p = sub$p[1L],
      alpha_method = sub$alpha_method[1L],
      gamma_method = sub$gamma_method[1L],
      bias_superpop_mean = mean(sub$bias_superpop),
      bias_realized_mean = mean(sub$bias_realized),
      moment_l2_mean = mean(sub$moment_l2),
      intercept_error_mean = mean(sub$intercept_error),
      source_intercept_mean = mean(sub$source_intercept),
      balance_gap_mean = mean(sub$balance_gap),
      max_weight_p95 = quantile(sub$max_weight, 0.95, names = FALSE),
      n_rows = nrow(sub),
      stringsAsFactors = FALSE
    )
  }))
  summary <- summary[order(summary$p, summary$alpha_method, summary$gamma_method), ]
  rownames(summary) <- NULL
  print(summary, row.names = FALSE, digits = 4)

  raw_vs_shift <- summary[
    summary$alpha_method %in% c("init", "cal") &
      summary$gamma_method %in% c("raw", "shift_pop", "shift_sample",
                                  "raw_opposite", "shift_pop_opposite",
                                  "shift_sample_opposite",
                                  "fit_l0", "fit_l001", "cal_l001"),
    c("p", "alpha_method", "gamma_method", "bias_superpop_mean",
      "bias_realized_mean", "moment_l2_mean", "intercept_error_mean",
      "source_intercept_mean"),
    drop = FALSE
  ]
  cat("\n[c2-gamma-normalization] Headline rows for misspecified fitted outcome basis:\n")
  print(raw_vs_shift, row.names = FALSE, digits = 4)

  expect_true(nrow(df) > 0L)
  expect_true(all(is.finite(df$bias_superpop)))
  expect_true(all(is.finite(df$moment_l2)))
  expect_true(all(c("raw", "shift_pop", "shift_sample",
                    "raw_opposite", "shift_pop_opposite", "shift_sample_opposite",
                    "fit_l0", "fit_l001", "cal_l001") %in%
                    unique(df$gamma_method)))
})
