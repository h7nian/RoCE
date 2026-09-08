#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0001 Weight-layer influence function for the soft-threshold aggregation rule
# Task:    Delta-method derivative of the production (lambda*t-1)_+ soft-threshold
#          weights through the inner-fold moments; finite-difference checks on the
#          seed 10013 bundle, then replay over the 100 saved v19 seeds.
#          Reads results/direct_tate_mc500_b5000/independent_inference_pilot_v19/,
#          writes diagnosis/out/weight_layer/. No nuisance refits, no resampling.

suppressPackageStartupMessages(library(RoCE))
`%||%` <- RoCE:::`%||%`
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")
# Empirical-mass bookkeeping (inner records, pooled centered moments, outer
# values/increments) is shared with the earlier quadratic candidate.
source("diagnosis/tate_common_weight/quadratic_bias_weights.R")
source("diagnosis/tate_common_weight/quadratic_weight_influence.R")

PILOT_ROOT <- "results/direct_tate_mc500_b5000/independent_inference_pilot_v19"
RHO_VALUES <- c(0, 0.5, 1, 1.5, 2, 2.5)
KINK_TOLERANCE <- 1e-6

.extract_tate_fit <- function(bundle, rho) {
  group <- bundle$group_result %||% bundle$grouped
  fit <- group$artifacts[[as.character(rho)]]$direct_tate_results$one_round_crossfit
  if (is.null(fit)) stop("bundle lacks the one_round_crossfit TATE fit for rho ", rho)
  fit
}

.extract_tate_row <- function(bundle, rho) {
  group <- bundle$group_result %||% bundle$grouped
  rows <- group$results[[as.character(rho)]]
  row <- rows[rows$method == "one_round_crossfit_ate", , drop = FALSE]
  if (nrow(row) != 1L) stop("expected one one_round_crossfit_ate row for rho ", rho)
  row
}

# Solve the production weight rule from a moment list, exactly as
# .compute_phase2_weights() does (same floors, no clipping).
.soft_threshold_weights <- function(moments, lambda) {
  eta <- optimize_weights(
    moments$avg_source_est,
    list(V_ot = moments$V_ot, V_t = moments$V_t, V_s = moments$V_s),
    moments$C_ot,
    list(n_t = moments$n_t, n_s = moments$n_s),
    lambda, moments$avg_target_est, moments$C_cross,
    clip_weights = FALSE
  )
  if (attr(eta, "psd_ridge") != 0) stop("PSD ridge active; derivative assumes ridge = 0")
  as.numeric(eta)
}

# Quadratic-form coefficients of N_all * Var(eta) = N_all * (c + 2 l'eta + eta'Q eta)
.variance_quadratic_form <- function(m) {
  K <- length(m$V_t)
  Q <- matrix(m$V_ot, K, K)
  Q <- Q - outer(m$C_ot, rep(1, K)) - outer(rep(1, K), m$C_ot)
  cross <- m$C_cross
  diag(cross) <- 0
  Q <- (Q + cross) / m$n_t
  diag(Q) <- diag(Q) + m$V_t / m$n_t + m$V_s / m$n_s
  list(Q = Q, l = (m$C_ot - m$V_ot) / m$n_t, N_all = m$n_t + sum(m$n_s))
}

# Wald ingredients for the penalty p_j = (lambda t_j - 1)_+ , t_j = |d_j| / s_j.
.wald_ingredients <- function(m, lambda) {
  d <- m$avg_target_est - m$avg_source_est
  s2 <- m$V_ot / m$n_t + m$V_t / m$n_t + m$V_s / m$n_s - 2 * m$C_ot / m$n_t
  if (any(s2 <= RoCE:::VARIANCE_MIN)) stop("discrepancy variance floor active")
  s <- sqrt(s2)
  t <- abs(d) / s
  list(d = d, s = s, t = t, active_penalty = lambda * t > 1,
       kink = abs(lambda * t - 1) < KINK_TOLERANCE)
}

# Per-observation indirect gradient of the estimator through one outer fold's
# weights. `data` is .quadratic_inner_moments() output; `increment` is the
# estimator derivative with respect to this fold's weights.
.soft_threshold_fold_influence <- function(data, eta, lambda, increment) {
  m <- data$moments
  K <- length(eta)
  form <- .variance_quadratic_form(m)
  wald <- .wald_ingredients(m, lambda)
  active <- which(abs(eta) > 1e-12)
  target_grad <- numeric(nrow(data$target_raw))
  source_grad <- lapply(m$n_s, function(n) numeric(n))
  if (length(active) == 0L) {
    return(list(target = target_grad, source = source_grad, kinks = sum(wald$kink)))
  }
  sign_eta <- sign(eta[active])
  residual <- form$Q[active, active, drop = FALSE] %*% eta[active] + form$l[active] +
    sign_eta * pmax(lambda * wald$t[active] - 1, 0) / (2 * form$N_all)
  if (max(abs(residual)) > 1e-8 * max(1, max(abs(form$l)))) {
    stop("stored weights do not satisfy the soft-threshold KKT conditions: ",
         max(abs(residual)))
  }
  h <- drop(solve(form$Q[active, active, drop = FALSE], increment[active]))
  W_t <- m$n_t
  c0 <- data$target_centered[, 1L]
  cs <- data$target_centered[, -1L, drop = FALSE]
  y0 <- data$target_raw[, 1L]
  ys <- data$target_raw[, -1L, drop = FALSE]
  dV_ot <- (c0^2 - m$V_ot) / W_t
  dC_ot <- sweep(c0 * cs, 2L, m$C_ot, "-") / W_t
  dV_t <- sweep(cs^2, 2L, m$V_t, "-") / W_t
  dmu_ot <- (y0 - m$avg_target_est) / W_t
  dmu_pred <- sweep(ys, 2L, m$mu_pred_ts, "-") / W_t
  penalty_scale <- lambda * wald$active_penalty / (2 * form$N_all)
  for (a in seq_along(active)) {
    j <- active[a]
    dQ_eta <- numeric(length(c0))
    for (b in seq_along(active)) {
      k <- active[b]
      dCc <- if (j == k) 0 else (cs[, j] * cs[, k] - m$C_cross[j, k]) / W_t
      dQ_eta <- dQ_eta + eta[k] * (dV_ot - dC_ot[, j] - dC_ot[, k] + dCc) / m$n_t
    }
    dQ_eta <- dQ_eta + eta[j] * dV_t[, j] / m$n_t
    dl <- (dC_ot[, j] - dV_ot) / m$n_t
    ds2 <- (dV_ot + dV_t[, j] - 2 * dC_ot[, j]) / m$n_t
    dd <- dmu_ot - dmu_pred[, j]
    dt <- sign(wald$d[j]) * dd / wald$s[j] - abs(wald$d[j]) * ds2 / (2 * wald$s[j]^3)
    v <- dQ_eta + dl + sign_eta[a] * penalty_scale[j] * dt
    target_grad <- target_grad - h[a] * v
    # Source j observations enter only through V_s,j and delta_j.
    c_s <- data$source_centered[[j]]
    r_s <- data$source_raw[[j]]
    dV_s <- (c_s^2 - m$V_s[j]) / m$n_s[j]
    ddelta <- (r_s - m$delta_ts[j]) / m$n_s[j]
    dt_s <- sign(wald$d[j]) * (-ddelta) / wald$s[j] -
      abs(wald$d[j]) * (dV_s / m$n_s[j]) / (2 * wald$s[j]^3)
    v_s <- eta[j] * dV_s / m$n_s[j] + sign_eta[a] * penalty_scale[j] * dt_s
    source_grad[[j]] <- source_grad[[j]] - h[a] * v_s
  }
  list(target = target_grad, source = source_grad, kinks = sum(wald$kink))
}

# Full per-observation gradient (direct + weight layer) and variances.
.weight_layer_influence <- function(fit) {
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  mass <- lapply(sizes, function(n) rep(1, n))
  inner <- fit$intermediates$inner_fold_info
  lambdas <- fit$fold_lambdas
  records <- lapply(inner, .quadratic_inner_records)
  centered <- lapply(records, .quadratic_inner_moments, mass = mass)
  for (k in seq_along(inner)) {
    for (field in c("V_ot", "V_t", "V_s", "C_ot", "C_cross", "avg_target_est", "avg_source_est")) {
      stopifnot(max(abs(inner[[k]][[field]] - centered[[k]]$moments[[field]])) < 1e-10)
    }
  }
  weights <- t(vapply(seq_along(inner), function(k) {
    .soft_threshold_weights(centered[[k]]$moments, lambdas[k])
  }, numeric(length(sizes) - 1L)))
  if (max(abs(weights - fit$fold_weights)) > 1e-10) {
    stop("re-solved weights differ from stored fold weights by ",
         max(abs(weights - fit$fold_weights)))
  }
  outer <- .quadratic_outer_values(fit, weights)
  direct <- lapply(seq_along(sizes), function(s) (outer$values[[s]] - mean(outer$values[[s]])) / sizes[s])
  indirect <- lapply(sizes, function(n) numeric(n))
  kinks <- 0L
  for (k in seq_along(inner)) {
    fold <- .soft_threshold_fold_influence(centered[[k]], weights[k, ], lambdas[k], outer$increments[k, ])
    kinks <- kinks + fold$kinks
    ids <- unlist(records[[k]]$target_ids, use.names = FALSE)
    stopifnot(!anyDuplicated(ids), !any(ids %in% fit$intermediates$fold_info[[k]]$target_idx))
    indirect[[1L]][ids] <- indirect[[1L]][ids] + fold$target
    for (j in seq_along(fold$source)) {
      ids <- unlist(records[[k]]$source_ids[[j]], use.names = FALSE)
      stopifnot(!anyDuplicated(ids), !any(ids %in% fit$intermediates$fold_info[[k]]$source_idx[[j]]))
      indirect[[j + 1L]][ids] <- indirect[[j + 1L]][ids] + fold$source[[j]]
    }
  }
  full <- Map(`+`, direct, indirect)
  stopifnot(max(abs(vapply(full, sum, numeric(1)))) < 1e-10,
            max(abs(vapply(indirect, sum, numeric(1)))) < 1e-10)
  fixed_variance <- sum(unlist(direct)^2)
  indirect_variance <- sum(unlist(indirect)^2)
  cross_term <- 2 * sum(unlist(direct) * unlist(indirect))
  baseline <- .common_tate_from_weights(fit, weights)
  estimate <- sum(vapply(outer$values, mean, numeric(1)))
  stopifnot(abs(baseline$variance - fixed_variance) < 1e-12,
            abs(baseline$estimate - estimate) < 1e-12,
            abs(fit$estimate - estimate) < 1e-12,
            abs(fit$variance - fixed_variance) < 1e-12)
  list(estimate = estimate, weights = weights, direct = direct, indirect = indirect,
       gradient = full, fixed_variance = fixed_variance,
       weight_layer_variance = sum(unlist(full)^2),
       indirect_variance = indirect_variance, direct_indirect_cross_term = cross_term,
       kinks = kinks, fold_source_cells = length(weights))
}

# The estimator as a functional of observation masses, re-solving the
# production weight rule; used only for finite-difference checks.
.weight_layer_functional <- function(fit, mass) {
  sizes <- c(fit$intermediates$sample_sizes$n_t, fit$intermediates$sample_sizes$n_source)
  inner <- fit$intermediates$inner_fold_info
  weights <- t(vapply(seq_along(inner), function(k) {
    records <- .quadratic_inner_records(inner[[k]])
    .soft_threshold_weights(.quadratic_inner_moments(records, mass)$moments, fit$fold_lambdas[k])
  }, numeric(length(sizes) - 1L)))
  outer <- .quadratic_outer_values(fit, weights)
  sum(vapply(seq_along(sizes), function(s) sum(mass[[s]] * outer$values[[s]]) / sum(mass[[s]]), numeric(1)))
}

.directional_checks <- function(bundle_dir, n_directions = 8L, epsilon = 1e-5) {
  bundle <- readRDS(file.path(bundle_dir, "artifacts.rds"))
  checks <- list()
  for (rho in RHO_VALUES) {
    fit <- .extract_tate_fit(bundle, rho)
    result <- .weight_layer_influence(fit)
    masses <- lapply(result$gradient, function(x) rep(1, length(x)))
    stopifnot(abs(.weight_layer_functional(fit, masses) - result$estimate) < 1e-12)
    for (j in seq_len(n_directions)) {
      direction <- lapply(seq_along(masses), function(s) sin(seq_along(masses[[s]]) * j + s))
      plus <- Map(function(x, d) x + epsilon * d, masses, direction)
      minus <- Map(function(x, d) x - epsilon * d, masses, direction)
      numerical <- (.weight_layer_functional(fit, plus) - .weight_layer_functional(fit, minus)) / (2 * epsilon)
      analytical <- sum(unlist(result$gradient) * unlist(direction))
      direct_only <- sum(unlist(result$direct) * unlist(direction))
      checks[[length(checks) + 1L]] <- data.frame(
        rho = rho, direction = j, numerical = numerical, analytical = analytical,
        absolute_error = abs(numerical - analytical),
        direct_only_error = abs(numerical - direct_only),
        kinks = result$kinks
      )
    }
  }
  do.call(rbind, checks)
}

.replay_pilot <- function(root, sim_ids) {
  rows <- list()
  for (sim_id in sim_ids) {
    directory <- file.path(root, sprintf("seed_%06d", sim_id))
    bundle <- readRDS(file.path(directory, "artifacts.rds"))
    for (rho in RHO_VALUES) {
      fit <- .extract_tate_fit(bundle, rho)
      row <- .extract_tate_row(bundle, rho)
      answer <- .weight_layer_influence(fit)
      stopifnot(abs(row$estimate - answer$estimate) < 1e-10,
                abs(row$se^2 - answer$fixed_variance) < 1e-10)
      se <- sqrt(answer$weight_layer_variance)
      rows[[length(rows) + 1L]] <- data.frame(
        sim_id = sim_id, rho = rho, estimate = answer$estimate, truth = row$truth,
        fixed_weight_se = row$se, weight_layer_se = se,
        fixed_weight_coverage = abs(answer$estimate - row$truth) <= stats::qnorm(0.975) * row$se,
        weight_layer_coverage = abs(answer$estimate - row$truth) <= stats::qnorm(0.975) * se,
        indirect_variance = answer$indirect_variance,
        direct_indirect_cross_term = answer$direct_indirect_cross_term,
        kinks = answer$kinks, fold_source_cells = answer$fold_source_cells
      )
    }
    rm(bundle)
    if (sim_id %% 10L == 0L) cat("replayed seeds through", sim_id, "\n")
  }
  do.call(rbind, rows)
}

.summarize_replay <- function(rows) {
  do.call(rbind, lapply(split(rows, rows$rho), function(x) {
    error <- x$estimate - x$truth
    data.frame(
      rho = x$rho[1], n = nrow(x), bias = mean(error), empirical_sd = stats::sd(error),
      rmse = sqrt(mean(error^2)),
      mean_fixed_se = mean(x$fixed_weight_se), mean_weight_layer_se = mean(x$weight_layer_se),
      fixed_se_to_sd = mean(x$fixed_weight_se) / stats::sd(error),
      weight_layer_se_to_sd = mean(x$weight_layer_se) / stats::sd(error),
      fixed_weight_coverage = mean(x$fixed_weight_coverage),
      weight_layer_coverage = mean(x$weight_layer_coverage),
      coverage_mcse = sqrt(mean(x$weight_layer_coverage) * (1 - mean(x$weight_layer_coverage)) / nrow(x)),
      negative_cross_terms = sum(x$direct_indirect_cross_term < 0),
      kink_fraction = sum(x$kinks) / sum(x$fold_source_cells)
    )
  }))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 1L) stop("usage: weight_layer.R OUTPUT_DIRECTORY")
  output <- file.path(args[1], "v1")
  if (file.exists(output)) stop("output already exists: ", output)
  checks <- .directional_checks(file.path(PILOT_ROOT, "seed_010013"))
  cat("maximum directional error:", max(checks$absolute_error), "\n")
  cat("maximum direct-only error:", max(checks$direct_only_error), "\n")
  rows <- .replay_pilot(PILOT_ROOT, 10001:10100)
  summary <- .summarize_replay(rows)
  print(summary, row.names = FALSE, digits = 4)
  code_files <- c("diagnosis/weight_layer/weight_layer.R",
                  "diagnosis/tate_common_weight/quadratic_bias_weights.R",
                  "diagnosis/tate_common_weight/quadratic_weight_influence.R")
  roce_write_atomic_directory(output, function(stage) {
    write.csv(checks, file.path(stage, "directional_checks.csv"), row.names = FALSE)
    write.csv(rows, file.path(stage, "replay_rows.csv"), row.names = FALSE)
    write.csv(summary, file.path(stage, "replay_summary.csv"), row.names = FALSE)
    write.csv(data.frame(path = code_files, sha256 = vapply(code_files, roce_sha256_file, "")),
              file.path(stage, "code_hashes.csv"), row.names = FALSE)
    writeLines(c(
      "task=weight_layer", "history_entry=0001",
      paste0("maximum_directional_error=", max(checks$absolute_error)),
      paste0("maximum_direct_only_error=", max(checks$direct_only_error)),
      "weight_rule=soft_threshold_production", "nuisance_fits_differentiated=FALSE",
      "resampling_draws=0", "new_mc_replications=0"
    ), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "weight layer prototype")
}

if (sys.nframe() == 0L) main()
