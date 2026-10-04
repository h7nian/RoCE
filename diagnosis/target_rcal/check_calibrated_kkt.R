#!/usr/bin/env Rscript
# Independently reconstruct gradients from saved data and coefficients, without
# calling the package's gradient, loss, or fitting implementations.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: check_calibrated_kkt.R FULL_FITS_RDS OUTPUT_CSV")
artifact <- readRDS(args[1L])
clamp <- function(x, radius) pmax(pmin(x, radius), -radius)
kkt_residual <- function(gradient, coefficients, lambda) {
  slopes <- coefficients[-1L]
  residual <- ifelse(abs(slopes) > 1e-7,
                     abs(gradient[-1L] + lambda * sign(slopes)),
                     pmax(abs(gradient[-1L]) - lambda, 0))
  max(abs(gradient[1L]), residual)
}
rows <- list()
for (arm_name in c("mu1", "mu0")) {
  arm <- artifact$reference$arm_results[[arm_name]]
  A_val <- if (arm_name == "mu1") 1L else 0L
  target <- artifact$data_split$t
  for (k1 in seq_len(arm$n_folds)) for (site in names(arm$fold_results[[k1]]$source_results)) {
    fit <- arm$fold_results[[k1]]$source_results[[site]]
    source <- artifact$data_split[[site]]
    gamma <- as.numeric(fit$gamma_s)
    alpha <- as.numeric(fit$alpha_ts)
    target_moments <- list()
    gamma_scores <- alpha_scores <- list()
    n_calibration <- 0L
    for (key in names(fit$per_k2_alpha)) {
      k2 <- as.integer(sub("k2_", "", key, fixed = TRUE))
      target_rows <- artifact$folds$target_folds[[k2]]$original_idx
      source_rows <- artifact$folds$source_folds[[site]][[k2]]$original_idx
      target_outcome_design <- cbind(1, target$W_outcome[target_rows, , drop = FALSE])
      target_site_design <- cbind(1, target$Z_site[target_rows, , drop = FALSE])
      outcome_design <- cbind(1, source$W_outcome[source_rows, , drop = FALSE])
      site_design <- cbind(1, source$Z_site[source_rows, , drop = FALSE])
      initial_alpha <- as.numeric(fit$per_k2_alpha[[key]])
      initial_gamma <- as.numeric(fit$per_k2_gamma[[key]])
      target_initial_mean <- plogis(clamp(drop(target_outcome_design %*% initial_alpha), arm$M_tau))
      source_initial_mean <- plogis(clamp(drop(outcome_design %*% initial_alpha), arm$M_tau))
      target_moments[[key]] <- colMeans(target_site_design *
                                         (target_initial_mean * (1 - target_initial_mean)))
      indicator <- as.numeric(source$A[source_rows] == A_val)
      gamma_weight <- exp(-clamp(drop(site_design %*% gamma), arm$M_tau))
      gamma_scores[[key]] <- colSums(site_design *
        (indicator * source_initial_mean * (1 - source_initial_mean) * gamma_weight))
      alpha_weight <- exp(-clamp(drop(site_design %*% initial_gamma), arm$M_tau))
      fitted_mean <- plogis(drop(outcome_design %*% alpha))
      alpha_scores[[key]] <- colSums(outcome_design *
        (indicator * alpha_weight * (fitted_mean - source$Y[source_rows])))
      n_calibration <- n_calibration + length(source_rows)
    }
    gamma_gradient <- Reduce(`+`, target_moments) / length(target_moments) -
      Reduce(`+`, gamma_scores) / n_calibration
    alpha_gradient <- Reduce(`+`, alpha_scores) / n_calibration
    for (model in c("calibrated_tilting", "calibrated_outcome")) {
      coefficients <- if (model == "calibrated_tilting") fit$gamma_s else fit$alpha_ts
      gradient <- if (model == "calibrated_tilting") gamma_gradient else alpha_gradient
      lambda <- as.numeric(attr(coefficients, "lambda_used"))
      stopifnot(length(lambda) == 1L, is.finite(lambda), lambda > 0)
      residual <- kkt_residual(gradient, as.numeric(coefficients), lambda)
      rows[[length(rows) + 1L]] <- data.frame(
        A_val, k1, site, model, n_calibration, lambda,
        reported_converged = isTRUE(attr(coefficients, "converged")),
        max_kkt_residual = residual, residual_over_lambda = residual / lambda
      )
    }
  }
}
rows <- do.call(rbind, rows)
write.csv(rows, args[2L], row.names = FALSE)
print(aggregate(cbind(max_kkt_residual, residual_over_lambda) ~ model + A_val,
                rows, max), row.names = FALSE)
