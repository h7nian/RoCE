#!/usr/bin/env Rscript

.target_tate_system <- function(block, theta, limits, derivatives = FALSE) {
  Z <- block$Z; W <- block$W; A <- block$A; Y <- block$Y
  pz <- ncol(Z); pw <- ncol(W); n <- nrow(W)
  index <- list(propensity = seq_len(pz), outcome_treated = pz+seq_len(pw),
                outcome_control = pz+pw+seq_len(pw))
  raw_p <- plogis(drop(Z %*% theta[index$propensity]))
  raw_m1 <- plogis(drop(W %*% theta[index$outcome_treated]))
  raw_m0 <- plogis(drop(W %*% theta[index$outcome_control]))
  clip <- function(x, lower, upper) pmax(lower, pmin(upper, x))
  p <- clip(raw_p, limits$ps[1], limits$ps[2])
  m1 <- clip(raw_m1, limits$outcome[1], limits$outcome[2])
  m0 <- clip(raw_m0, limits$outcome[1], limits$outcome[2])
  score <- m1-m0+A/p*(Y-m1)-(1-A)/(1-p)*(Y-m0)
  moment_rows <- cbind(Z*(A-raw_p), W*(A*(Y-raw_m1)), W*((1-A)*(Y-raw_m0)))
  moments <- colMeans(moment_rows)
  J <- D <- NULL
  if (derivatives) {
    J <- matrix(0, length(theta), length(theta))
    D <- numeric(length(theta))
    J[index$propensity, index$propensity] <- -crossprod(Z, Z*(raw_p*(1-raw_p)))/n
    J[index$outcome_treated, index$outcome_treated] <- -crossprod(W, W*(A*raw_m1*(1-raw_m1)))/n
    J[index$outcome_control, index$outcome_control] <- -crossprod(W, W*((1-A)*raw_m0*(1-raw_m0)))/n
    dp <- raw_p*(1-raw_p)*(raw_p > limits$ps[1] & raw_p < limits$ps[2])
    dm1 <- raw_m1*(1-raw_m1)*(raw_m1 > limits$outcome[1] & raw_m1 < limits$outcome[2])
    dm0 <- raw_m0*(1-raw_m0)*(raw_m0 > limits$outcome[1] & raw_m0 < limits$outcome[2])
    D[index$propensity] <- colMeans(Z*(-A*(Y-m1)/p^2-(1-A)*(Y-m0)/(1-p)^2)*dp)
    D[index$outcome_treated] <- colMeans(W*((1-A/p)*dm1))
    D[index$outcome_control] <- -colMeans(W*((1-(1-A)/(1-p))*dm0))
  }
  stopifnot(all(is.finite(score)), all(is.finite(moments)))
  list(score = score, estimate = mean(score), moment_rows = moment_rows, moments = moments,
       jacobian = J, gradient = D, index = index, predictions = list(p = p, m1 = m1, m0 = m0),
       clipping = c(propensity = sum(p != raw_p), outcome_treated = sum(m1 != raw_m1),
                    outcome_control = sum(m0 != raw_m0)))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) stop("usage: probe_target_nuisance_system.R V19_LIBRARY FIT_STATE_BUNDLE SEED_BUNDLE OUTPUT")
  lib <- normalizePath(args[1], mustWork = TRUE)
  states_root <- normalizePath(args[2], mustWork = TRUE)
  seed_root <- normalizePath(args[3], mustWork = TRUE); output <- args[4]
  if (file.exists(output)) stop("target-system output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/sparse_moment_projection.R")
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"),
    .weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(seed_root, "sha256.txt")) ==
      "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  checked_rds <- function(directory, name) {
    lines <- readLines(file.path(directory, "sha256.txt"))
    expected <- substr(lines[substring(lines, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, name)) == expected)
    readRDS(file.path(directory, name))
  }
  fits <- checked_rds(states_root, "target_cv_fit_states.rds")
  bundle <- checked_rds(seed_root, "artifacts.rds")
  artifact <- bundle$group_result$artifacts[["0"]]
  fitted_tate <- artifact$direct_tate_results$one_round_crossfit
  ids <- lapply(fitted_tate$intermediates$fold_info, `[[`, "target_idx")
  train_ids <- unlist(ids[2:5], use.names = FALSE); eval_ids <- ids[[1L]]
  stopifnot(!any(train_ids %in% eval_ids), !anyDuplicated(c(train_ids, eval_ids)))
  make_block <- function(index) list(
    W = cbind(1, artifact$data_split$t$W_outcome[index, , drop = FALSE]),
    Z = cbind(1, artifact$data_split$t$Z_site[index, , drop = FALSE]),
    A = artifact$data_split$t$A[index], Y = artifact$data_split$t$Y[index], ids = index)
  training <- make_block(train_ids); evaluation <- make_block(eval_ids)
  keys <- c(propensity = "PS_1_1_NA", outcome_treated = "OR_1_1_NA", outcome_control = "OR_0_1_NA")
  stopifnot(all(keys %in% names(fits)))
  theta <- unlist(lapply(fits[keys], function(x) as.numeric(stats::coef(x, s = "lambda.min"))), use.names = FALSE)
  limits <- list(ps = c(RoCE:::PROP_SCORE_LOWER, RoCE:::PROP_SCORE_UPPER),
                 outcome = c(RoCE:::OUTCOME_PRED_LOWER, RoCE:::OUTCOME_PRED_UPPER))
  system <- .target_tate_system(training, theta, limits, TRUE)
  heldout <- .target_tate_system(evaluation, theta, limits)
  saved1 <- fitted_tate$arm_results$mu1$fold_results[[1L]]$target_only
  saved0 <- fitted_tate$arm_results$mu0$fold_results[[1L]]$target_only
  prediction_error <- max(abs(heldout$predictions$p-saved1$prop_scores),
    abs(heldout$predictions$p-saved0$prop_scores), abs(heldout$predictions$m1-saved1$m_pred),
    abs(heldout$predictions$m0-saved0$m_pred))
  estimate_error <- abs(heldout$estimate-(saved1$estimate-saved0$estimate))
  stopifnot(length(theta) == 603L, prediction_error < 1e-12, estimate_error < 1e-12)
  checks <- do.call(rbind, lapply(1:12, function(k) {
    direction <- sin(seq_along(theta)*k); direction <- direction/sqrt(sum(direction^2))
    epsilon <- 1e-6
    plus <- .target_tate_system(training, theta+epsilon*direction, limits)
    minus <- .target_tate_system(training, theta-epsilon*direction, limits)
    data.frame(direction = k,
      jacobian_error = max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(system$jacobian %*% direction))),
      score_gradient_error = abs((plus$estimate-minus$estimate)/(2*epsilon)-sum(system$gradient*direction)))
  }))
  stopifnot(max(checks$jacobian_error) < 1e-7, max(checks$score_gradient_error) < 1e-7)
  kkt <- do.call(rbind, lapply(names(keys), function(block) {
    index <- system$index[[block]]
    selected <- if (block == "propensity") rep(TRUE, nrow(training$W)) else
      training$A == as.integer(block == "outcome_treated")
    X <- if (block == "propensity") training$Z[, -1, drop = FALSE] else training$W[selected, -1, drop = FALSE]
    feature_sd <- sqrt(colMeans(sweep(X, 2L, colMeans(X), "-")^2))
    penalty <- c(0, fits[[keys[[block]]]]$lambda.min*feature_sd*mean(selected))
    gradient <- -system$moments[index]; coefficients <- theta[index]
    active <- coefficients != 0
    residual <- pmax(abs(gradient)-penalty, 0)
    residual[active] <- abs(gradient[active]+penalty[active]*sign(coefficients[active]))
    data.frame(block, maximum_standardized_glmnet_kkt_error = max(residual))
  }))
  stopifnot(max(kkt$maximum_standardized_glmnet_kkt_error) < 1e-3)
  attempts <- lapply(c(1, .5, .1), function(fraction) {
    penalty <- fraction*max(abs(system$gradient))
    a <- numeric(length(theta)); block_fits <- list()
    failure <- tryCatch({
      for (block in names(system$index)) {
        index <- system$index[[block]]
        answer <- .solve_sparse_moment_projection(-system$jacobian[index, index, drop = FALSE],
          -system$gradient[index], penalty, tolerance = 1e-8, max_iterations = 2000L)
        a[index] <- answer$coefficients; block_fits[[block]] <- answer
      }
      NULL
    }, error = identity)
    corrected <- if (is.null(failure)) heldout$score-drop(heldout$moment_rows %*% a) else rep(NA_real_, length(heldout$score))
    list(summary = data.frame(fraction, penalty, solved = is.null(failure),
      coefficient_l1 = sum(abs(a)), original_target_fold_tate = heldout$estimate,
      adjusted_target_fold_tate = mean(corrected), original_score_sd = sd(heldout$score),
      adjusted_score_sd = sd(corrected), adjoint_residual = max(abs(drop(t(system$jacobian) %*% a)-system$gradient)),
      failure = if (is.null(failure)) "" else conditionMessage(failure)),
      coefficients = a, block_fits = block_fits, failure = failure, heldout_score = corrected)
  })
  summaries <- do.call(rbind, lapply(attempts, `[[`, "summary"))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(checks, file.path(stage, "directional_checks.csv"), row.names = FALSE)
    write.csv(kkt, file.path(stage, "glmnet_kkt_checks.csv"), row.names = FALSE)
    write.csv(summaries, file.path(stage, "projection_probes.csv"), row.names = FALSE)
    saveRDS(list(training = training, evaluation = evaluation, theta = theta, limits = limits,
      system = system, attempts = attempts), file.path(stage, "target_system.rds"))
    writeLines(c("target_nuisance_system_probe=complete", "sim_id=10013", "outer_fold=1",
      "nuisance_dimension=603", "shared_propensity_block_count=1",
      paste0("maximum_prediction_error=", prediction_error), paste0("fold_tate_error=", estimate_error),
      paste0("maximum_jacobian_error=", max(checks$jacobian_error)),
      paste0("maximum_score_gradient_error=", max(checks$score_gradient_error)),
      paste0("training_", names(system$clipping), "_clipped=", system$clipping),
      paste0("evaluation_", names(heldout$clipping), "_clipped=", heldout$clipping),
      "holdout_used_for_penalty_selection=FALSE", "new_mc_replications=0", "resampling_draws=0",
      "nuisance_refits_by_this_probe=0", "inference_validated=FALSE",
      paste0("source_states_manifest_sha256=", roce_sha256_file(file.path(states_root, "sha256.txt"))),
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/probe_target_nuisance_system.R"))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "target nuisance system probe")
  print(kkt, row.names = FALSE); print(summaries, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
