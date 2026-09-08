#!/usr/bin/env Rscript

# No-refit p=100 matrix/derivative gate. This does not publish a new estimator,
# choose a production penalty, calculate a new CI, or add MC replications.
.saved_projection_block <- function(data, index, arm) {
  list(W = cbind(1, data$W_outcome[index, , drop = FALSE]),
       Z = cbind(1, data$Z_site[index, , drop = FALSE]),
       A = as.numeric(data$A[index] == arm), Y = data$Y[index], index = index)
}

.saved_projection_state <- function(bundle, arm = 1L) {
  if (!is.numeric(arm) || length(arm) != 1L || is.na(arm) || !arm %in% 0:1) {
    stop("arm must be 0 or 1")
  }
  arm <- as.integer(arm)
  artifact <- bundle$group_result$artifacts[["0"]]
  fit <- artifact$direct_tate_results$one_round_crossfit
  data <- artifact$data_split
  source_fit <- fit$arm_results[[paste0("mu", arm)]]$fold_results[[1L]]$source_results$s1
  info <- fit$intermediates$fold_info
  source_index <- match("s1", names(data)[names(data) != "t"])
  target_ids <- lapply(info, `[[`, "target_idx")
  source_ids <- lapply(info, function(x) x$source_idx[[source_index]])
  keys <- paste0("k2_", 2:5)
  stopifnot(identical(names(source_fit$per_k2_alpha), keys),
            identical(names(source_fit$per_k2_gamma), keys),
            identical(sort(unlist(target_ids, use.names = FALSE)), seq_len(data$t$n)),
            identical(sort(unlist(source_ids, use.names = FALSE)), seq_len(data$s1$n)))
  pieces <- setNames(lapply(2:5, function(k) {
    complement <- setdiff(1:5, c(1L, k))
    list(target_training = .saved_projection_block(data$t, unlist(target_ids[complement]), arm),
         source_training = .saved_projection_block(data$s1, unlist(source_ids[complement]), arm),
         target_calibration = .saved_projection_block(data$t, target_ids[[k]], arm),
         source_calibration = .saved_projection_block(data$s1, source_ids[[k]], arm))
  }), keys)
  values <- c(setNames(lapply(source_fit$per_k2_alpha, as.numeric), paste0("alpha_initial_", keys)),
              setNames(lapply(source_fit$per_k2_gamma, as.numeric), paste0("gamma_initial_", keys)),
              list(alpha_final = as.numeric(source_fit$alpha_ts),
                   gamma_final = as.numeric(source_fit$gamma_s)))
  sizes <- lengths(values)
  ends <- cumsum(sizes)
  indices <- setNames(lapply(seq_along(sizes), function(i)
    seq.int(ends[i]-sizes[i]+1L, ends[i])), names(values))
  stopifnot(length(values$alpha_final) == ncol(data$t$W_outcome)+1L,
    length(values$gamma_final) == ncol(data$t$Z_site)+1L,
    ncol(data$t$W_outcome) == ncol(data$s1$W_outcome),
    ncol(data$t$Z_site) == ncol(data$s1$Z_site))
  list(pieces = pieces, keys = keys, indices = indices, theta = unlist(values, use.names = FALSE),
       target = .saved_projection_block(data$t, unlist(target_ids[2:5]), arm),
       source = .saved_projection_block(data$s1, unlist(source_ids[2:5]), arm),
       M_fit = 5, M_inference = 5, logistic_clip = RoCE:::LOGISTIC_CLIP,
       final_penalties = c(alpha_final = attr(source_fit$alpha_ts, "lambda_used"),
                           gamma_final = attr(source_fit$gamma_s, "lambda_used")))
}

.evaluate_saved_projection <- function(state, theta, derivatives = FALSE) {
  index <- state$indices
  alpha <- theta[index$alpha_final]; gamma <- theta[index$gamma_final]
  clip <- function(x, bound) pmax(-bound, pmin(bound, x))
  gram <- function(left, weight, right) crossprod(left, weight*right)
  count_source <- nrow(state$source$W)
  count_target <- nrow(state$target$W)
  moments <- numeric(length(theta))
  J <- if (derivatives) matrix(0, length(theta), length(theta)) else NULL
  if (derivatives) {
    H_alpha_final <- matrix(0, length(alpha), length(alpha))
    H_gamma_final <- matrix(0, length(gamma), length(gamma))
  }
  boundary_distance <- Inf
  truncation_counts <- c(initial_outcome = 0L, initial_tilting = 0L, final_tilting = 0L)
  for (key in state$keys) {
    source_denominator <- if (is.null(state$source_calibration_denominators)) count_source else
      state$source_calibration_denominators[[key]]
    stopifnot(length(source_denominator) == 1L, is.finite(source_denominator), source_denominator > 0)
    piece <- state$pieces[[key]]
    at <- index[[paste0("alpha_initial_", key)]]
    gt <- index[[paste0("gamma_initial_", key)]]
    a0 <- theta[at]; g0 <- theta[gt]
    tt <- piece$target_training; st <- piece$source_training
    tc <- piece$target_calibration; sc <- piece$source_calibration
    mt <- plogis(drop(tt$W %*% a0))
    wt <- exp(-drop(st$Z %*% g0))
    moments[at] <- colMeans(tt$W*(tt$A*(tt$Y-mt)))
    moments[gt] <- colMeans(tt$Z)-colMeans(st$Z*(st$A*wt))
    alpha0_t <- drop(tc$W %*% a0); alpha0_s <- drop(sc$W %*% a0)
    gamma0_s <- drop(sc$Z %*% g0)
    # The target gradient summary is not M_fit-truncated in the actual code.
    # Only the source-side derivative plug-in in the calibrated loss is.
    stopifnot(max(abs(alpha0_t)) < state$logistic_clip-1)
    u_t <- plogis(alpha0_t)
    u_s <- plogis(clip(alpha0_s, state$M_fit))
    d_t <- u_t*(1-u_t); d_s <- u_s*(1-u_s)
    w0 <- exp(-clip(gamma0_s, state$M_fit))
    w <- exp(-drop(sc$Z %*% gamma))
    m <- plogis(drop(sc$W %*% alpha))
    moments[index$alpha_final] <- moments[index$alpha_final]+
      drop(crossprod(sc$W, sc$A*w0*(sc$Y-m)))/source_denominator
    moments[index$gamma_final] <- moments[index$gamma_final]+
      colMeans(tc$Z*d_t)/length(state$keys)-
      drop(crossprod(sc$Z, sc$A*w*d_s))/source_denominator
    boundary_distance <- min(boundary_distance,
      abs(abs(c(alpha0_s, gamma0_s))-state$M_fit))
    truncation_counts["initial_outcome"] <- truncation_counts["initial_outcome"]+
      sum(abs(alpha0_s) > state$M_fit)
    truncation_counts["initial_tilting"] <- truncation_counts["initial_tilting"]+
      sum(sc$A == 1 & abs(gamma0_s) > state$M_fit)
    if (derivatives) {
      J[at, at] <- -gram(tt$W, tt$A*mt*(1-mt), tt$W)/nrow(tt$W)
      J[gt, gt] <- gram(st$Z, st$A*wt, st$Z)/nrow(st$Z)
      H_alpha_final <- H_alpha_final+gram(sc$W, sc$A*w0*m*(1-m), sc$W)/source_denominator
      H_gamma_final <- H_gamma_final+gram(sc$Z, sc$A*w*d_s, sc$Z)/source_denominator
      J[index$alpha_final, gt] <- -gram(sc$W,
        sc$A*w0*(sc$Y-m)*(abs(gamma0_s) < state$M_fit), sc$Z)/source_denominator
      J[index$gamma_final, at] <- gram(tc$Z,
        d_t*(1-2*u_t), tc$W)/
        (length(state$keys)*nrow(tc$W))-
        gram(sc$Z, sc$A*w*d_s*(1-2*u_s)*(abs(alpha0_s) < state$M_fit), sc$W)/source_denominator
    }
  }
  target_prediction <- plogis(drop(state$target$W %*% alpha))
  source_prediction <- plogis(drop(state$source$W %*% alpha))
  final_log_weight <- drop(state$source$Z %*% gamma)
  final_weight <- exp(-clip(final_log_weight, state$M_inference))
  score <- mean(target_prediction)+mean(state$source$A*final_weight*(state$source$Y-source_prediction))
  boundary_distance <- min(boundary_distance, abs(abs(final_log_weight)-state$M_inference))
  truncation_counts["final_tilting"] <- sum(state$source$A == 1 & abs(final_log_weight) > state$M_inference)
  D <- numeric(length(theta))
  if (derivatives) {
    J[index$alpha_final, index$alpha_final] <- -H_alpha_final
    J[index$gamma_final, index$gamma_final] <- H_gamma_final
    D[index$alpha_final] <- colMeans(state$target$W*(target_prediction*(1-target_prediction)))-
      colMeans(state$source$W*(state$source$A*final_weight*source_prediction*(1-source_prediction)))
    D[index$gamma_final] <- -colMeans(state$source$Z*(state$source$A*final_weight*
      (state$source$Y-source_prediction)*(abs(final_log_weight) < state$M_inference)))
  }
  stopifnot(all(is.finite(moments)), is.finite(score), boundary_distance > 1e-5)
  list(moments = moments, jacobian = J, score_gradient = D, score = score,
       boundary_distance = boundary_distance, truncation_counts = truncation_counts)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) %in% c(3L, 4L)) stop("usage: probe_saved_nuisance_projection.R V19_LIBRARY SEED_BUNDLE OUTPUT [ARM]")
  arm <- if (length(args) == 4L) suppressWarnings(as.numeric(args[4])) else 1L
  lib <- normalizePath(args[1], mustWork = TRUE)
  directory <- normalizePath(args[2], mustWork = TRUE)
  output <- args[3]
  if (file.exists(output)) stop("matrix probe output already exists")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  source("diagnosis/tate_common_weight/sparse_moment_projection.R")
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"),
    .weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
    roce_sha256_file(file.path(directory, "sha256.txt")) ==
      "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  lines <- readLines(file.path(directory, "sha256.txt"))
  expected <- substr(lines[substring(lines, 67L) == "artifacts.rds"], 1L, 64L)
  stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, "artifacts.rds")) == expected)
  state <- .saved_projection_state(readRDS(file.path(directory, "artifacts.rds")), arm)
  evaluated <- .evaluate_saved_projection(state, state$theta, TRUE)
  native <- RoCE:::calculate_correction_term_cpp(state$source$Z[, -1, drop = FALSE],
    state$source$A, state$source$Y, state$theta[state$indices$gamma_final],
    state$theta[state$indices$alpha_final], state$source$W[, -1, drop = FALSE], 5, 1L, 1L, 1L)
  native_score <- mean(RoCE:::predict_glm_cpp(state$target$W[, -1, drop = FALSE],
    state$theta[state$indices$alpha_final], 1L, 1L))+native$delta_ts
  stopifnot(abs(native_score-evaluated$score) < 1e-12, !native$clip_diagnostics$any_safety_clipped)
  checks <- do.call(rbind, lapply(1:12, function(k) {
    direction <- sin(seq_along(state$theta)*k); direction <- direction/sqrt(sum(direction^2))
    epsilon <- 1e-6
    plus <- .evaluate_saved_projection(state, state$theta+epsilon*direction)
    minus <- .evaluate_saved_projection(state, state$theta-epsilon*direction)
    data.frame(direction = k,
      jacobian_error = max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(evaluated$jacobian %*% direction))),
      score_gradient_error = abs((plus$score-minus$score)/(2*epsilon)-sum(evaluated$score_gradient*direction)))
  }))
  stopifnot(max(checks$jacobian_error) < 1e-7, max(checks$score_gradient_error) < 1e-7)
  block_checks <- lapply(c("alpha_final", "gamma_final"), function(block) {
    index <- state$indices[[block]]
    loss_gradient <- if (block == "alpha_final") -evaluated$moments[index] else evaluated$moments[index]
    penalty <- c(0, rep(state$final_penalties[[block]], length(index)-1L))
    coefficients <- state$theta[index]; active <- coefficients != 0
    violation <- pmax(abs(loss_gradient)-penalty, 0)
    violation[active] <- abs(loss_gradient[active]+penalty[active]*sign(coefficients[active]))
    data.frame(block, maximum_native_kkt_error = max(violation))
  })
  block_checks <- do.call(rbind, block_checks)
  stopifnot(max(block_checks$maximum_native_kkt_error) < 1e-4)
  J <- evaluated$jacobian; D <- evaluated$score_gradient; index <- state$indices
  attempts <- lapply(c(1, .5, .1), function(fraction) {
    penalty <- fraction*max(abs(D))
    fits <- list(); coefficients <- numeric(length(D))
    fit_block <- function(name, rhs, sign) {
      pos <- index[[name]]
      fit <- .solve_sparse_moment_projection(sign*J[pos, pos, drop = FALSE], rhs,
                                              penalty, tolerance = 1e-8, max_iterations = 2000L)
      coefficients[pos] <<- fit$coefficients
      fits[[name]] <<- fit
    }
    failure <- tryCatch({
      fit_block("alpha_final", -D[index$alpha_final], -1)
      fit_block("gamma_final", D[index$gamma_final], 1)
      for (key in state$keys) {
        at <- paste0("alpha_initial_", key); gt <- paste0("gamma_initial_", key)
        fit_block(at, drop(crossprod(J[index$gamma_final, index[[at]], drop = FALSE], coefficients[index$gamma_final])), -1)
        fit_block(gt, -drop(crossprod(J[index$alpha_final, index[[gt]], drop = FALSE], coefficients[index$alpha_final])), 1)
      }
      NULL
    }, error = identity)
    residual <- max(abs(drop(t(J) %*% coefficients)-D))
    if (is.null(failure)) stopifnot(residual <= penalty+1e-5)
    list(summary = data.frame(fraction, penalty, solved = is.null(failure),
      completed_blocks = length(fits), adjoint_residual = residual,
      coefficient_l1 = sum(abs(coefficients)),
      failure = if (is.null(failure)) "" else conditionMessage(failure)),
      block_fits = fits, coefficients = coefficients, failure = failure)
  })
  summaries <- do.call(rbind, lapply(attempts, `[[`, "summary"))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(checks, file.path(stage, "directional_derivatives.csv"), row.names = FALSE)
    write.csv(block_checks, file.path(stage, "native_calibration_kkt.csv"), row.names = FALSE)
    write.csv(summaries, file.path(stage, "solver_attempts.csv"), row.names = FALSE)
    saveRDS(list(state = state, evaluated = evaluated, attempts = attempts), file.path(stage, "matrix_probe.rds"))
    writeLines(c("saved_p100_matrix_probe=complete", "sim_id=10013", "rho=0", "outer_fold=1", "source=s1", paste0("arm=", arm),
      paste0("nuisance_dimension=", length(state$theta)),
      paste0("maximum_jacobian_error=", max(checks$jacobian_error)),
      paste0("maximum_score_gradient_error=", max(checks$score_gradient_error)),
      paste0("native_score_error=", abs(native_score-evaluated$score)),
      paste0("minimum_truncation_boundary_distance=", evaluated$boundary_distance),
      paste0(names(evaluated$truncation_counts), "_truncations=", evaluated$truncation_counts),
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R")),
      paste0("solver_sha256=", roce_sha256_file("diagnosis/tate_common_weight/sparse_moment_projection.R")),
      paste0("artifact_sha256=", expected),
      "penalty_policy_selected=FALSE", "new_mc_replications=0", "nuisance_refits=0",
      "primary_estimator_se_ci_changed=FALSE", "inference_validated=FALSE"), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "saved nuisance projection matrix gate")
  print(checks); print(block_checks); print(summaries)
}

if (sys.nframe() == 0L) main()
