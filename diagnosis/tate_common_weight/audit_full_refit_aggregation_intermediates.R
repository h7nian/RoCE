#!/usr/bin/env Rscript

.afri_max_error <- function(x, y) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (length(x) != length(y) || any(!is.finite(c(x, y)))) return(Inf)
  if (length(x) == 0L) return(0)
  max(abs(x - y))
}

.afri_scalar <- function(x, name, type = c("numeric", "logical", "character")) {
  type <- match.arg(type)
  valid_type <- switch(type, numeric = is.numeric(x), logical = is.logical(x),
                       character = is.character(x))
  if (!valid_type || length(x) != 1L || !is.null(dim(x)) || is.object(x) ||
      is.na(x) ||
      (type == "numeric" && !is.finite(x)) ||
      (type == "character" && !nzchar(x))) {
    stop(name, " must be one nonmissing ", type, " scalar.", call. = FALSE)
  }
  x
}

.afri_validate_draw_pair <- function(csv, saved) {
  if (!is.list(csv) || length(csv) != 2L || !is.list(saved) ||
      length(saved) != 2L || any(vapply(csv, nrow, integer(1L)) != 1L)) {
    stop("draw pair must contain two one-row result records.", call. = FALSE)
  }
  required <- c("sim_id", "rho", "draw_id", "identity", "status")
  if (any(vapply(csv, function(row) !all(required %in% names(row)), logical(1L)))) {
    stop("draw result is missing pair identity fields.", call. = FALSE)
  }
  sim_id <- vapply(csv, function(row) .afri_scalar(row$sim_id, "sim_id"), numeric(1L))
  rho <- vapply(csv, function(row) .afri_scalar(row$rho, "rho"), numeric(1L))
  draw_id <- vapply(csv, function(row) .afri_scalar(row$draw_id, "draw_id"), numeric(1L))
  identity <- vapply(csv, function(row) .afri_scalar(row$identity, "identity", "logical"), logical(1L))
  status <- vapply(csv, function(row) .afri_scalar(row$status, "status", "character"), character(1L))
  source_bundle <- vapply(saved, function(object) {
    .afri_scalar(object$source_bundle, "source_bundle", "character")
  }, character(1L))
  for (i in seq_len(2L)) {
    summary <- saved[[i]]$summary
    if (!is.data.frame(summary) || nrow(summary) != 1L ||
        !all(required %in% names(summary)) ||
        .afri_scalar(summary$sim_id, "saved sim_id") != sim_id[[i]] ||
        .afri_scalar(summary$rho, "saved rho") != rho[[i]] ||
        .afri_scalar(summary$draw_id, "saved draw_id") != draw_id[[i]] ||
        !identical(.afri_scalar(summary$identity, "saved identity", "logical"),
                   identity[[i]]) ||
        .afri_scalar(summary$status, "saved status", "character") != status[[i]]) {
      stop("draw CSV and RDS summary identity fields disagree.", call. = FALSE)
    }
  }
  if (any(sim_id < 1 | sim_id != floor(sim_id) |
          sim_id > .Machine$integer.max) ||
      any(draw_id < 0 | draw_id != floor(draw_id) |
          draw_id > .Machine$integer.max) ||
      sim_id[[1L]] != sim_id[[2L]] || rho[[1L]] != rho[[2L]] ||
      source_bundle[[1L]] != source_bundle[[2L]] ||
      draw_id[[1L]] != 0 || !identity[[1L]] ||
      draw_id[[2L]] <= 0 || identity[[2L]] || any(status != "completed")) {
    stop("draw pair is not a completed same-simulation/rho/source identity-positive pair.",
         call. = FALSE)
  }
  list(sim_id = sim_id[[1L]], rho = rho[[1L]],
       source_bundle = source_bundle[[1L]])
}

.afri_select_source_fit <- function(source_saved, rho) {
  rho <- .afri_scalar(rho, "rho")
  rho_key <- format(rho, scientific = FALSE, trim = TRUE)
  artifact <- source_saved$group_result$artifacts[[rho_key]]
  fit <- artifact$direct_tate_results$one_round_crossfit
  if (is.null(fit)) {
    stop("source calibration bundle is missing rho artifact: ", rho_key,
         call. = FALSE)
  }
  fit
}

.afri_manual_fold <- function(fit, label, k) {
  v <- fit$intermediates$inner_fold_info[[k]]
  eta <- as.numeric(fit$fold_weights[k, ])
  nt <- as.numeric(v$n_t); ns <- as.numeric(v$n_s)
  N <- nt + sum(ns); K <- length(eta); S <- sum(eta)
  cross <- (as.matrix(v$C_cross) + t(as.matrix(v$C_cross))) / 2
  D <- v$V_t / nt + v$V_s / ns
  variance_manual <- (1 - S)^2 * v$V_ot / nt +
    sum(eta^2 * D) +
    2 * (1 - S) * sum(eta * v$C_ot) / nt +
    if (K > 1L) sum(outer(eta, eta) * cross) / nt else 0
  variance_api <- RoCE::calculate_aggregated_variance(
    eta, list(V_ot = v$V_ot, V_t = v$V_t, V_s = v$V_s),
    v$C_ot, list(n_t = nt, n_s = ns), v$C_cross
  )
  discrepancy_variance <- v$V_ot / nt + v$V_t / nt + v$V_s / ns -
    2 * v$C_ot / nt
  discrepancy_se <- sqrt(pmax(discrepancy_variance, RoCE:::VARIANCE_MIN))
  wald <- abs(v$avg_target_est - v$avg_source_est) / discrepancy_se
  penalty <- pmax(fit$fold_lambdas[[k]] * wald - 1, 0)

  H <- 2 * diag(D, K) + 2 * v$V_ot / nt * matrix(1, K, K) -
    2 / nt * (outer(v$C_ot, rep(1, K)) +
      outer(rep(1, K), v$C_ot)) + 2 / nt * cross
  eigenvalues <- eigen((H + t(H)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values
  margin <- 1e-10 * max(abs(eigenvalues))
  expected_ridge <- max(0, margin - min(eigenvalues))
  ridge <- as.numeric(fit$fold_weight_psd_ridge[[k]])
  dotC <- sum(eta * v$C_ot)
  gradient <- vapply(seq_len(K), function(j) {
    -2 * (1 - S) * v$V_ot / nt + 2 * eta[[j]] * D[[j]] +
      2 * v$C_ot[[j]] * (1 - S) / nt - 2 * dotC / nt +
      2 * (sum(cross[j, ] * eta) - cross[j, j] * eta[[j]]) / nt +
      ridge * eta[[j]]
  }, numeric(1L))
  kkt <- ifelse(
    abs(eta) > 1e-12,
    abs(N * gradient + penalty * sign(eta)),
    pmax(abs(N * gradient) - penalty, 0)
  )
  rerun <- RoCE::optimize_weights(
    v$avg_source_est,
    list(V_ot = v$V_ot, V_t = v$V_t, V_s = v$V_s),
    v$C_ot, list(n_t = nt, n_s = ns), fit$fold_lambdas[[k]],
    v$avg_target_est, v$C_cross, clip_weights = FALSE
  )
  objective <- N * variance_manual + N * ridge / 2 * sum(eta^2) +
    sum(penalty * abs(eta))
  data.frame(
    fit = label, fold = k, n_t = nt, n_s_min = min(ns), n_s_max = max(ns),
    lambda = fit$fold_lambdas[[k]],
    discrepancy_variance_min = min(discrepancy_variance),
    discrepancy_variance_max = max(discrepancy_variance),
    wald_reconstruction_error = .afri_max_error(
      wald, fit$fold_wald_statistics[k, ]
    ),
    penalty_reconstruction_error = .afri_max_error(
      penalty, fit$fold_penalty_coefficients[k, ]
    ),
    variance_api_reconstruction_error = abs(variance_manual - variance_api),
    normalized_objective = objective,
    normalized_kkt_residual_max = max(kkt),
    optimizer_rerun_weight_error = .afri_max_error(rerun, eta),
    psd_ridge = ridge, psd_ridge_reconstruction_error =
      abs(ridge - expected_ridge),
    minimum_signed_weight = min(eta), maximum_signed_weight = max(eta),
    target_anchor_weight = 1 - sum(eta),
    stringsAsFactors = FALSE
  )
}

.afri_fold_source_diagnostics <- function(fit, label, k) {
  v <- fit$intermediates$inner_fold_info[[k]]
  eta <- as.numeric(fit$fold_weights[k, ])
  nt <- as.numeric(v$n_t); ns <- as.numeric(v$n_s)
  N <- nt + sum(ns); K <- length(eta); S <- sum(eta)
  cross <- (as.matrix(v$C_cross) + t(as.matrix(v$C_cross))) / 2
  ridge <- as.numeric(fit$fold_weight_psd_ridge[[k]])
  D <- v$V_t / nt + v$V_s / ns
  dotC <- sum(eta * v$C_ot)
  discrepancy_variance <- v$V_ot / nt + D - 2 * v$C_ot / nt
  discrepancy_se <- sqrt(pmax(discrepancy_variance, RoCE:::VARIANCE_MIN))
  wald <- abs(v$avg_target_est - v$avg_source_est) / discrepancy_se
  penalty <- pmax(fit$fold_lambdas[[k]] * wald - 1, 0)
  gradient <- vapply(seq_len(K), function(j) {
    -2 * (1 - S) * v$V_ot / nt + 2 * eta[[j]] * D[[j]] +
      2 * v$C_ot[[j]] * (1 - S) / nt - 2 * dotC / nt +
      2 * (sum(cross[j, ] * eta) - cross[j, j] * eta[[j]]) / nt +
      ridge * eta[[j]]
  }, numeric(1L))
  kkt <- ifelse(
    abs(eta) > 1e-12,
    abs(N * gradient + penalty * sign(eta)),
    pmax(abs(N * gradient) - penalty, 0)
  )
  data.frame(
    fit = label, fold = k, source = names(fit$weights),
    avg_target_estimate = v$avg_target_est,
    avg_source_estimate = as.numeric(v$avg_source_est),
    discrepancy = as.numeric(v$avg_source_est - v$avg_target_est),
    V_ot = v$V_ot, V_t = as.numeric(v$V_t), V_s = as.numeric(v$V_s),
    C_ot = as.numeric(v$C_ot), n_t = nt, n_s = ns,
    discrepancy_variance = discrepancy_variance,
    discrepancy_se = discrepancy_se, reconstructed_wald = wald,
    stored_wald = as.numeric(fit$fold_wald_statistics[k, ]),
    reconstructed_penalty = penalty,
    stored_penalty = as.numeric(fit$fold_penalty_coefficients[k, ]),
    signed_weight = eta, smooth_gradient = gradient,
    normalized_kkt_residual = kkt,
    stringsAsFactors = FALSE
  )
}

.afri_compare_combined <- function(observed, rebuilt) {
  fields <- c("V_ot", "V_t", "V_s", "C_ot", "C_cross",
              "fold_target_estimate", "fold_source_estimates",
              "avg_target_est", "avg_source_est", "mu_pred_ts", "delta_ts",
              "varphi_ot")
  available <- fields[fields %in% names(observed) & fields %in% names(rebuilt)]
  errors <- vapply(available, function(name)
    .afri_max_error(observed[[name]], rebuilt[[name]]), numeric(1L))
  for (name in intersect(c("zeta_components", "xi_components"),
                         intersect(names(observed), names(rebuilt)))) {
    errors <- c(errors, stats::setNames(max(vapply(
      seq_along(observed[[name]]), function(j) .afri_max_error(
        observed[[name]][[j]], rebuilt[[name]][[j]]
      ), numeric(1L))), name))
  }
  max(errors)
}

.afri_arm_audit <- function(fit, label) {
  K <- length(fit$weights)
  rows <- lapply(seq_len(fit$n_folds), function(k) {
    fold_rebuilt <- RoCE:::.combine_tate_fold_info(
      fit$arm_results$mu1$intermediates$fold_info[[k]],
      fit$arm_results$mu0$intermediates$fold_info[[k]], K
    )
    inner_rebuilt <- RoCE:::.combine_tate_inner_fold_info(
      fit$arm_results$mu1$intermediates$inner_fold_info[[k]],
      fit$arm_results$mu0$intermediates$inner_fold_info[[k]], K
    )
    fold_estimate_manual <- fit$intermediates$fold_info[[k]]$
      fold_target_estimate + sum(fit$fold_weights[k, ] * (
        fit$intermediates$fold_info[[k]]$fold_source_estimates -
          fit$intermediates$fold_info[[k]]$fold_target_estimate
      ))
    data.frame(
      fit = label, fold = k,
      fold_arm_contrast_error = .afri_compare_combined(
        fit$intermediates$fold_info[[k]], fold_rebuilt
      ),
      inner_arm_contrast_error = .afri_compare_combined(
        fit$intermediates$inner_fold_info[[k]], inner_rebuilt
      ),
      fold_aggregation_identity_error = abs(
        fold_estimate_manual - fit$fold_aggregated_estimates[[k]]
      ), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

.afri_tail_row <- function(values, fit, fold, site, component) {
  values <- as.numeric(values); centered <- values - mean(values)
  data.frame(
    fit = fit, fold = fold, site = site, component = component,
    n = length(values), mean = mean(values), sd = stats::sd(values),
    q95_abs_centered = unname(stats::quantile(abs(centered), .95)),
    q99_abs_centered = unname(stats::quantile(abs(centered), .99)),
    max_abs_centered = max(abs(centered)),
    stringsAsFactors = FALSE
  )
}

.afri_component_tails <- function(fit, label) {
  sources <- names(fit$weights); rows <- list(); index <- 0L
  for (k in seq_len(fit$n_folds)) {
    info <- fit$intermediates$fold_info[[k]]
    index <- index + 1L
    rows[[index]] <- .afri_tail_row(
      info$varphi_ot, label, k, "t", "target_varphi_tau"
    )
    for (j in seq_along(sources)) {
      index <- index + 1L
      rows[[index]] <- .afri_tail_row(
        info$zeta_components[[j]], label, k, sources[[j]], "target_zeta_tau"
      )
      index <- index + 1L
      rows[[index]] <- .afri_tail_row(
        info$xi_components[[j]], label, k, sources[[j]], "source_xi_tau"
      )
    }
  }
  do.call(rbind, rows)
}

.afri_final_tails <- function(fit, label) {
  sizes <- fit$intermediates$sample_sizes
  group_sizes <- c(t = sizes$n_t, sizes$n_source)
  starts <- cumsum(c(1L, head(group_sizes, -1L)))
  rows <- lapply(seq_along(group_sizes), function(i) {
    index <- starts[[i]]:(starts[[i]] + group_sizes[[i]] - 1L)
    values <- fit$all_phi_agg[index]; centered <- values - mean(values)
    data.frame(
      fit = label, site = names(group_sizes)[[i]], n = length(values),
      mean = mean(values), sd = stats::sd(values),
      q99_abs_centered = unname(stats::quantile(abs(centered), .99)),
      max_abs_centered = max(abs(centered)),
      variance_contribution = sum(centered^2) / length(fit$all_phi_agg)^2,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$variance_contribution_fraction <-
    result$variance_contribution / sum(result$variance_contribution)
  result
}

audit_full_refit_aggregation_main <- function(
    args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) stop(paste(
    "usage: audit_full_refit_aggregation_intermediates.R",
    "IDENTITY_DRAW_DIR DRAW1_DIR PACKAGE_LIB OUTPUT_DIR"), call. = FALSE)
  draw_dirs <- normalizePath(args[1:2], mustWork = TRUE)
  package_lib <- normalizePath(args[[3L]], mustWork = TRUE)
  output <- file.path(normalizePath(dirname(args[[4L]]), mustWork = FALSE),
                      basename(args[[4L]]))
  if (file.exists(output) || dir.exists(output))
    stop("audit output already exists.", call. = FALSE)
  .libPaths(c(package_lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  if (!identical(dirname(normalizePath(find.package("RoCE"))), package_lib))
    stop("RoCE is not loaded from PACKAGE_LIB.", call. = FALSE)
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  for (directory in draw_dirs) {
    lines <- readLines(file.path(directory, "sha256.txt"), warn = FALSE)
    files <- substring(lines, 67L); hashes <- substr(lines, 1L, 64L)
    if (length(lines) != 2L || !setequal(files, c("result.csv", "draw.rds")) ||
        any(vapply(file.path(directory, files), roce_sha256_file,
                   character(1L)) != hashes))
      stop("draw checksum validation failed: ", directory, call. = FALSE)
  }
  saved <- lapply(file.path(draw_dirs, "draw.rds"), readRDS)
  csv <- lapply(file.path(draw_dirs, "result.csv"), function(path)
    read.csv(path, stringsAsFactors = FALSE, na.strings = c("", "NA")))
  pair <- .afri_validate_draw_pair(csv, saved)
  package_hash <- .weight_calibration_installed_package_fingerprint(
    package_lib, roce_sha256_file
  )
  if (any(vapply(csv, function(x) x$package_fingerprint != package_hash,
                 logical(1L))))
    stop("draw/package fingerprint mismatch.", call. = FALSE)
  fits <- list(identity = saved[[1L]]$result$fitted,
               draw1 = saved[[2L]]$result$fitted)
  lapply(fits, RoCE:::.validate_weight_bootstrap_result)
  source_bundle <- normalizePath(pair$source_bundle, mustWork = TRUE)
  source_checksum <- readLines(file.path(source_bundle, "sha256.txt"),
                               warn = FALSE)
  source_files <- substring(source_checksum, 67L)
  source_hashes <- substr(source_checksum, 1L, 64L)
  if (length(source_checksum) != 4L ||
      any(vapply(file.path(source_bundle, source_files), roce_sha256_file,
                 character(1L)) != source_hashes))
    stop("source calibration bundle checksum failed.", call. = FALSE)
  source_saved <- readRDS(file.path(source_bundle, "artifacts.rds"))
  source_fit <- .afri_select_source_fit(source_saved, pair$rho)
  RoCE:::.validate_weight_bootstrap_result(source_fit)
  identity_source_error <- max(c(
    abs(fits$identity$estimate - source_fit$estimate),
    abs(fits$identity$se - source_fit$se),
    .afri_max_error(fits$identity$fold_weights, source_fit$fold_weights),
    .afri_max_error(fits$identity$fold_wald_statistics,
                    source_fit$fold_wald_statistics),
    .afri_max_error(fits$identity$fold_penalty_coefficients,
                    source_fit$fold_penalty_coefficients),
    .afri_max_error(fits$identity$all_phi_agg, source_fit$all_phi_agg)
  ))
  rm(source_saved)
  gc(verbose = FALSE)

  fold_audit <- do.call(rbind, unlist(lapply(names(fits), function(label)
    lapply(seq_len(fits[[label]]$n_folds), function(k)
      .afri_manual_fold(fits[[label]], label, k))), recursive = FALSE))
  fold_source_audit <- do.call(rbind, unlist(lapply(names(fits), function(label)
    lapply(seq_len(fits[[label]]$n_folds), function(k)
      .afri_fold_source_diagnostics(fits[[label]], label, k))),
    recursive = FALSE))
  arm_audit <- do.call(rbind, lapply(names(fits), function(label)
    .afri_arm_audit(fits[[label]], label)))

  reference <- fits$identity
  draw1 <- fits$draw1
  counts <- saved[[2L]]$result$multiplicities
  sources <- names(reference$weights)
  reweighted_inner <- lapply(
    reference$intermediates$inner_fold_info,
    RoCE:::.reweight_tate_inner_info,
    target_multiplier = counts$t,
    source_multipliers = counts[sources],
    source_names = sources
  )
  matched_phase2 <- RoCE:::.compute_phase2_weights(
    n_folds = reference$n_folds, inner_fold_info = reweighted_inner,
    fold_info = reference$intermediates$fold_info,
    lambda_selection = reference$aggregation_lambda_selection,
    K = length(sources), verbose = FALSE,
    lambda_rule = reference$aggregation_lambda_rule,
    screening_rule = reference$aggregation_screening_rule
  )
  matched_saved <- saved[[2L]]$result$matched_fixed_nuisance
  matched_audit <- data.frame(
    metric = c("fold_weights", "estimate_original_weights",
               "estimate_relearned_weights"),
    maximum_error = c(
      .afri_max_error(matched_phase2$fold_weights,
                      matched_saved$fold_weights),
      abs(RoCE:::.weight_bootstrap_outer_estimate(
        reference$intermediates$fold_info, reference$fold_weights,
        counts$t, counts[sources], sources,
        reference$intermediates$sample_sizes$n_t,
        reference$intermediates$sample_sizes$n_source
      ) - matched_saved$estimate_original_weights),
      abs(RoCE:::.weight_bootstrap_outer_estimate(
        reference$intermediates$fold_info, matched_phase2$fold_weights,
        counts$t, counts[sources], sources,
        reference$intermediates$sample_sizes$n_t,
        reference$intermediates$sample_sizes$n_source
      ) - matched_saved$estimate_relearned_weights)
    ), stringsAsFactors = FALSE
  )
  input_fields <- c("V_ot", "V_t", "V_s", "C_ot", "C_cross",
                    "avg_target_est", "avg_source_est")
  input_comparison <- do.call(rbind, lapply(seq_len(reference$n_folds),
    function(k) do.call(rbind, lapply(input_fields, function(field) data.frame(
      fold = k, field = field,
      matched_minus_original_max = .afri_max_error(
        reweighted_inner[[k]][[field]],
        reference$intermediates$inner_fold_info[[k]][[field]]
      ),
      refit_minus_matched_max = .afri_max_error(
        draw1$intermediates$inner_fold_info[[k]][[field]],
        reweighted_inner[[k]][[field]]
      ), stringsAsFactors = FALSE
    )))))
  component_tails <- do.call(rbind, lapply(names(fits), function(label)
    .afri_component_tails(fits[[label]], label)))
  final_tails <- do.call(rbind, lapply(names(fits), function(label)
    .afri_final_tails(fits[[label]], label)))
  overall <- data.frame(
    metric = c("identity_refit_vs_source_fit", "estimate_mean_phi",
               "variance_site_stratified",
               "arm_contrast", "wald", "penalty", "optimizer_rerun",
               "optimizer_kkt", "matched_fixed_nuisance"),
    maximum_error = c(
      identity_source_error,
      max(vapply(fits, function(f) abs(f$estimate - mean(f$all_phi_agg)),
                 numeric(1L))),
      max(vapply(fits, function(f) abs(f$variance -
        RoCE:::.multisite_pseudovalue_variance(
          f$all_phi_agg, c(f$intermediates$sample_sizes$n_t,
                           f$intermediates$sample_sizes$n_source))), numeric(1L))),
      max(c(arm_audit$fold_arm_contrast_error,
            arm_audit$inner_arm_contrast_error)),
      max(fold_audit$wald_reconstruction_error),
      max(fold_audit$penalty_reconstruction_error),
      max(fold_audit$optimizer_rerun_weight_error),
      max(fold_audit$normalized_kkt_residual_max),
      max(matched_audit$maximum_error)
    ), stringsAsFactors = FALSE
  )
  if (any(overall$maximum_error[
        overall$metric != "optimizer_kkt"
      ] > 1e-10) ||
      overall$maximum_error[overall$metric == "optimizer_kkt"] > 1e-5)
    stop("aggregation intermediate audit failed tolerance.", call. = FALSE)
  script_hash <- roce_sha256_file(
    "diagnosis/tate_common_weight/audit_full_refit_aggregation_intermediates.R"
  )
  roce_write_atomic_directory(output, function(stage) {
    tables <- list(overall_audit.csv = overall, fold_optimizer_audit.csv = fold_audit,
      fold_source_diagnostics.csv = fold_source_audit,
      arm_contrast_audit.csv = arm_audit, matched_fixed_audit.csv = matched_audit,
      optimizer_input_comparison.csv = input_comparison,
      influence_component_tails.csv = component_tails,
      final_pseudovalue_tails.csv = final_tails)
    for (name in names(tables))
      write.csv(tables[[name]], file.path(stage, name), row.names = FALSE)
    writeLines(c("aggregation_intermediate_audit=passed",
      paste0("package_fingerprint=", package_hash),
      paste0("audit_script_fingerprint=", script_hash),
      "weight_domain=unconstrained_signed_RK; no simplex/nonnegative projection",
      "interpretation=implementation and intermediate-variable audit only; no coverage claim"),
      file.path(stage, "metadata.txt"))
    files <- c(names(tables), "metadata.txt")
    hashes <- vapply(file.path(stage, files), roce_sha256_file, character(1L))
    writeLines(paste(hashes, files, sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "full-refit aggregation intermediate audit")
  print(overall)
  invisible(overall)
}

if (sys.nframe() == 0L) audit_full_refit_aggregation_main()
