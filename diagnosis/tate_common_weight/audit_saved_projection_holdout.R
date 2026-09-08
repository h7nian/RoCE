#!/usr/bin/env Rscript

# Evaluate the already fitted probe coefficients on untouched outer fold 1.
# All three probe penalties are retained; no truth/coverage-based selection.
.projection_holdout_scores <- function(state, block, site, coefficients) {
  stopifnot(site %in% c("target", "source"), length(coefficients) == length(state$theta))
  theta <- state$theta; index <- state$indices
  W <- block$W; Z <- block$Z; A <- block$A; Y <- block$Y
  prediction <- plogis(drop(W %*% theta[index$alpha_final]))
  log_weight <- drop(Z %*% theta[index$gamma_final])
  clip <- function(x, bound) pmax(-bound, pmin(bound, x))
  original <- if (site == "target") prediction else
    A*exp(-clip(log_weight, state$M_inference))*(Y-prediction)
  terms <- matrix(0, nrow(W), length(index), dimnames = list(NULL, names(index)))
  direct_means <- setNames(numeric(length(index)), names(index))
  mean_initial_weight <- mean_initial_derivative <- numeric(nrow(W))
  source_total <- nrow(state$source$W)
  for (key in state$keys) {
    alpha_name <- paste0("alpha_initial_", key)
    gamma_name <- paste0("gamma_initial_", key)
    a0 <- theta[index[[alpha_name]]]; g0 <- theta[index[[gamma_name]]]
    initial_prediction <- plogis(drop(W %*% a0))
    if (site == "target") {
      alpha_multiplier <- A*(Y-initial_prediction)
      gamma_multiplier <- rep(1, nrow(W))
      mean_initial_derivative <- mean_initial_derivative+
        initial_prediction*(1-initial_prediction)/length(state$keys)
    } else {
      alpha_multiplier <- numeric(nrow(W))
      gamma_multiplier <- -A*exp(-drop(Z %*% g0))
      source_fraction <- nrow(state$pieces[[key]]$source_calibration$W)/source_total
      mean_initial_weight <- mean_initial_weight+
        source_fraction*exp(-clip(drop(Z %*% g0), state$M_fit))
      initial_clipped <- plogis(clip(drop(W %*% a0), state$M_fit))
      mean_initial_derivative <- mean_initial_derivative+
        source_fraction*initial_clipped*(1-initial_clipped)
    }
    terms[, alpha_name] <- alpha_multiplier*drop(W %*% coefficients[index[[alpha_name]]])
    terms[, gamma_name] <- gamma_multiplier*drop(Z %*% coefficients[index[[gamma_name]]])
    direct_means[alpha_name] <- sum(colMeans(W*alpha_multiplier)*coefficients[index[[alpha_name]]])
    direct_means[gamma_name] <- sum(colMeans(Z*gamma_multiplier)*coefficients[index[[gamma_name]]])
  }
  alpha_multiplier <- if (site == "target") numeric(nrow(W)) else
    A*mean_initial_weight*(Y-prediction)
  gamma_multiplier <- if (site == "target") mean_initial_derivative else
    -A*exp(-log_weight)*mean_initial_derivative
  terms[, "alpha_final"] <- alpha_multiplier*drop(W %*% coefficients[index$alpha_final])
  terms[, "gamma_final"] <- gamma_multiplier*drop(Z %*% coefficients[index$gamma_final])
  direct_means["alpha_final"] <- sum(colMeans(W*alpha_multiplier)*coefficients[index$alpha_final])
  direct_means["gamma_final"] <- sum(colMeans(Z*gamma_multiplier)*coefficients[index$gamma_final])
  identity_error <- max(abs(colMeans(terms)-direct_means))
  corrected <- original-rowSums(terms)
  stopifnot(identity_error < 1e-10, all(is.finite(corrected)))
  list(original = original, corrected = corrected, terms = terms, identity_error = identity_error)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) stop("usage: audit_saved_projection_holdout.R V19_LIBRARY MATRIX_PROBE SEED_BUNDLE OUTPUT")
  lib <- normalizePath(args[1], mustWork = TRUE)
  probe <- normalizePath(args[2], mustWork = TRUE)
  seed <- normalizePath(args[3], mustWork = TRUE)
  output <- args[4]
  if (file.exists(output)) stop("holdout output already exists")
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/run_weight_bootstrap_calibration.R")
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = lib))
  stopifnot(normalizePath(find.package("RoCE")) == file.path(lib, "RoCE"),
    .weight_calibration_installed_package_fingerprint(lib, roce_sha256_file) ==
      "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7")
  read_checked_rds <- function(directory, name, expected_manifest = NULL) {
    manifest <- file.path(directory, "sha256.txt")
    if (!is.null(expected_manifest)) stopifnot(roce_sha256_file(manifest) == expected_manifest)
    lines <- readLines(manifest)
    expected <- substr(lines[substring(lines, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(directory, name)) == expected)
    readRDS(file.path(directory, name))
  }
  matrix_probe <- read_checked_rds(probe, "matrix_probe.rds")
  bundle <- read_checked_rds(seed, "artifacts.rds",
    "f7f3bcdf901a395062fd8864905a4bdaa37dcf4d7be34cd3fae1165a3ee09575")
  state <- matrix_probe$state
  fixture <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", fixture)
  stopifnot(isTRUE(all.equal(state, fixture$.saved_projection_state(bundle), tolerance = 0)))
  artifact <- bundle$group_result$artifacts[["0"]]
  f <- artifact$direct_tate_results$one_round_crossfit
  info <- f$intermediates$fold_info[[1L]]
  heldout <- list(
    target = fixture$.saved_projection_block(artifact$data_split$t, info$target_idx, 1L),
    source = fixture$.saved_projection_block(artifact$data_split$s1, info$source_idx[[1L]], 1L))
  for (site in names(heldout)) {
    stopifnot(!any(heldout[[site]]$index %in% state[[site]]$index))
    for (piece in state$pieces) {
      for (part in c("training", "calibration")) {
        stopifnot(!any(heldout[[site]]$index %in% piece[[paste(site, part, sep = "_")]]$index))
      }
    }
  }
  source_fit <- f$arm_results$mu1$fold_results[[1L]]$source_results$s1
  results <- lapply(matrix_probe$attempts, function(attempt) {
    if (!is.null(attempt$failure)) stop("probe contains an unsolved coefficient attempt; retain and review it before holdout use")
    scores <- lapply(names(heldout), function(site)
      .projection_holdout_scores(state, heldout[[site]], site, attempt$coefficients))
    names(scores) <- names(heldout)
    original_mean <- sum(vapply(scores, function(x) mean(x$original), numeric(1)))
    corrected_mean <- sum(vapply(scores, function(x) mean(x$corrected), numeric(1)))
    stopifnot(abs(original_mean-source_fit$mu_ts) < 1e-12)
    if (attempt$summary$fraction == 1) {
      stopifnot(all(attempt$coefficients == 0), original_mean == corrected_mean)
    }
    rows <- do.call(rbind, lapply(names(scores), function(site) {
      x <- scores[[site]]; centered <- x$corrected-mean(x$corrected)
      original_centered <- x$original-mean(x$original)
      ss <- sum(centered^2)
      data.frame(fraction = attempt$summary$fraction, site,
        n = length(centered), mean_original = mean(x$original), mean_corrected = mean(x$corrected),
        mean_adjustment = mean(x$corrected-x$original),
        original_score_sd = sd(x$original), corrected_score_sd = sd(x$corrected),
        max_abs_centered_original = max(abs(original_centered)),
        max_abs_centered_corrected = max(abs(centered)),
        top1_ss_fraction = if (ss == 0) NA_real_ else max(centered^2)/ss,
        block_mean_identity_error = x$identity_error)
    }))
    blocks <- do.call(rbind, lapply(names(scores), function(site)
      data.frame(fraction = attempt$summary$fraction, site,
        block = colnames(scores[[site]]$terms), mean_adjustment = -colMeans(scores[[site]]$terms))))
    list(summary = data.frame(fraction = attempt$summary$fraction,
      penalty = attempt$summary$penalty, coefficient_l1 = sum(abs(attempt$coefficients)),
      original_fold_mean = original_mean, adjusted_fold_mean = corrected_mean),
      sites = rows, blocks = blocks, scores = scores)
  })
  summaries <- do.call(rbind, lapply(results, `[[`, "summary"))
  sites <- do.call(rbind, lapply(results, `[[`, "sites"))
  blocks <- do.call(rbind, lapply(results, `[[`, "blocks"))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(summaries, file.path(stage, "fold_summary.csv"), row.names = FALSE)
    write.csv(sites, file.path(stage, "site_score_stability.csv"), row.names = FALSE)
    write.csv(blocks, file.path(stage, "block_adjustments.csv"), row.names = FALSE)
    saveRDS(results, file.path(stage, "heldout_scores.rds"))
    writeLines(c("saved_projection_holdout_check=complete", "scope=single_source_single_arm_outer_fold",
      "sim_id=10013", "rho=0", "outer_fold=1", "source=s1", "arm=1",
      "holdout_ids_disjoint_from_all_fitting_and_projection_rows=TRUE",
      "all_three_predeclared_probe_fractions_retained=TRUE", "penalty_selected=FALSE",
      "truth_or_coverage_used_for_selection=FALSE", "new_mc_replications=0", "nuisance_refits=0",
      "primary_estimator_se_ci_changed=FALSE", "inference_validated=FALSE",
      paste0("probe_manifest_sha256=", roce_sha256_file(file.path(probe, "sha256.txt"))),
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/audit_saved_projection_holdout.R"))),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "saved nuisance projection holdout check")
  print(summaries, row.names = FALSE); print(sites, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
