#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L) stop("usage: audit_source_validation_split.R LIBRARY SEED FIT_BUNDLE OUTPUT")
  .libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
  suppressPackageStartupMessages(library(RoCE, lib.loc = args[1]))
  source("scripts/slurm/result_provenance.R")
  source("scripts/slurm/atomic_output.R")
  source("diagnosis/tate_common_weight/source_working_basis.R")
  equations <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/probe_saved_nuisance_projection.R", equations)
  holdout <- new.env(parent = globalenv())
  sys.source("diagnosis/tate_common_weight/audit_saved_projection_holdout.R", holdout)
  checked <- function(root, name) {
    manifest <- readLines(file.path(root, "sha256.txt"))
    expected <- substr(manifest[substring(manifest, 67L) == name], 1L, 64L)
    stopifnot(length(expected) == 1L, roce_sha256_file(file.path(root, name)) == expected)
    readRDS(file.path(root, name))
  }
  artifact <- checked(args[2], "artifacts.rds")$group_result$artifacts[["0"]]
  captured <- checked(args[3], "source_validation_states.rds")
  working_config <- if (is.null(captured$working_config)) "C1" else captured$working_config
  data <- .source_working_basis(artifact$data_split, working_config)
  info <- artifact$direct_tate_results$one_round_crossfit$intermediates$fold_info
  ids <- list(target = lapply(info, `[[`, "target_idx"), source = lapply(info, function(x) x$source_idx[[1L]]))
  fold_map <- captured$global_fold_map
  stopifnot(length(fold_map) %in% c(4L,5L),
    if(length(fold_map)==4L) setequal(fold_map,2:5) else identical(fold_map,1:5),
    length(captured$splits) == length(fold_map)-1L)
  validation_fold <- fold_map[1L]
  training_folds <- fold_map[-1L]
  for (site in names(ids)) {
    stopifnot(identical(captured$validation_ids[[site]], ids[[site]][[validation_fold]]),
      identical(captured$excluded_outer_ids[[site]], ids[[site]][[1L]]))
    for (split in captured$splits) {
      k <- split$global_calibration_fold
      stopifnot(k %in% training_folds, fold_map[split$local_calibration_fold] == k,
        identical(split[[paste0(site, "_calibration_ids")]], ids[[site]][[k]]),
        identical(split[[paste0(site, "_training_ids")]],
          unlist(ids[[site]][setdiff(training_folds, k)], use.names = FALSE)))
    }
  }
  states <- reports <- list()
  for (arm in 0:1) {
    result <- captured$results[[paste0("mu", arm)]]
    stopifnot(!result$failed, is.null(result$failure))
    fit <- result$result
    keys <- paste0("k2_", seq.int(2L,length(fold_map)))
    stopifnot(identical(names(fit$per_k2_alpha), keys), identical(names(fit$per_k2_gamma), keys))
    pieces <- setNames(lapply(captured$splits, function(split) {
      answer <- list()
      for (site in names(ids)) for (part in c("training", "calibration")) {
        name <- paste(site, part, sep = "_")
        answer[[name]] <- equations$.saved_projection_block(data[[if (site == "target") "t" else "s1"]],
          split[[paste0(name, "_ids")]], arm)
      }
      answer
    }), keys)
    values <- c(setNames(lapply(fit$per_k2_alpha, as.numeric), paste0("alpha_initial_", keys)),
      setNames(lapply(fit$per_k2_gamma, as.numeric), paste0("gamma_initial_", keys)),
      list(alpha_final = as.numeric(fit$alpha_ts), gamma_final = as.numeric(fit$gamma_s)))
    sizes <- lengths(values); ends <- cumsum(sizes)
    indices <- setNames(lapply(seq_along(sizes), function(k) seq.int(ends[k]-sizes[k]+1L, ends[k])), names(values))
    state <- list(pieces = pieces, keys = keys, indices = indices, theta = unlist(values, use.names = FALSE),
      target = equations$.saved_projection_block(data$t, unlist(ids$target[training_folds], use.names = FALSE), arm),
      source = equations$.saved_projection_block(data$s1, unlist(ids$source[training_folds], use.names = FALSE), arm),
      M_fit = 5, M_inference = 5, logistic_clip = RoCE:::LOGISTIC_CLIP)
    evaluated <- equations$.evaluate_saved_projection(state, state$theta, TRUE)
    errors <- vapply(1:6, function(k) {
      direction <- sin(seq_along(state$theta)*k); direction <- direction/sqrt(sum(direction^2))
      epsilon <- 1e-6
      plus <- equations$.evaluate_saved_projection(state, state$theta+epsilon*direction)
      minus <- equations$.evaluate_saved_projection(state, state$theta-epsilon*direction)
      c(jacobian = max(abs((plus$moments-minus$moments)/(2*epsilon)-drop(evaluated$jacobian %*% direction))),
        score = abs((plus$score-minus$score)/(2*epsilon)-sum(evaluated$score_gradient*direction)))
    }, numeric(2L))
    kkt <- vapply(c("alpha_final", "gamma_final"), function(name) {
      index <- indices[[name]]
      gradient <- evaluated$moments[index]*if (name == "alpha_final") -1 else 1
      parameter <- if (name == "alpha_final") fit$alpha_ts else fit$gamma_s
      lambda <- attr(parameter, "lambda_used")
      stopifnot(length(lambda) == 1L, is.finite(lambda))
      penalty <- c(0, rep(lambda, length(index)-1L))
      active <- parameter != 0
      residual <- pmax(abs(gradient)-penalty, 0)
      residual[active] <- abs(gradient[active]+penalty[active]*sign(parameter[active]))
      max(residual)
    }, numeric(1L))
    heldout_mean <- sum(vapply(names(ids), function(site) {
      block <- equations$.saved_projection_block(data[[if (site == "target") "t" else "s1"]], ids[[site]][[validation_fold]], arm)
      mean(holdout$.projection_holdout_scores(state, block, site, numeric(length(state$theta)))$original)
    }, numeric(1L)))
    stopifnot(max(errors) < 1e-7, max(kkt) < 1e-4, abs(heldout_mean-fit$mu_ts) < 1e-12)
    states[[paste0("mu", arm)]] <- list(state = state, evaluated = evaluated, derivative_errors = errors, kkt = kkt)
    reports[[paste0("mu", arm)]] <- data.frame(arm, dimension = length(state$theta),
      jacobian_error = max(errors[1, ]), score_gradient_error = max(errors[2, ]),
      maximum_final_kkt_error = max(kkt), heldout_point_error = abs(heldout_mean-fit$mu_ts),
      warning_count = length(result$warnings))
  }
  reports <- do.call(rbind, reports)
  roce_write_atomic_directory(args[4], function(stage) {
    saveRDS(states, file.path(stage, "validation_equation_states.rds"))
    write.csv(reports, file.path(stage, "equation_checks.csv"), row.names = FALSE)
    writeLines(c("inference_validated=FALSE", "projection_policy_selected=FALSE",
      paste0("script_sha256=", roce_sha256_file("diagnosis/tate_common_weight/audit_source_validation_split.R")),
      paste0("fit_manifest_sha256=", roce_sha256_file(file.path(args[3], "sha256.txt")))), file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "), file.path(stage, "sha256.txt"))
  }, caller = "source validation equation audit")
  print(reports, row.names = FALSE)
}
if (sys.nframe() == 0L) main()
