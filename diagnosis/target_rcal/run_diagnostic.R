#!/usr/bin/env Rscript
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 4L || !(args[1L] %in% c("smoke", "target", "full")) ||
      !(args[2L] %in% c("C1", "C2", "C3"))) {
    stop("usage: run_diagnostic.R smoke|target|full C1|C2|C3 SEED OUTPUT_ROOT")
  }
  mode <- args[1L]
  configuration <- args[2L]
  seed <- as.integer(args[3L])
  rcal_backend <- match.arg(Sys.getenv("ROCE_RCAL_BACKEND", "native"),
                            c("native", "reference"))
  stopifnot(length(seed) == 1L, !is.na(seed), seed > 0L)
  library_path <- Sys.getenv("ROCE_PROJECT_LIB")
  stopifnot(nzchar(library_path))
  .libPaths(c(library_path, .libPaths()))
  suppressPackageStartupMessages(library(RoCE))
  stopifnot(normalizePath(dirname(getNamespaceInfo("RoCE", "path"))) ==
              normalizePath(library_path))
  source("diagnosis/target_rcal/target_rcal_helpers.R")
  source("diagnosis/target_rcal/native_rcal_helpers.R")
  source("scripts/slurm/result_provenance.R")
  output <- file.path(args[4L], paste(mode, rcal_backend, sep = "_"),
                      configuration, sprintf("seed_%04d", seed))
  if (dir.exists(output)) stop("Refusing to overwrite an existing diagnostic seed.")
  dir.create(output, recursive = TRUE)
  warnings <- character()
  on.exit(writeLines(warnings, file.path(output, "warnings.txt")), add = TRUE)
  code_files <- c("diagnosis/target_rcal/target_rcal_helpers.R",
                  "diagnosis/target_rcal/run_diagnostic.R",
                  "diagnosis/target_rcal/native_rcal_helpers.R")
  code_hashes <- vapply(code_files, roce_sha256_file, character(1L))
  started <- proc.time()[["elapsed"]]
  withCallingHandlers({
    set.seed(seed)
    data <- RoCE:::generate_simulation_data(
      n_total = 3000L, K = 2L, p = 100L, config = configuration,
      estimand_type = "superpopulation", outcome_type = "binary", dgp_type = "face",
      ate_deviation = 0, n_deviated_sites = 0L, deviation_mechanism = "treated_arm",
      warn_ignored = FALSE
    )
    data_split <- RoCE:::split_data_by_site(data)
    stopifnot(all(vapply(data_split, function(site) site$n, integer(1L)) == 1000L),
              ncol(data_split$t$W_outcome) == 200L)
    folds <- RoCE:::build_crossfit_folds(data_split, 10L)
    truth <- target_truth_predictions(data)
    reference <- if (mode == "full") RoCE::run_tate_crossfit(
      data_split, n_folds = 10L, communication_mode = "one_round", verbose = FALSE,
      nlambda_init = 100L, family = "binomial", precomputed_folds = folds,
      nuisance_lambda_rule = "min", screening_rule = "soft_penalty", n_cores = 2L
    ) else NULL
    variants <- list(mu1 = vector("list", 10L), mu0 = vector("list", 10L))
    diagnostic_rows <- list()
    for (arm_name in names(variants)) {
      A_val <- if (arm_name == "mu1") 1L else 0L
      for (k1 in if (mode == "smoke") 1L else seq_len(10L)) {
        baseline <- if (is.null(reference)) NULL else
          reference$arm_results[[arm_name]]$fold_results[[k1]]$target_only
        outer <- target_fold_variants(folds$target_folds, truth, k1, A_val,
                                      baseline = baseline, rcal_backend = rcal_backend)
        diagnostic_rows[[length(diagnostic_rows) + 1L]] <- attr(outer, "diagnostics")
        inner <- list()
        if (mode == "full") for (k2 in setdiff(seq_len(10L), k1)) {
          k2_name <- paste0("k2_", k2)
          baseline <- reference$arm_results[[arm_name]]$fold_results[[k1]]$
            target_only_inner[[k2_name]]
          inner[[k2_name]] <- target_fold_variants(
            folds$target_folds, truth, k1, A_val, k2 = k2, baseline = baseline,
            rcal_backend = rcal_backend
          )
          diagnostic_rows[[length(diagnostic_rows) + 1L]] <-
            attr(inner[[k2_name]], "diagnostics")
        }
        variants[[arm_name]][[k1]] <- list(outer = outer, inner = inner)
        cat(sprintf("completed %s %s seed=%d arm=%d fold=%d elapsed=%.1fs\n",
                    mode, configuration, seed, A_val, k1,
                    proc.time()[["elapsed"]] - started))
        flush.console()
      }
    }
    variant_names <- names(variants$mu1[[1L]]$outer)
    rows <- lapply(variant_names, function(variant) {
      arm_values <- lapply(variants, function(arm) {
        do.call(rbind, lapply(Filter(Negate(is.null), arm), function(fold) {
          fit <- fold$outer[[variant]]
          cbind(phi = fit$estimate + fit$varphi_ot, remainder = fit$expected_remainder)
        }))
      })
      contrast <- arm_values$mu1[, "phi"] - arm_values$mu0[, "phi"]
      estimate <- mean(contrast)
      se <- sqrt(mean((contrast - estimate)^2) / length(contrast))
      full <- if (mode == "full")
        reaggregate_target_variant(data_split, folds, reference, variants, variant) else NULL
      if (variant == "lasso" && mode == "full") {
        stopifnot(isTRUE(all.equal(full$estimate, reference$estimate, tolerance = 1e-12)),
                  isTRUE(all.equal(full$se, reference$se, tolerance = 1e-10)),
                  isTRUE(all.equal(full$fold_weights, reference$fold_weights, tolerance = 1e-10)))
      }
      truth_tate <- data$mu1_true - data$mu0_true
      data.frame(configuration, seed, mode, rcal_backend, variant,
                 intervention_scope = "target_anchor_only", n_site = 1000L, p = 100L,
                 n_folds = 10L, nlambda = 100L, target_estimate = estimate,
                 target_bias = estimate - truth_tate, target_se = se,
                 mu1_bias = mean(arm_values$mu1[, "phi"]) - data$mu1_true,
                 mu0_bias = mean(arm_values$mu0[, "phi"]) - data$mu0_true,
                 target_remainder = mean(arm_values$mu1[, "remainder"] -
                                           arm_values$mu0[, "remainder"]),
                 full_estimate = if (is.null(full)) NA_real_ else full$estimate,
                 full_bias = if (is.null(full)) NA_real_ else full$estimate - truth_tate,
                 full_se = if (is.null(full)) NA_real_ else full$se,
                 anchor_weight = if (is.null(full)) NA_real_ else 1 - sum(full$weights))
    })
    write.csv(do.call(rbind, rows), file.path(output, "results.csv"), row.names = FALSE)
    write.csv(do.call(rbind, diagnostic_rows), file.path(output, "fit_diagnostics.csv"),
              row.names = FALSE)
    if (mode == "full") saveRDS(list(reference = reference, variants = variants,
                                     folds = folds, data_split = data_split, truth = truth),
                                file.path(output, "fits.rds"))
    stopifnot(identical(code_hashes, vapply(code_files, roce_sha256_file, character(1L))))
    writeLines(c(paste0("library=", normalizePath(library_path)),
                 paste0("RCAL=", packageVersion("RCAL")),
                 paste0("helper_sha256=", code_hashes[[1L]]),
                 paste0("driver_sha256=", code_hashes[[2L]]),
                 paste0("native_helper_sha256=", code_hashes[[3L]]),
                 paste0("rcal_backend=", rcal_backend),
                 paste0("elapsed_seconds=", proc.time()[["elapsed"]] - started)),
               file.path(output, "metadata.txt"))
  }, warning = function(warning) {
    warnings <<- c(warnings, conditionMessage(warning))
    invokeRestart("muffleWarning")
  })
  writeLines(warnings, file.path(output, "warnings.txt"))
  writeLines("completed", file.path(output, "COMPLETED"))
}
if (sys.nframe() == 0L) main()
