#!/usr/bin/env Rscript

# Independent read-only recalculation of the completed, prespecified n=100
# checkpoint. Does not refit, select seeds, or call the original summarizer.
main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) stop("usage: review_independent_n100.R PILOT_ROOT OUTPUT")
  root <- normalizePath(args[[1L]], mustWork = TRUE)
  output <- args[[2L]]
  if (file.exists(output)) stop("review output already exists")
  source("scripts/slurm/atomic_output.R")
  sha <- function(path) {
    value <- system2("sha256sum", shQuote(path), stdout = TRUE)
    if (!is.null(attr(value, "status")) || length(value) != 1L ||
        !grepl("^[0-9a-f]{64}  ", value)) stop("sha256sum failed: ", path)
    substr(value, 1L, 64L)
  }
  verify <- function(directory, expected) {
    manifest <- file.path(directory, "sha256.txt")
    lines <- readLines(manifest, warn = FALSE)
    names <- substring(lines, 67L)
    stopifnot(length(lines) == length(expected), !anyDuplicated(names),
              setequal(names, expected),
              setequal(list.files(directory, all.files = TRUE, no.. = TRUE),
                       c(expected, "sha256.txt")),
              all(grepl("^[0-9a-f]{64}  ", lines)),
              identical(unname(vapply(file.path(directory, names), sha, "")),
                        substr(lines, 1L, 64L)))
    sha(manifest)
  }
  metadata <- function(path) {
    lines <- readLines(path, warn = FALSE)
    keys <- sub("=.*$", "", lines)
    stopifnot(!anyDuplicated(keys), all(grepl("=", lines, fixed = TRUE)))
    setNames(sub("^[^=]*=", "", lines), keys)
  }
  summary_dir <- file.path(root, "summaries/n100")
  summary_hash <- verify(summary_dir, c(
    "rho_summary.csv", "all_raw_rows.csv", "all_inference_rows.csv",
    "input_references.csv", "attempt_status.csv", "metadata.txt"
  ))
  references <- read.csv(file.path(summary_dir, "input_references.csv"))
  attempts <- read.csv(file.path(summary_dir, "attempt_status.csv"))
  stopifnot(identical(references$task_id, 1:100),
            identical(references$sim_id, 10001:10100),
            identical(attempts$task_id, 1:100), all(attempts$status == "complete"))
  raw <- inference <- checks <- qc <- vector("list", 100L)
  for (i in 1:100) {
    seed <- sprintf("seed_%06d", 10000L + i)
    directory <- file.path(root, seed)
    audit <- file.path(root, "audits", seed)
    stopifnot(references$bundle_directory[i] == directory,
              references$seed_audit_directory[i] == audit)
    raw_hash <- verify(directory, c("results.csv", "inference_audit.csv",
                                    "diagnostic_qc.csv", "artifacts.rds", "metadata.txt"))
    audit_hash <- verify(audit, c("audit_passed.txt", "summary.csv",
                                  "numerical_checks.csv", "nuisance_diagnostics.csv"))
    stopifnot(raw_hash == references$bundle_checksum_manifest_hash[i],
              audit_hash == references$seed_audit_checksum_manifest_hash[i])
    gate <- metadata(file.path(audit, "audit_passed.txt"))
    stopifnot(gate[["independent_inference_seed_audit"]] == "passed",
              gate[["bundle_checksum_manifest_fingerprint"]] == raw_hash,
              gate[["package_fingerprint"]] ==
                "2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7",
              gate[["workflow_fingerprint"]] ==
                "9214165abf36fe03fbd72d271158a3cc7b9f07f153d067ba47ba88712a6dada4",
              gate[["manifest_fingerprint"]] ==
                "1dc5b92cf549b417a71eadf1cfe4be50bf8a8219e147e171d53f073ff763a16d")
    raw[[i]] <- read.csv(file.path(directory, "results.csv"))
    inference[[i]] <- read.csv(file.path(directory, "inference_audit.csv"))
    checks[[i]] <- read.csv(file.path(audit, "numerical_checks.csv"))
    qc[[i]] <- read.csv(file.path(directory, "diagnostic_qc.csv"))
    stopifnot(nrow(raw[[i]]) == 84L, nrow(inference[[i]]) == 6L,
              nrow(checks[[i]]) == 72L,
              all(raw[[i]]$sim_id == 10000L + i),
              all(inference[[i]]$sim_id == 10000L + i),
              !any(qc[[i]]$implementation_failed),
              all(checks[[i]]$maximum_error >= 0),
              all(is.finite(checks[[i]]$maximum_error)))
    if (i %% 10L == 0L) cat("verified seed bundles:", i, "\n")
  }
  raw <- do.call(rbind, raw)
  inference <- do.call(rbind, inference)
  checks <- do.call(rbind, checks)
  qc <- do.call(rbind, qc)
  canonical <- function(x) {
    x <- x[order(x$sim_id, x$rho, x$method), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  stopifnot(isTRUE(all.equal(canonical(raw), canonical(read.csv(
    file.path(summary_dir, "all_raw_rows.csv"))), tolerance = 1e-12)),
    isTRUE(all.equal(canonical(inference), canonical(read.csv(
      file.path(summary_dir, "all_inference_rows.csv"))), tolerance = 1e-12)))
  kkt <- grepl("kkt", checks$check, ignore.case = TRUE)
  stopifnot(any(kkt), all(checks$maximum_error[kkt] <= 1e-5),
            all(checks$maximum_error[!kkt] <= 1e-10))
  official <- read.csv(file.path(summary_dir, "rho_summary.csv"))
  stopifnot(identical(official$rho, c(0, .5, 1, 1.5, 2, 2.5)))
  mcse <- function(x) sd(x) / sqrt(length(x))
  exact <- function(x) {
    k <- sum(x); n <- length(x)
    c(if (k == 0) 0 else qbeta(.025, k, n-k+1),
      if (k == n) 1 else qbeta(.975, k+1, n-k))
  }
  comparisons <- target_summary <- vector("list", 6L)
  for (j in seq_len(6L)) {
    rho <- official$rho[j]
    soft <- canonical(raw[raw$rho == rho & raw$method == "one_round_crossfit_ate", ])
    target <- canonical(raw[raw$rho == rho & raw$method == "target_only_ate", ])
    inf <- canonical(inference[inference$rho == rho, ])
    stopifnot(identical(soft$sim_id, 10001:10100),
              identical(target$sim_id, soft$sim_id), identical(inf$sim_id, soft$sim_id))
    e <- soft$estimate - soft$truth
    t <- target$estimate - target$truth
    coverage <- soft$ci_lower <= soft$truth & soft$ci_upper >= soft$truth
    target_coverage <- target$ci_lower <= target$truth & target$ci_upper >= target$truth
    bse <- soft$se_weight_relearn_bootstrap
    bcov <- abs(e) <= qnorm(.975) * bse
    stopifnot(all(coverage == soft$coverage), all(bcov == inf$weight_relearn_coverage),
              max(abs(inf$estimate - soft$estimate)) < 1e-12,
              max(abs(inf$analytic_se - soft$se)) < 1e-12,
              max(abs(inf$weight_relearn_bootstrap_se - bse)) < 1e-12)
    ci <- exact(coverage); bci <- exact(bcov)
    recomputed <- c(n_independent_seeds = 100, mean_error = mean(e),
      mean_error_mcse = mcse(e), empirical_sd = sd(e), rmse = sqrt(mean(e^2)),
      target_rmse = sqrt(mean(t^2)),
      paired_mse_difference_soft_minus_target = mean(e^2-t^2),
      paired_mse_difference_mcse = mcse(e^2-t^2), mean_analytic_se = mean(soft$se),
      mean_fixed_weight_bootstrap_se = mean(soft$se_fixed_weight_bootstrap),
      mean_weight_relearn_bootstrap_se = mean(bse),
      analytic_se_to_empirical_sd = mean(soft$se)/sd(e),
      fixed_weight_se_to_empirical_sd = mean(soft$se_fixed_weight_bootstrap)/sd(e),
      weight_relearn_se_to_empirical_sd = mean(bse)/sd(e),
      mean_analytic_ci_length = mean(soft$ci_upper-soft$ci_lower),
      mean_weight_relearn_ci_length = mean(2*qnorm(.975)*bse),
      analytic_coverage = mean(coverage), analytic_coverage_exact_lower = ci[1],
      analytic_coverage_exact_upper = ci[2], weight_relearn_coverage = mean(bcov),
      weight_relearn_coverage_exact_lower = bci[1], weight_relearn_coverage_exact_upper = bci[2],
      analytic_coverage_minus_target = mean(coverage-target_coverage),
      analytic_coverage_minus_target_mcse = mcse(coverage-target_coverage),
      weight_relearn_coverage_minus_target = mean(bcov-target_coverage),
      weight_relearn_coverage_minus_target_mcse = mcse(bcov-target_coverage),
      weight_relearn_coverage_minus_analytic = mean(bcov-coverage),
      weight_relearn_coverage_minus_analytic_mcse = mcse(bcov-coverage))
    errors <- abs(recomputed - unlist(official[j, names(recomputed)], use.names = FALSE))
    stopifnot(all(is.finite(errors)), max(errors) < 1e-10)
    comparisons[[j]] <- data.frame(rho, field = names(recomputed), error = errors)
    target_summary[[j]] <- data.frame(rho, bias = mean(t), bias_mcse = mcse(t),
      empirical_sd = sd(t), rmse = sqrt(mean(t^2)), mean_se = mean(target$se),
      coverage = mean(target_coverage))
  }
  comparisons <- do.call(rbind, comparisons)
  target_summary <- do.call(rbind, target_summary)
  stopifnot(!any(official$inference_validated), all(official$statistical_review_required),
            !any(official$warning_capture_complete), all(is.na(official$warning_count)))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(comparisons, file.path(stage, "statistical_recalculation.csv"), row.names = FALSE)
    write.csv(target_summary, file.path(stage, "target_summary.csv"), row.names = FALSE)
    writeLines(c("independent_n100_integrity_and_statistical_review=passed",
      paste0("summary_checksum_manifest_fingerprint=", summary_hash),
      paste0("review_script_sha256=", sha("diagnosis/tate_common_weight/review_independent_n100.R")),
      "raw_bundles_verified=100", "audit_bundles_verified=100",
      "raw_rows=8400", "inference_rows=600", "numerical_checks=7200",
      paste0("max_kkt_error=", max(checks$maximum_error[kkt])),
      paste0("max_non_kkt_error=", max(checks$maximum_error[!kkt])),
      paste0("max_statistical_recalculation_error=", max(comparisons$error)),
      paste0("hard_nuisance_failures=", sum(qc$nuisance_nonconverged_fits_total, na.rm = TRUE)),
      paste0("weight_bootstrap_failures=", sum(qc$weight_bootstrap_failures_total, na.rm = TRUE)),
      "primary_influence_tail_review=not_performed_by_this_script",
      "inference_validated=FALSE", "warning_capture_complete=FALSE"),
      file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, sha, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "independent n100 checkpoint review")
  print(official[, c("rho", "rmse", "target_rmse", "analytic_coverage",
                       "weight_relearn_coverage")], row.names = FALSE)
  print(target_summary, row.names = FALSE)
}

if (sys.nframe() == 0L) main()
