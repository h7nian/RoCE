#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2L) {
    stop(
      "usage: select_grouped_cutoff.R SUMMARY_DIRECTORY OUTPUT_DIRECTORY",
      call. = FALSE
    )
  }
  summary_directory <- normalizePath(args[[1L]], mustWork = TRUE)
  output_directory <- args[[2L]]
  source(file.path("scripts", "slurm", "atomic_output.R"))
  source(file.path("scripts", "slurm", "result_provenance.R"))

  audit_path <- file.path(summary_directory, "audit_passed.txt")
  summary_path <- file.path(summary_directory, "cutoff_summary.csv")
  if (!file.exists(audit_path) || !file.exists(summary_path)) {
    stop("cutoff selection requires an audited grouped summary.", call. = FALSE)
  }
  audit_lines <- readLines(audit_path, warn = FALSE)
  require_audit_value <- function(key, expected = NULL) {
    prefix <- paste0(key, "=")
    values <- sub(prefix, "", audit_lines[startsWith(audit_lines, prefix)])
    if (length(values) != 1L || is.na(values) || !nzchar(values)) {
      stop("cutoff audit has no unique value for ", key, ".", call. = FALSE)
    }
    if (!is.null(expected) && !identical(values, expected)) {
      stop(
        "cutoff audit mismatch for ", key, ": expected ", expected,
        ", observed ", values, ".", call. = FALSE
      )
    }
    values
  }
  require_audit_value("grouped_cutoff_diagnostic", "passed")
  require_audit_value("expected_replications", "10")
  require_audit_value("simulation_ids", "501--510")
  require_audit_value(
    "production_primary_cutoff_identity", "not_run_disjoint_pilot"
  )
  pilot_primary_cutoff <- suppressWarnings(as.numeric(
    require_audit_value("primary_cutoff")
  ))
  if (length(pilot_primary_cutoff) != 1L ||
      !is.finite(pilot_primary_cutoff) || pilot_primary_cutoff <= 0) {
    stop("cutoff audit has an invalid tested-package primary cutoff.",
         call. = FALSE)
  }
  package_fingerprint <- require_audit_value("package_fingerprint")
  workflow_fingerprint <- require_audit_value("workflow_fingerprint")
  manifest_fingerprint <- require_audit_value("manifest_fingerprint")
  summary_workflow_fingerprint <- require_audit_value(
    "summary_workflow_fingerprint"
  )
  summary_fingerprint <- require_audit_value("cutoff_summary_sha256")
  fingerprint_values <- c(
    package_fingerprint, workflow_fingerprint, manifest_fingerprint,
    summary_workflow_fingerprint, summary_fingerprint
  )
  if (any(!grepl("^[0-9a-f]{64}$", tolower(fingerprint_values)))) {
    stop("cutoff audit contains a malformed SHA-256 fingerprint.",
         call. = FALSE)
  }
  observed_summary_fingerprint <- roce_sha256_file(summary_path)
  if (!identical(
      tolower(summary_fingerprint), observed_summary_fingerprint
  )) {
    stop("cutoff summary fingerprint differs from its audit gate.",
         call. = FALSE)
  }
  selection_workflow_fingerprint <- tolower(trimws(Sys.getenv(
    "ROCE_CUTOFF_SELECTION_WORKFLOW_FINGERPRINT", ""
  )))
  if (!grepl("^[0-9a-f]{64}$", selection_workflow_fingerprint)) {
    stop(
      "ROCE_CUTOFF_SELECTION_WORKFLOW_FINGERPRINT must be SHA-256.",
      call. = FALSE
    )
  }

  summary <- read.csv(summary_path, stringsAsFactors = FALSE)
  required_columns <- c(
    "config", "p", "K", "rho", "cutoff", "n_replications", "rmse",
    "coverage", "coverage_mcse", "target_only_rmse"
  )
  missing_columns <- setdiff(required_columns, names(summary))
  if (length(missing_columns) > 0L) {
    stop(
      "cutoff summary is missing: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  expected_rhos <- c(0, 1, 2.5)
  expected_cutoffs <- c(1, 1.5, 2, 2.5, 3)
  if (anyNA(summary[required_columns]) ||
      nrow(summary) != length(expected_rhos) * length(expected_cutoffs) ||
      !identical(sort(unique(summary$rho)), expected_rhos) ||
      !identical(sort(unique(summary$cutoff)), expected_cutoffs) ||
      any(summary$config != "C1") || any(summary$p != 100L) ||
      any(summary$K != 2L) || any(summary$n_replications != 10L) ||
      any(!is.finite(unlist(summary[c(
        "rho", "cutoff", "rmse", "coverage", "coverage_mcse",
        "target_only_rmse"
      )]))) ||
      any(summary$rmse < 0) || any(summary$target_only_rmse < 0) ||
      any(summary$coverage < 0 | summary$coverage > 1) ||
      any(abs(summary$coverage * 10 - round(summary$coverage * 10)) > 1e-12) ||
      any(summary$coverage_mcse < 0) ||
      any(abs(
        summary$coverage_mcse -
          sqrt(summary$coverage * (1 - summary$coverage) / 10)
      ) > 1e-12) ||
      anyDuplicated(summary[c("rho", "cutoff")])) {
    stop("cutoff summary does not match the audited disjoint pilot grid.",
         call. = FALSE)
  }
  target_rmse_by_rho <- split(summary$target_only_rmse, summary$rho)
  if (!all(vapply(target_rmse_by_rho, function(values) {
    length(unique(values)) == 1L
  }, logical(1L)))) {
    stop("target-only RMSE changed across the cutoff grid.", call. = FALSE)
  }

  rho_zero <- summary[summary$rho == 0, , drop = FALSE]
  nontransport <- summary[summary$rho %in% c(1, 2.5), , drop = FALSE]
  rho_zero_min_rmse <- min(rho_zero$rmse)
  rho_zero_tolerance <- 1.05 * rho_zero_min_rmse
  candidates <- do.call(rbind, lapply(expected_cutoffs, function(cutoff) {
    zero <- rho_zero[rho_zero$cutoff == cutoff, , drop = FALSE]
    shifted <- nontransport[nontransport$cutoff == cutoff, , drop = FALSE]
    data.frame(
      cutoff = cutoff,
      aggregation_lambda = 1 / cutoff,
      rho0_rmse = zero$rmse,
      rho0_relative_to_best = zero$rmse / rho_zero_min_rmse,
      rho0_efficiency_eligible = zero$rmse <= rho_zero_tolerance,
      worst_nontransport_rmse = max(shifted$rmse),
      mean_nontransport_rmse = mean(shifted$rmse),
      min_nontransport_coverage_descriptive = min(shifted$coverage),
      stringsAsFactors = FALSE
    )
  }))
  eligible <- candidates[candidates$rho0_efficiency_eligible, , drop = FALSE]
  if (nrow(eligible) == 0L) {
    stop("no cutoff passes the rho-zero efficiency tolerance.", call. = FALSE)
  }
  eligible <- eligible[order(
    eligible$worst_nontransport_rmse,
    eligible$mean_nontransport_rmse,
    eligible$cutoff
  ), , drop = FALSE]
  selected_cutoff <- eligible$cutoff[[1L]]
  candidates$selected <- candidates$cutoff == selected_cutoff
  candidates$selection_rule <- paste0(
    "rho0_rmse<=1.05*best;minimize_worst_rmse_at_rho_1_and_2.5;",
    "tie_break_mean_rmse_then_smaller_cutoff;coverage_not_used"
  )
  candidates$pilot_simulation_ids <- "501--510"
  candidates$pilot_primary_cutoff <- pilot_primary_cutoff
  candidates$package_fingerprint <- package_fingerprint
  candidates$simulation_workflow_fingerprint <- workflow_fingerprint
  candidates$selection_workflow_fingerprint <-
    selection_workflow_fingerprint
  candidates$manifest_fingerprint <- manifest_fingerprint
  candidates$summary_workflow_fingerprint <- summary_workflow_fingerprint
  candidates$input_summary_fingerprint <- summary_fingerprint

  roce_write_atomic_directory(
    output_directory,
    writer = function(staging_directory) {
      decision_path <- file.path(staging_directory, "cutoff_decision.csv")
      write.csv(
        candidates,
        decision_path,
        row.names = FALSE
      )
      decision_fingerprint <- roce_sha256_file(decision_path)
      writeLines(c(
        "cutoff_selection=passed",
        paste0("selected_cutoff=", selected_cutoff),
        paste0("selected_aggregation_lambda=", 1 / selected_cutoff),
        "pilot_simulation_ids=501--510",
        paste0("pilot_primary_cutoff=", pilot_primary_cutoff),
        paste0("package_fingerprint=", package_fingerprint),
        paste0("simulation_workflow_fingerprint=", workflow_fingerprint),
        paste0(
          "selection_workflow_fingerprint=",
          selection_workflow_fingerprint
        ),
        paste0("manifest_fingerprint=", manifest_fingerprint),
        paste0(
          "summary_workflow_fingerprint=", summary_workflow_fingerprint
        ),
        paste0("input_summary_fingerprint=", summary_fingerprint),
        paste0("cutoff_decision_fingerprint=", decision_fingerprint),
        paste0("selection_rule=", candidates$selection_rule[[1L]])
      ), file.path(staging_directory, "cutoff_selection_passed.txt"))
    },
    caller = "disjoint-pilot cutoff selection"
  )
  print(candidates, row.names = FALSE)
  message("[pass] selected cutoff ", selected_cutoff)
}

main()
