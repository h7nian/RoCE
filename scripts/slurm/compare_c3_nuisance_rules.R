#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
source(file.path("scripts", "slurm", "atomic_output.R"))
if (length(args) != 3L) {
  stop(
    paste0(
      "usage: compare_c3_nuisance_rules.R MIN_RESULT_DIR ",
      "ONE_SE_RESULT_DIR OUTPUT_DIR"
    ),
    call. = FALSE
  )
}

result_directories <- vapply(
  args[1:2], normalizePath, character(1L), mustWork = TRUE
)
output_directory <- args[[3L]]
if (dir.exists(output_directory) || file.exists(output_directory)) {
  stop("comparison output already exists: ", output_directory, call. = FALSE)
}

methods <- c(
  "fitted_fitted", "oracle_outcome", "oracle_propensity", "oracle_both"
)
required_columns <- c(
  "sim_id", "method", "estimate", "se", "coverage", "truth",
  "experiment", "dgp_type", "outcome_family", "heterogeneity_type",
  "estimand_type", "estimand_scope", "config", "p", "K", "rho",
  "n_site", "n_folds", "nlambda_init", "nuisance_lambda_rule",
  "package_library", "package_fingerprint", "workflow_fingerprint"
)

# The diagnostic was renamed from diagnose_c3_target_remainder.R to
# diagnose_target_remainder.R; archived runs keep the old file and experiment
# labels, so both are accepted here.
REMAINDER_RAW_FILES <- c("target_remainder_raw.csv", "c3_target_remainder_raw.csv")
REMAINDER_EXPERIMENTS <- c(
  "target_remainder_diagnostic", "c3_target_remainder_diagnostic"
)

remainder_raw_path <- function(directory) {
  candidates <- file.path(directory, REMAINDER_RAW_FILES)
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0L) {
    stop("missing C3 raw result in: ", directory, call. = FALSE)
  }
  found[[1L]]
}

read_result <- function(directory, expected_rule) {
  path <- remainder_raw_path(directory)
  result <- utils::read.csv(path, stringsAsFactors = FALSE)
  missing_columns <- setdiff(required_columns, names(result))
  if (length(missing_columns) > 0L) {
    stop(
      "C3 result is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  expected_design <-
    result$experiment %in% REMAINDER_EXPERIMENTS &
    result$dgp_type == "face" & result$outcome_family == "binomial" &
    result$heterogeneity_type == "none" &
    result$estimand_type == "superpopulation" &
    result$estimand_scope == "tate" &
    result$config == "C3" & result$p == 100L & result$K == 4L &
    result$rho == 0 & result$n_site == 1000L & result$n_folds == 5L &
    result$nlambda_init == 100L &
    result$nuisance_lambda_rule == expected_rule
  if (!all(expected_design)) {
    stop(
      "C3 result does not match the expected p=100, K=4, rho=0 design/rule.",
      call. = FALSE
    )
  }
  if (any(!is.finite(as.matrix(result[c("estimate", "se", "truth")]))) ||
      any(result$se <= 0) || anyNA(result$coverage)) {
    stop("C3 result has invalid estimates, SEs, truths, or coverage.",
         call. = FALSE)
  }
  method_parts <- split(result$sim_id, result$method)
  if (!identical(sort(names(method_parts)), sort(methods)) ||
      !all(vapply(method_parts, function(ids) {
        identical(sort(as.integer(ids)), seq_len(500L))
      }, logical(1L)))) {
    stop("C3 result must contain methods x sim_id 1:500 exactly once.",
         call. = FALSE)
  }
  if (anyDuplicated(result[c("sim_id", "method")])) {
    stop("C3 result contains duplicate sim_id/method rows.", call. = FALSE)
  }
  libraries <- unique(as.character(result$package_library))
  fingerprints <- unique(tolower(as.character(result$package_fingerprint)))
  workflow_fingerprints <- unique(tolower(as.character(
    result$workflow_fingerprint
  )))
  if (length(libraries) != 1L || !nzchar(libraries) ||
      length(fingerprints) != 1L ||
      !grepl("^[0-9a-f]{64}$", fingerprints) ||
      length(workflow_fingerprints) != 1L ||
      !grepl("^[0-9a-f]{64}$", workflow_fingerprints)) {
    stop("C3 result has invalid or mixed package/workflow provenance.",
         call. = FALSE)
  }
  result
}

minimum <- read_result(result_directories[[1L]], "min")
one_se <- read_result(result_directories[[2L]], "1se")
if (!identical(
      unique(minimum$package_fingerprint),
      unique(one_se$package_fingerprint)
    ) || !identical(
      unique(minimum$workflow_fingerprint),
      unique(one_se$workflow_fingerprint)
    )) {
  stop("paired nuisance-rule diagnostics must share one code build/workflow.",
       call. = FALSE)
}
paired <- merge(
  minimum, one_se,
  by = c("sim_id", "method"),
  suffixes = c("_min", "_1se"),
  sort = TRUE
)
if (nrow(paired) != 2000L ||
    !identical(paired$truth_min, paired$truth_1se)) {
  stop("min and one-SE results are not an exact paired 500-replication design.",
       call. = FALSE)
}

oracle_both <- paired[paired$method == "oracle_both", , drop = FALSE]
oracle_alignment_fields <- c("estimate", "se", "coverage", "truth")
oracle_alignment_passed <- all(vapply(oracle_alignment_fields, function(field) {
  identical(
    oracle_both[[paste0(field, "_min")]],
    oracle_both[[paste0(field, "_1se")]]
  )
}, logical(1L)))
if (!oracle_alignment_passed) {
  stop(
    "double-oracle rows changed across rules; DGP/fold pairing is not exact.",
    call. = FALSE
  )
}

paired$estimate_difference_1se_minus_min <-
  paired$estimate_1se - paired$estimate_min
paired$coverage_gained <- !paired$coverage_min & paired$coverage_1se
paired$coverage_lost <- paired$coverage_min & !paired$coverage_1se

summary_parts <- lapply(methods, function(method_name) {
  rows <- paired[paired$method == method_name, , drop = FALSE]
  difference <- rows$estimate_difference_1se_minus_min
  data.frame(
    method = method_name,
    n_paired_replications = nrow(rows),
    mean_estimate_difference_1se_minus_min = mean(difference),
    difference_mcse = stats::sd(difference) / sqrt(nrow(rows)),
    bias_min = mean(rows$estimate_min - rows$truth_min),
    bias_1se = mean(rows$estimate_1se - rows$truth_1se),
    rmse_min = sqrt(mean((rows$estimate_min - rows$truth_min)^2)),
    rmse_1se = sqrt(mean((rows$estimate_1se - rows$truth_1se)^2)),
    coverage_min = mean(rows$coverage_min),
    coverage_1se = mean(rows$coverage_1se),
    coverage_gained = sum(rows$coverage_gained),
    coverage_lost = sum(rows$coverage_lost),
    paired_estimate_correlation = stats::cor(
      rows$estimate_min, rows$estimate_1se
    ),
    stringsAsFactors = FALSE
  )
})
paired_summary <- do.call(rbind, summary_parts)

interaction_remainder <- function(result) {
  estimates <- lapply(methods, function(method_name) {
    rows <- result[result$method == method_name, c("sim_id", "estimate")]
    names(rows)[[2L]] <- method_name
    rows
  })
  wide <- Reduce(function(left, right) {
    merge(left, right, by = "sim_id", sort = TRUE)
  }, estimates)
  wide$fitted_fitted - wide$oracle_outcome -
    wide$oracle_propensity + wide$oracle_both
}
remainder_min <- interaction_remainder(minimum)
remainder_1se <- interaction_remainder(one_se)
remainder_difference <- remainder_1se - remainder_min
interaction_summary <- data.frame(
  rule = c("min", "1se", "1se_minus_min"),
  mean_interaction_remainder = c(
    mean(remainder_min), mean(remainder_1se), mean(remainder_difference)
  ),
  mcse = c(
    stats::sd(remainder_min),
    stats::sd(remainder_1se),
    stats::sd(remainder_difference)
  ) / sqrt(500L),
  stringsAsFactors = FALSE
)

provenance <- data.frame(
  rule = c("min", "1se"),
  result_directory = result_directories,
  package_library = c(
    unique(minimum$package_library), unique(one_se$package_library)
  ),
  package_fingerprint = c(
    unique(minimum$package_fingerprint),
    unique(one_se$package_fingerprint)
  ),
  workflow_fingerprint = c(
    unique(minimum$workflow_fingerprint),
    unique(one_se$workflow_fingerprint)
  ),
  raw_result_md5 = unname(tools::md5sum(vapply(
    result_directories, remainder_raw_path, character(1L)
  ))),
  stringsAsFactors = FALSE
)

roce_write_atomic_directory(
  output_directory,
  writer = function(staging_directory) {
    utils::write.csv(
      paired, file.path(staging_directory, "paired_replications.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      paired_summary, file.path(staging_directory, "paired_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      interaction_summary,
      file.path(staging_directory, "interaction_remainder_summary.csv"),
      row.names = FALSE
    )
    utils::write.csv(
      provenance, file.path(staging_directory, "result_provenance.csv"),
      row.names = FALSE
    )
    writeLines(
      c(
        "paired_design=passed",
        "double_oracle_exact_alignment=passed",
        "n_paired_replications=500",
        "primary_simulation_nuisance_rule=min",
        "rule_choice_not_retuned_on_diagnostic_coverage=true",
        paste0(
          "package_fingerprint=",
          unique(minimum$package_fingerprint)
        ),
        paste0(
          "workflow_fingerprint=",
          unique(minimum$workflow_fingerprint)
        ),
        paste0(
          "fitted_coverage_min=",
          paired_summary$coverage_min[
            paired_summary$method == "fitted_fitted"
          ]
        ),
        paste0(
          "fitted_coverage_1se=",
          paired_summary$coverage_1se[
            paired_summary$method == "fitted_fitted"
          ]
        )
      ),
      file.path(staging_directory, "comparison_audit_passed.txt")
    )
  },
  caller = "paired C3 nuisance-rule comparison"
)

print(paired_summary, row.names = FALSE, digits = 7)
print(interaction_summary, row.names = FALSE, digits = 7)
message("[pass] paired C3 nuisance-rule comparison written to ",
        output_directory)
}

main()
