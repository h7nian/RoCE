#!/usr/bin/env Rscript

source("diagnosis/tate_common_weight/summarize_full_refit_draws.R")
provenance <- new.env(parent = baseenv())
sys.source("scripts/slurm/result_provenance.R", provenance)

make_root <- function() {
  root <- tempfile("full_refit_summary_fixture_")
  dir.create(root)
  root
}

write_draw <- function(root, directory_name, draw_id, status = "completed",
                       workflow = strrep("b", 64L), identity_shift = 0,
                       rds_shift = 0, reference_se = 0.03,
                       protocol = c("row_level", "origin_grouped"),
                       partition_records = NULL) {
  protocol <- match.arg(protocol)
  source_bundle <- file.path(root, "source_bundle")
  if (!dir.exists(source_bundle)) {
    dir.create(source_bundle)
    writeLines(paste0("package_fingerprint=", strrep("a", 64L)),
               file.path(source_bundle, "metadata.txt"))
    metadata_hash <- provenance$roce_sha256_file(
      file.path(source_bundle, "metadata.txt")
    )
    writeLines(paste(metadata_hash, "metadata.txt", sep = "  "),
               file.path(source_bundle, "sha256.txt"))
  }
  source_bundle_hash <- provenance$roce_sha256_file(
    file.path(source_bundle, "sha256.txt")
  )
  bundle <- file.path(root, directory_name)
  dir.create(bundle)
  identity <- draw_id == 0L
  base <- if (identity) 0.2 + identity_shift else 0.2 + draw_id / 100
  succeeded <- status == "completed"
  row <- data.frame(
    sim_id = 1L, rho = 0, draw_id = draw_id,
    resample_seed = 200000L + draw_id, identity = identity,
    status = status,
    failure_message = if (succeeded) NA_character_ else "synthetic failure",
    failure_stage = if (succeeded) NA_character_ else "nuisance_refit",
    estimate_refit_relearned_weights = if (succeeded) base else NA_real_,
    estimate_refit_original_weights = if (succeeded) {
      base - if (identity) identity_shift else 0.01
    } else NA_real_,
    estimate_fixed_nuisance_original_weights = if (succeeded) {
      base - if (identity) identity_shift else 0.02
    } else NA_real_,
    estimate_fixed_nuisance_relearned_weights = if (succeeded) {
      base - if (identity) identity_shift else 0.005
    } else NA_real_,
    estimate_target_refit = if (succeeded) 0.18 else NA_real_,
    estimate_reference = 0.2, se_analytic_reference = reference_se,
    elapsed_seconds = 1, nuisance_refit = TRUE,
    inference_validated = FALSE,
    conditional_on_saved_design_and_outer_partition = TRUE,
    nuisance_cv_groups_duplicate_origins = protocol == "origin_grouped",
    nuisance_cv_duplicate_origin_leakage_possible = protocol == "row_level",
    resampling_scheme = "site_by_original_outer_fold_multinomial",
    package_fingerprint = if (protocol == "origin_grouped") {
      strrep("c", 64L)
    } else strrep("a", 64L),
    workflow_fingerprint = workflow,
    source_bundle_fingerprint = source_bundle_hash,
    stringsAsFactors = FALSE
  )
  if (protocol == "origin_grouped") {
    row$nuisance_cv_partition <- "origin_grouped"
    row$source_package_fingerprint <- strrep("a", 64L)
  }
  write.csv(row, file.path(bundle, "result.csv"), row.names = FALSE, na = "")
  result <- if (succeeded) {
    list(
      estimate_refit_relearned_weights =
        row$estimate_refit_relearned_weights,
      estimate_refit_original_weights = row$estimate_refit_original_weights,
      matched_fixed_nuisance = list(
        estimate_original_weights =
          row$estimate_fixed_nuisance_original_weights,
        estimate_relearned_weights =
          row$estimate_fixed_nuisance_relearned_weights
      ),
      fitted = list(
        estimate = row$estimate_refit_relearned_weights,
        se = reference_se,
        target_only = list(estimate = row$estimate_target_refit)
      ),
      nuisance_cv_grouping = if (protocol == "origin_grouped") "origin" else NULL,
      nuisance_cv_partitions = if (protocol == "origin_grouped") {
        if (is.null(partition_records)) list(list(
          caller = "synthetic", n_folds = 3L,
          cv_group_id = c(1L, 1L, 2L, 2L, 3L, 3L),
          cv_fold_id = c(1L, 1L, 2L, 2L, 3L, 3L), valid = TRUE
        )) else partition_records
      } else NULL
    )
  } else {
    list(
      failure_message = row$failure_message,
      failure_stage = row$failure_stage
    )
  }
  saved_summary <- row
  saved_summary$estimate_reference <-
    saved_summary$estimate_reference + rds_shift
  saveRDS(list(
    result = result, summary = saved_summary, warnings = character(),
    source_bundle = source_bundle
  ), file.path(bundle, "draw.rds"))
  files <- c("result.csv", "draw.rds")
  hashes <- vapply(file.path(bundle, files), provenance$roce_sha256_file,
                   character(1L))
  writeLines(paste(hashes, files, sep = "  "),
             file.path(bundle, "sha256.txt"))
  invisible(bundle)
}

root <- make_root()
write_draw(root, "draw0", 0L)
write_draw(root, "draw1", 1L)
write_draw(root, "draw2", 2L)
output <- tempfile("full_refit_summary_success_")
summarize_full_refit_draws_main(c(root, "1", "0", "0:2", output))
summary <- read.csv(file.path(output, "resampling_summary.csv"))
audit <- read.csv(file.path(output, "draw_audit.csv"), na.strings = c("", "NA"))
stopifnot(
  nrow(summary) == 1L, nrow(audit) == 3L,
  identical(audit$draw_id, 0:2),
  summary$n_attempted_nonidentity == 2L,
  summary$n_successful_nonidentity == 2L,
  summary$n_failed_nonidentity == 0L,
  summary$diagnostic_draw_set_complete,
  !summary$inference_validated,
  summary$small_B_warning,
  summary$paired_identity_error < 1e-12
)
old <- setwd(output)
checksum <- system2("sha256sum", c("-c", "sha256.txt"),
                    stdout = TRUE, stderr = TRUE)
setwd(old)
stopifnot(is.null(attr(checksum, "status")), all(grepl("OK$", checksum)))

failed_root <- make_root()
write_draw(failed_root, "draw0", 0L)
write_draw(failed_root, "draw1", 1L)
write_draw(failed_root, "draw2", 2L, status = "failed")
failed_output <- tempfile("full_refit_summary_failed_")
summarize_full_refit_draws_main(c(
  failed_root, "1", "0", "0:2", failed_output
))
failed_summary <- read.csv(file.path(failed_output, "resampling_summary.csv"))
failed_audit <- read.csv(file.path(failed_output, "draw_audit.csv"),
                         na.strings = c("", "NA"))
stopifnot(
  failed_summary$n_failed_nonidentity == 1L,
  !failed_summary$diagnostic_draw_set_complete,
  !failed_summary$inference_validated,
  failed_audit$status[failed_audit$draw_id == 2L] == "failed",
  failed_audit$failure_message[failed_audit$draw_id == 2L] ==
    "synthetic failure"
)

missing_output <- tempfile("full_refit_summary_missing_")
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  root, "1", "0", "0:3", missing_output
)), silent = TRUE), "try-error"), !file.exists(missing_output))

mixed_root <- make_root()
write_draw(mixed_root, "draw0", 0L)
write_draw(mixed_root, "draw1", 1L)
write_draw(mixed_root, "draw2", 2L, workflow = strrep("d", 64L))
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  mixed_root, "1", "0", "0:2", tempfile("mixed_output_")
)), silent = TRUE), "try-error"))

duplicate_root <- make_root()
write_draw(duplicate_root, "draw0", 0L)
write_draw(duplicate_root, "draw1a", 1L)
write_draw(duplicate_root, "draw1b", 1L)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  duplicate_root, "1", "0", "0:1", tempfile("duplicate_output_")
)), silent = TRUE), "try-error"))

identity_root <- make_root()
write_draw(identity_root, "draw0", 0L, identity_shift = 0.01)
write_draw(identity_root, "draw1", 1L)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  identity_root, "1", "0", "0:1", tempfile("identity_output_")
)), silent = TRUE), "try-error"))

disagreement_root <- make_root()
write_draw(disagreement_root, "draw0", 0L)
write_draw(disagreement_root, "draw1", 1L, rds_shift = 0.01)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  disagreement_root, "1", "0", "0:1", tempfile("disagreement_output_")
)), silent = TRUE), "try-error"))

nonfinite_se_root <- make_root()
write_draw(nonfinite_se_root, "draw0", 0L, reference_se = Inf)
write_draw(nonfinite_se_root, "draw1", 1L, reference_se = Inf)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  nonfinite_se_root, "1", "0", "0:1", tempfile("nonfinite_se_output_")
)), silent = TRUE), "try-error"))

# Exercise both decompositions with a varying baseline and changes that are
# correlated with it; constant-change examples cannot detect omitted baseline
# covariance terms.
A <- c(0, 1, 3, 4)
weight <- c(0.2, 0.8, -0.1, 1.5)
nuisance <- c(-0.3, 0.4, 1.1, 0.2)
interaction <- c(0.1, -0.2, 0.5, -0.4)
decomposition_fixture <- data.frame(
  estimate_fixed_nuisance_original_weights = A,
  estimate_fixed_nuisance_relearned_weights = A + weight,
  estimate_refit_original_weights = A + nuisance,
  estimate_refit_relearned_weights = A + weight + nuisance + interaction
)
decomposition <- .frs_covariance_decomposition(decomposition_fixture)
value <- function(scope, term) decomposition$value[
  decomposition$decomposition_scope == scope & decomposition$term == term
]
stopifnot(
  abs(value("change_D_minus_A", "change_decomposition_error")) < 1e-12,
  abs(value("full_D", "full_decomposition_error")) < 1e-12,
  value("full_D", "variance_baseline_A") > 0,
  abs(value("full_D", "twice_cov_A_weight")) > 0,
  abs(value("full_D", "reconstructed_variance_D") -
      stats::var(decomposition_fixture$estimate_refit_relearned_weights)) <
    1e-12
)
one_draw <- .frs_covariance_decomposition(decomposition_fixture[1L, ])
zero_draw <- .frs_covariance_decomposition(
  decomposition_fixture[FALSE, , drop = FALSE]
)
stopifnot(all(is.na(one_draw$value)), all(is.na(zero_draw$value)))

grouped_root <- make_root()
write_draw(grouped_root, "draw0", 0L, protocol = "origin_grouped")
write_draw(grouped_root, "draw1", 1L, protocol = "origin_grouped")
grouped_output <- tempfile("full_refit_summary_grouped_")
summarize_full_refit_draws_main(c(
  grouped_root, "1", "0", "0:1", grouped_output
))
grouped_summary <- read.csv(file.path(
  grouped_output, "resampling_summary.csv"
))
grouped_audit <- read.csv(file.path(grouped_output, "draw_audit.csv"))
stopifnot(
  identical(grouped_summary$nuisance_cv_partition, "origin_grouped"),
  all(grouped_audit$nuisance_cv_partition == "origin_grouped"),
  all(grouped_audit$nuisance_cv_partition_records == 1L),
  all(grouped_audit$nuisance_cv_explicit_partition_records == 1L)
)

crossing_record <- list(list(
  caller = "synthetic", n_folds = 3L,
  cv_group_id = c(1L, 1L, 2L, 2L, 3L, 3L),
  cv_fold_id = c(1L, 2L, 2L, 2L, 3L, 3L), valid = TRUE
))
crossing_root <- make_root()
write_draw(crossing_root, "draw0", 0L, protocol = "origin_grouped")
write_draw(crossing_root, "draw1", 1L, protocol = "origin_grouped",
           partition_records = crossing_record)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  crossing_root, "1", "0", "0:1", tempfile("crossing_partition_")
)), silent = TRUE), "try-error"))

empty_partition_root <- make_root()
write_draw(empty_partition_root, "draw0", 0L, protocol = "origin_grouped")
write_draw(empty_partition_root, "draw1", 1L, protocol = "origin_grouped",
           partition_records = list())
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  empty_partition_root, "1", "0", "0:1", tempfile("empty_partition_")
)), silent = TRUE), "try-error"))

false_valid_record <- list(list(
  caller = "synthetic", n_folds = 3L,
  cv_group_id = c(1L, 1L, 2L, 2L, 3L, 3L),
  cv_fold_id = c(1L, 1L, 2L, 2L, 3L, 3L), valid = FALSE
))
false_valid_root <- make_root()
write_draw(false_valid_root, "draw0", 0L, protocol = "origin_grouped")
write_draw(false_valid_root, "draw1", 1L, protocol = "origin_grouped",
           partition_records = false_valid_record)
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  false_valid_root, "1", "0", "0:1", tempfile("false_valid_partition_")
)), silent = TRUE), "try-error"))

mixed_protocol_root <- make_root()
write_draw(mixed_protocol_root, "draw0", 0L, protocol = "row_level")
write_draw(mixed_protocol_root, "draw1", 1L, protocol = "origin_grouped")
stopifnot(inherits(try(summarize_full_refit_draws_main(c(
  mixed_protocol_root, "1", "0", "0:1", tempfile("mixed_protocol_")
)), silent = TRUE), "try-error"))

cat("summarize_full_refit_draws synthetic tests: PASS\n")
