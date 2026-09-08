#!/usr/bin/env Rscript
# Static method/reference alignment audit for RoCE.
#
# This script is intentionally read-only. It checks that docs/main.tex follows
# the local FACE/SMMAL references and that the current implementation follows
# docs/main.tex for aggregation cross-validation and two-layer cross-fitting.

args <- commandArgs(trailingOnly = TRUE)
repo_root <- if (length(args) >= 1L && nzchar(args[[1L]])) args[[1L]] else "."
repo_root <- normalizePath(repo_root, mustWork = TRUE)

output_root <- Sys.getenv(
  "METHOD_AUDIT_OUTPUT_ROOT",
  file.path("diagnosis", "method_alignment", "audit_results")
)
output_dir <- file.path(repo_root, output_root)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

job_id <- Sys.getenv("SLURM_JOB_ID", "local")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_prefix <- file.path(output_dir, sprintf("%s_job%s", timestamp, job_id))

path_in_repo <- function(...) file.path(repo_root, ...)

read_text <- function(path) {
  if (!file.exists(path)) return("")
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

compact_text <- function(text) gsub("[[:space:]]+", " ", text)
nospace_text <- function(text) gsub("[[:space:]]+", "", text)

extract_pdf_text <- function(path) {
  if (!file.exists(path)) return("")
  extracted_path <- sub("\\.pdf$", ".extracted.txt", path, ignore.case = TRUE)
  if (file.exists(extracted_path)) {
    return(read_text(extracted_path))
  }
  gs <- Sys.which("gs")
  if (!nzchar(gs) && file.exists("/usr/bin/gs")) gs <- "/usr/bin/gs"
  if (!nzchar(gs)) return("")
  out <- tryCatch(
    system2(
      gs,
      c("-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
        "-sOutputFile=-", path),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) character()
  )
  paste(out, collapse = "\n")
}

has_regex <- function(text, pattern) {
  grepl(pattern, text, perl = TRUE, ignore.case = TRUE)
}

rows <- list()
add_row <- function(area, check, status, severity, file, detail) {
  rows[[length(rows) + 1L]] <<- data.frame(
    area = area,
    check = check,
    status = status,
    severity = severity,
    file = file,
    detail = detail,
    stringsAsFactors = FALSE
  )
}

check_pattern <- function(area, check, text, pattern, file,
                          detail_ok, detail_bad,
                          severity = "critical",
                          compact = FALSE, nospace = FALSE) {
  target <- if (nospace) nospace_text(text) else if (compact) compact_text(text) else text
  ok <- has_regex(target, pattern)
  add_row(
    area, check,
    if (ok) "pass" else "fail",
    severity,
    file,
    if (ok) detail_ok else detail_bad
  )
  invisible(ok)
}

method_path <- if (file.exists(path_in_repo("method.tex"))) {
  path_in_repo("method.tex")
} else {
  path_in_repo("docs", "main.tex")
}
method_rel <- if (basename(method_path) == "method.tex") "method.tex" else "docs/main.tex"

main_tex <- read_text(method_path)
face_text <- extract_pdf_text(path_in_repo("docs", "FACE.pdf"))
smmal_text <- extract_pdf_text(path_in_repo("docs", "SMMAL.pdf"))
agg_code <- read_text(path_in_repo("R", "cross_fitting_aggregation.R"))
cf_code <- read_text(path_in_repo("R", "cross_fitting_algorithms.R"))
fit_code <- read_text(path_in_repo("R", "model_fitting.R"))
cv_code <- read_text(path_in_repo("src", "cv_utils.h"))
dr_code <- read_text(path_in_repo("src", "density_ratio.cpp"))
oracle_code <- read_text(path_in_repo("R", "estimators_oracle.R"))

add_row(
  "method_file",
  "method_source",
  if (nzchar(main_tex)) "pass" else "fail",
  if (nzchar(main_tex)) "info" else "critical",
  method_rel,
  if (identical(method_rel, "docs/main.tex")) {
    "No repository-level method.tex was found; using docs/main.tex as the method source."
  } else {
    "Found repository-level method.tex."
  }
)

if (!nzchar(face_text)) {
  add_row(
    "FACE_reference", "pdf_text_extraction", "skipped", "warning",
    "docs/FACE.pdf",
    "FACE.pdf exists, but text extraction was unavailable on this Slurm node."
  )
} else {
  check_pattern(
    "FACE_reference", "sample_split_validation", face_text,
    "sample[[:space:]]*splitting|trainingandvalidationdatasets",
    "docs/FACE.pdf",
    "FACE reference contains a sample-splitting training/validation rule for aggregation tuning.",
    "Could not find FACE sample-splitting validation language in extracted PDF text.",
    severity = "warning", nospace = TRUE
  )
  check_pattern(
    "FACE_reference", "validation_summaries", face_text,
    "validationsummaries|validationdatasets.*summarystatistics|summarystatistics.*validationdatasets",
    "docs/FACE.pdf",
    "FACE reference evaluates training-set weights using validation-set summaries.",
    "Could not find FACE validation-summary language in extracted PDF text.",
    severity = "warning", nospace = TRUE
  )
  check_pattern(
    "FACE_reference", "five_fold_cv", face_text,
    "5-foldcross-validation|5foldcrossvalidation",
    "docs/FACE.pdf",
    "FACE reference states a 5-fold cross-validation implementation for lambda.",
    "Could not find FACE 5-fold CV language in extracted PDF text.",
    severity = "info", nospace = TRUE
  )
}

if (!nzchar(smmal_text)) {
  add_row(
    "SMMAL_reference", "pdf_text_extraction", "skipped", "warning",
    "docs/SMMAL.pdf",
    "SMMAL.pdf exists, but text extraction was unavailable on this Slurm node."
  )
} else {
  check_pattern(
    "SMMAL_reference", "two_level_crossfit", smmal_text,
    "two[-[:space:]]*(level|layer)[[:space:]]*cross|two(level|layer)cross",
    "docs/SMMAL.pdf",
    "SMMAL reference describes the two-layer/two-level cross-fitting idea.",
    "Could not find two-layer/two-level cross-fitting language in extracted SMMAL text.",
    severity = "warning", nospace = TRUE
  )
  check_pattern(
    "SMMAL_reference", "out_of_two_folds", smmal_text,
    "out-of-two-fold|outoftwofold",
    "docs/SMMAL.pdf",
    "SMMAL reference trains initial nuisance fits on out-of-two-fold data.",
    "Could not find out-of-two-fold initial fitting language in extracted SMMAL text.",
    severity = "warning", nospace = TRUE
  )
  check_pattern(
    "SMMAL_reference", "calibrated_losses", smmal_text,
    "calibratedloss|calibratedlosses",
    "docs/SMMAL.pdf",
    "SMMAL reference defines calibrated losses after initial nuisance estimation.",
    "Could not find calibrated-loss language in extracted SMMAL text.",
    severity = "warning", nospace = TRUE
  )
}

check_pattern(
  "main_tex", "face_style_validation", main_tex,
  "FACE-style validation for the aggregation penalty",
  method_rel,
  "docs/main.tex contains the FACE-style aggregation validation subsection.",
  "docs/main.tex is missing the FACE-style aggregation validation subsection.",
  compact = TRUE
)
check_pattern(
  "main_tex", "smmal_fold_summed", main_tex,
  "Fold-summed SMMAL-style calibration",
  method_rel,
  "docs/main.tex contains the fold-summed SMMAL-style calibration subsection.",
  "docs/main.tex is missing the fold-summed SMMAL-style calibration subsection.",
  compact = TRUE
)
check_pattern(
  "main_tex", "inner_training_validation", main_tex,
  "T_\\{k_1,-m\\}|inner-training|inner validation fold",
  method_rel,
  "docs/main.tex separates inner-training summaries from validation-fold scoring.",
  "docs/main.tex does not clearly separate inner training and validation folds.",
  compact = TRUE
)

check_pattern(
  "implementation", "aggregation_inner_cv_selector", agg_code,
  "select_aggregation_lambda_inner_cv",
  "R/cross_fitting_aggregation.R",
  "Aggregation implementation uses the inner-CV selector.",
  "Aggregation implementation does not contain select_aggregation_lambda_inner_cv()."
)
check_pattern(
  "implementation", "aggregation_train_val_split", agg_code,
  "components\\[-m\\][\\s\\S]*components\\[\\[m\\]\\]",
  "R/cross_fitting_aggregation.R",
  "Aggregation lambda CV trains weights on components[-m] and validates on components[[m]].",
  "Aggregation lambda CV does not visibly train on components[-m] and validate on components[[m]]."
)
check_pattern(
  "implementation", "aggregation_validation_objective", agg_code,
  "n_val \\* var_term|n_val\\*var_term",
  "R/cross_fitting_aggregation.R",
  "Validation objective multiplies validation variance by validation sample size.",
  "Validation objective does not show the FACE N_val-stabilized criterion."
)
check_pattern(
  "implementation", "aggregation_final_outer_training_weights", agg_code,
  "optimize_weights\\([\\s\\S]*inner_v\\$avg_source_est[\\s\\S]*inner_v\\$avg_target_est",
  "R/cross_fitting_aggregation.R",
  "Final outer-fold weights are recomputed from averaged inner training components.",
  "Final outer-fold weights are not visibly recomputed from inner training components."
)

check_pattern(
  "implementation", "secondary_folds_exclude_outer", cf_code,
  "secondary_folds <- setdiff\\(1:n_folds, k1\\)",
  "R/cross_fitting_algorithms.R",
  "Two-layer cross-fitting excludes the outer fold from secondary calibration folds.",
  "Could not find secondary_folds <- setdiff(1:n_folds, k1)."
)
check_pattern(
  "implementation", "initial_fits_out_of_two_folds", cf_code,
  "training_folds <- setdiff\\(1:n_folds, c\\(k1, k2\\)\\)",
  "R/cross_fitting_algorithms.R",
  "Initial nuisance fits use training folds excluding both k1 and k2.",
  "Could not find out-of-two-fold training_folds construction."
)
check_pattern(
  "implementation", "fold_summed_block_design", cf_code,
  "\\.make_plugin_block_design[\\s\\S]*alpha_plugin_block[\\s\\S]*gamma_plugin_block",
  "R/cross_fitting_algorithms.R",
  "Fold-specific plug-ins are represented with block designs in one calibrated optimization.",
  "Could not find block-design plug-in construction for fold-summed calibration."
)
check_pattern(
  "implementation", "calibrated_nuisance_fits", cf_code,
  "fit_unified_density_ratio[\\s\\S]*calibrated = TRUE[\\s\\S]*fit_unified_outcome[\\s\\S]*calibrated = TRUE",
  "R/cross_fitting_algorithms.R",
  "Both final gamma and alpha use calibrated nuisance losses.",
  "Could not confirm calibrated=TRUE for both final gamma and final alpha."
)

check_pattern(
  "nuisance_cv", "calibrated_gamma_cv_function", fit_code,
  "select_lambda_cv_calibrated_density_ratio_cpp",
  "R/model_fitting.R",
  "Calibrated gamma wrapper uses the calibrated density-ratio CV routine.",
  "Calibrated gamma wrapper does not call select_lambda_cv_calibrated_density_ratio_cpp()."
)
check_pattern(
  "nuisance_cv", "calibrated_alpha_cv_function", fit_code,
  "select_lambda_cv_calibrated_outcome_cpp",
  "R/model_fitting.R",
  "Calibrated alpha wrapper uses the calibrated outcome CV routine.",
  "Calibrated alpha wrapper does not call select_lambda_cv_calibrated_outcome_cpp()."
)

dr_training_full_scale <- has_regex(
  dr_code,
  "grad_acc\\.add\\([^\\n]+/ n\\)"
) || has_regex(cv_code, "grad_acc\\.add\\([^\\n]+/ n_total\\)")
add_row(
  "nuisance_cv",
  "density_ratio_training_full_source_scale",
  if (dr_training_full_scale) "pass" else "fail",
  if (dr_training_full_scale) "info" else "warning",
  "src/density_ratio.cpp; src/cv_utils.h",
  if (dr_training_full_scale) {
    "Density-ratio training gradient is on the full-source empirical scale."
  } else {
    "Could not confirm full-source empirical scaling in density-ratio training gradients."
  }
)

dr_val_has_source_scale <- has_regex(cv_code, "density_ratio_val_loss\\([^\\)]*(source_scale|arm_fraction|n_source|n_total)")
dr_val_uses_arm_mean <- grepl("exp_neg_g * psi_prime_val(j) / n_val", cv_code, fixed = TRUE)
add_row(
  "nuisance_cv",
  "density_ratio_validation_scale",
  if (dr_val_has_source_scale && !dr_val_uses_arm_mean) "pass" else "needs_followup",
  "warning",
  "src/cv_utils.h",
  if (dr_val_has_source_scale && !dr_val_uses_arm_mean) {
    "Density-ratio validation loss explicitly carries the source empirical scale."
  } else {
    paste(
      "Density-ratio validation appears to average over the treated validation arm.",
      "This is not the same empirical scale as E_s[I(A=a)f] in main.tex and should remain under diagnosis."
    )
  }
)

oracle_uses_raw_gamma <- has_regex(oracle_code, "gamma_true <- gamma_params\\[\\[gamma_key\\]\\][\\s\\S]*w_i <- exp\\(-eta_dr\\)")
add_row(
  "oracle_benchmark",
  "oracle_gamma_convention",
  if (oracle_uses_raw_gamma) "needs_source_fix_candidate" else "pass",
  "warning",
  "R/estimators_oracle.R",
  if (oracle_uses_raw_gamma) {
    paste(
      "Oracle DR currently uses raw DGP joint gamma directly.",
      "Earlier diagnosis shows source-conditional intercept calibration improves oracle alignment."
    )
  } else {
    "Oracle DR no longer appears to use raw DGP joint gamma directly."
  }
)

audit <- if (length(rows)) do.call(rbind, rows) else data.frame()
csv_path <- paste0(output_prefix, "_method_alignment_audit.csv")
write.csv(audit, csv_path, row.names = FALSE)

summary <- aggregate(
  check ~ area + status + severity,
  data = audit,
  FUN = length
)
names(summary)[names(summary) == "check"] <- "n_checks"
summary_path <- paste0(output_prefix, "_method_alignment_summary.csv")
write.csv(summary, summary_path, row.names = FALSE)

cat(sprintf("[method-audit] wrote %s\n", csv_path))
cat(sprintf("[method-audit] wrote %s\n", summary_path))
print(summary, row.names = FALSE)

critical_fail <- audit$status == "fail" & audit$severity == "critical"
if (any(critical_fail)) {
  stop(sprintf(
    "method_alignment_audit: %d critical alignment checks failed; inspect %s",
    sum(critical_fail), csv_path
  ), call. = FALSE)
}

cat("[method-audit] DONE\n")
