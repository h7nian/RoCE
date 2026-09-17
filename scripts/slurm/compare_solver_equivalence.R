#!/usr/bin/env Rscript
# Compare paired replicate outputs produced under the two nuisance solvers
# (ROCE_NUISANCE_SOLVER=coordinate_descent vs the proximal-Newton default):
# same seed, config, K and fold count. Expects <ROOT>/cd/raw and <ROOT>/newton/raw.
#
# usage: compare_solver_equivalence.R [ROOT]   (default results/.../solver_equiv_20260917)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else
  "results/direct_tate_mc500_b5000/solver_equiv_20260917"
a_files <- list.files(file.path(root, "cd", "raw"), pattern = "\\.csv$", full.names = TRUE)
if (length(a_files) == 0L) stop("no cd outputs under ", root, call. = FALSE)

report_methods <- c("one_round_crossfit_ate", "target_only_ate",
                    "one_round_crossfit_ate_quadratic_bias")
key_columns <- c("estimate", "se", "bias", "coverage", "ci_width")

for (a_file in a_files) {
  b_file <- file.path(root, "newton", "raw", basename(a_file))
  if (!file.exists(b_file)) { cat(basename(a_file), ": newton output missing\n"); next }
  a <- read.csv(a_file, stringsAsFactors = FALSE)
  b <- read.csv(b_file, stringsAsFactors = FALSE)
  cat("====", basename(a_file), "====\n")
  cat(sprintf("elapsed: cd %.0fs  newton %.0fs  (%.1fx)\n",
              a$pilot_elapsed_seconds[[1L]], b$pilot_elapsed_seconds[[1L]],
              a$pilot_elapsed_seconds[[1L]] / b$pilot_elapsed_seconds[[1L]]))
  stage_columns <- intersect(
    c("data_generation_seconds", "one_round_face_wall_seconds",
      "tate_contrast_stage_seconds", "comparison_mu1_seconds",
      "simulation_elapsed_seconds"),
    intersect(names(a), names(b))
  )
  for (sc in stage_columns) {
    cat(sprintf("  %-32s cd %8.0fs  newton %8.0fs\n", sc,
                as.numeric(a[[sc]][[1L]]), as.numeric(b[[sc]][[1L]])))
  }
  for (m in intersect(report_methods, a$method)) {
    ra <- a[a$method == m, ]; rb <- b[b$method == m, ]
    cat(sprintf("  %-40s", m))
    for (k in key_columns) {
      cat(sprintf(" %s: %.6g|%.6g", k, as.numeric(ra[[k]]), as.numeric(rb[[k]])))
    }
    cat("\n")
  }
  shared <- intersect(names(a), names(b))
  numeric_cols <- shared[vapply(shared, function(cc) is.numeric(a[[cc]]) && is.numeric(b[[cc]]), logical(1L))]
  numeric_cols <- setdiff(numeric_cols, grep("seconds|elapsed|allocated|threads|workers|cores", shared, value = TRUE))
  diffs <- vapply(numeric_cols, function(cc) {
    x <- as.numeric(a[[cc]]); y <- as.numeric(b[[cc]])
    if (length(x) != length(y)) return(NA_real_)
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) return(0)
    max(abs(x[ok] - y[ok]) / pmax(1, abs(x[ok])))
  }, numeric(1L))
  diffs <- sort(diffs[is.finite(diffs)], decreasing = TRUE)
  cat(sprintf("  numeric columns compared: %d; columns with relative diff > 1e-8: %d\n",
              length(diffs), sum(diffs > 1e-8)))
  top <- head(diffs[diffs > 1e-8], 15)
  if (length(top)) for (nm in names(top)) cat(sprintf("    %-60s %.3e\n", nm, top[[nm]]))
}
