#!/usr/bin/env Rscript
# ============================================================================
# diagnosis/bias/summarize.R
# ----------------------------------------------------------------------------
# Aggregates all per-setting *_summary.csv files in diagnosis/bias/results/
# into a single long-format table grouped by axis (n / p / K / shift).
# Focuses on the four most relevant methods for the bias question:
#   - one_round_crossfit   (our method, single-round)
#   - two_round_crossfit   (our method, two-round)
#   - federated_dr         (baseline; ~unbiased in C1)
#   - oracle_dr            (oracle DR; sanity reference)
#
# Output: diagnosis/bias/results/_summary_table.csv  +  console table.
#
# Usage:
#   Rscript diagnosis/bias/summarize.R
# ============================================================================

out_dir <- "diagnosis/bias/results"
files <- list.files(out_dir, pattern = "_summary\\.csv$", full.names = TRUE)
files <- files[!grepl("_summary_table\\.csv$", files)]

if (length(files) == 0L) {
  stop(sprintf("summarize.R: no *_summary.csv files in %s.", out_dir))
}

target_methods <- c("one_round_crossfit", "two_round_crossfit",
                    "federated_dr", "oracle_dr", "pooled_dr")

rows <- list()
for (f in files) {
  df <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) {
    warning(sprintf("summarize.R: skipping unreadable/empty %s", f))
    next
  }
  keep <- df$method %in% target_methods
  if (!any(keep)) {
    warning(sprintf("summarize.R: %s has no rows for target methods.", f))
    next
  }
  rows[[length(rows) + 1L]] <- df[keep, , drop = FALSE]
}

if (length(rows) == 0L) {
  stop("summarize.R: no usable rows after filtering for target methods.")
}

all_rows <- do.call(rbind, rows)

# Derive an axis label from the tag for easier reading.
axis_of <- function(tag) {
  if (startsWith(tag, "n_") || tag == "anchor") "n"     else
  if (startsWith(tag, "p_")) "p"                         else
  if (startsWith(tag, "K_")) "K"                         else
  if (startsWith(tag, "shift_")) "shift"                 else NA_character_
}
all_rows$axis <- vapply(all_rows$tag, axis_of, character(1))

select_cols <- c("axis", "tag", "method",
                 "n_total", "K", "shift",
                 "bias_mean", "bias_sd", "rmse",
                 "se_mean", "coverage", "n_success", "elapsed_min")
present_cols <- intersect(select_cols, names(all_rows))
out <- all_rows[, present_cols, drop = FALSE]

# Sort: axis, then method, then by the varied parameter inside each axis.
out <- out[order(out$axis, out$tag, out$method), , drop = FALSE]

table_path <- file.path(out_dir, "_summary_table.csv")
write.csv(out, file = table_path, row.names = FALSE)
cat(sprintf("[summarize] wrote %s (%d rows)\n", table_path, nrow(out)))

cat("\n========================  diagnosis/bias summary  ========================\n")
print_axis <- function(axis_name) {
  sub <- out[out$axis == axis_name, , drop = FALSE]
  if (nrow(sub) == 0L) return(invisible())
  cat(sprintf("\n---- axis = %s ----\n", axis_name))
  show_cols <- c("tag", "method", "n_total", "K", "shift",
                 "bias_mean", "bias_sd", "se_mean", "coverage", "n_success")
  show_cols <- intersect(show_cols, names(sub))
  fmt <- sub
  for (col in c("bias_mean", "bias_sd", "se_mean", "coverage")) {
    if (col %in% names(fmt)) fmt[[col]] <- signif(as.numeric(fmt[[col]]), 4)
  }
  print(fmt[, show_cols], row.names = FALSE)
}
for (ax in c("n", "p", "K", "shift")) print_axis(ax)
cat("\n==========================================================================\n")
