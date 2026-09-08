#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(
    "usage: build_direct_tate_rho_group_manifest.R PRIMARY_MANIFEST.csv OUTPUT.csv",
    call. = FALSE
  )
}

primary_path <- normalizePath(args[[1L]], mustWork = TRUE)
output_path <- args[[2L]]
if (file.exists(output_path)) {
  stop("refusing to overwrite existing rho-group manifest: ", output_path,
       call. = FALSE)
}
primary <- read.csv(primary_path, stringsAsFactors = FALSE)
rho_grid <- c(0, 0.5, 1, 1.5, 2, 2.5)
# deviation_mechanism stays last so the grouped column positions read by
# submit_rho_group_direct_tate.sh are unchanged.
required <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
  "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap", "M_tau",
  "M_tau_inference", "methods", "deviation_mechanism"
)
missing <- setdiff(required, names(primary))
if (length(missing) > 0L) {
  stop("primary manifest is missing: ", paste(missing, collapse = ", "),
       call. = FALSE)
}
grouped_experiments <- c("negative_transfer", "shared_shift")
if (nrow(primary) == 0L || any(primary$p != 100L) ||
    length(unique(primary$experiment)) != 1L ||
    !primary$experiment[[1L]] %in% grouped_experiments) {
  stop(
    "rho grouping is restricted to one p=100 negative-transfer or ",
    "shared-shift manifest.", call. = FALSE
  )
}

group_key <- interaction(
  primary$config, primary$K, primary$sim_id,
  drop = TRUE, lex.order = TRUE
)
groups <- split(primary, group_key)
stable_columns <- setdiff(required, c("task_id", "rho"))
rows <- lapply(groups, function(group) {
  group <- group[order(group$rho), , drop = FALSE]
  if (!identical(as.numeric(group$rho), rho_grid)) {
    stop(
      "each config/K/sim group must contain rho={0,0.5,1,1.5,2,2.5}.",
      call. = FALSE
    )
  }
  stable <- vapply(stable_columns, function(column) {
    length(unique(group[[column]])) == 1L
  }, logical(1L))
  if (!all(stable)) {
    stop(
      "rho group differs in invariant manifest field(s): ",
      paste(names(stable)[!stable], collapse = ", "), ".", call. = FALSE
    )
  }
  row <- group[1L, stable_columns, drop = FALSE]
  row$rho_values <- paste(format(rho_grid, trim = TRUE), collapse = ";")
  row$primary_task_ids <- paste(as.integer(group$task_id), collapse = ";")
  row
})
grouped <- do.call(rbind, rows)
grouped <- grouped[order(grouped$config, grouped$K, grouped$sim_id), , drop = FALSE]
grouped$group_task_id <- seq_len(nrow(grouped))
grouped <- grouped[c("group_task_id", setdiff(names(grouped), "group_task_id"))]
rownames(grouped) <- NULL

expected_groups <- nrow(primary) / length(rho_grid)
if (expected_groups %% 1 != 0 || nrow(grouped) != expected_groups ||
    anyDuplicated(unlist(strsplit(grouped$primary_task_ids, ";", fixed = TRUE))) ||
    length(unlist(strsplit(grouped$primary_task_ids, ";", fixed = TRUE))) !=
      nrow(primary)) {
  stop("rho-group manifest is not an exact partition of the primary manifest.",
       call. = FALSE)
}

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary <- tempfile(
  pattern = ".rho_group_manifest_", tmpdir = dirname(output_path),
  fileext = ".csv"
)
write.csv(grouped, temporary, row.names = FALSE)
if (!file.rename(temporary, output_path)) {
  unlink(temporary)
  stop("failed to atomically write rho-group manifest.", call. = FALSE)
}
message(sprintf(
  "[pass] wrote %d rho-group tasks covering %d p=100 primary rows",
  nrow(grouped), nrow(primary)
))
