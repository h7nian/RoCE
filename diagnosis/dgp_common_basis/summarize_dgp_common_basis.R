#!/usr/bin/env Rscript
# HISTORY: 2026-09-07 #0002 Common working basis and X-dagger misspecification: strength pilot
# Task:    Collect the 13 pilot cells (diagnosis/out/dgp_common_basis/<cell>/), read the
#          saved fits for the inference truncation and convergence diagnostics, and
#          apply the pre-registered omega selection rule of HISTORY #0002 §4.
#          Writes diagnosis/out/dgp_common_basis/selection/.

suppressPackageStartupMessages(library(RoCE))
source("scripts/slurm/result_provenance.R")
source("scripts/slurm/atomic_output.R")

CELL_ROOT <- "diagnosis/out/dgp_common_basis"
R_SQUARED_RANGE <- c(0.60, 0.90)
TRUNCATION_FRACTION_MAX <- 0.01
# The C1 mechanism itself reaches a raw density-ratio weight of about 29 at radius
# 5 with no truncation (HISTORY #0002, Rule 26 gate re-evaluation), so the
# absolute bound of 20 was replaced by a bound relative to the C1 value.
MAX_RAW_WEIGHT_RELATIVE_TO_C1 <- 1.5
# The DGP clips true propensities to [POSITIVITY_LOWER, POSITIVITY_UPPER] in every
# configuration, so criterion (c) is read as: no more clipping than the C1 mechanism.
PROPENSITY_CLIP_FRACTION_MARGIN <- 0.01
C1_DIFFERENCE_MAX <- 0.002

.nonconverged_count <- function(diagnostics) {
  if (!is.data.frame(diagnostics)) return(NA_integer_)
  flags <- diagnostics[grepl("_converged$", names(diagnostics))]
  sum(vapply(flags, function(column) sum(!as.logical(column)), integer(1L)))
}

.cell_summary <- function(directory) {
  summary <- read.csv(file.path(directory, "summary.csv"))
  fit <- readRDS(file.path(directory, "fit.rds"))$fit
  clip <- fit$clip_diagnostics$total
  summary$inference_logit_truncation_fraction <- clip$logit_truncation_fraction
  summary$max_raw_weight <- clip$max_raw_weight
  summary$max_abs_logit <- clip$max_abs_logit
  summary$nonconverged_fits <- .nonconverged_count(fit$nuisance_fit_diagnostics)
  summary$misspecified_r_squared <- switch(
    as.character(summary$config),
    C1 = NA_real_,
    C2 = summary$r_squared_outcome,
    C3 = summary$r_squared_propensity,
    C4 = min(summary$r_squared_outcome, summary$r_squared_propensity)
  )
  summary
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  output <- file.path(CELL_ROOT, "selection")
  if (file.exists(output)) stop("output already exists: ", output)
  cells <- list.dirs(CELL_ROOT, recursive = FALSE, full.names = TRUE)
  cells <- cells[file.exists(file.path(cells, "summary.csv"))]
  table <- do.call(rbind, lapply(cells, .cell_summary))
  table <- table[order(table$config, table$omega), ]
  table$r_squared_ok <- is.na(table$misspecified_r_squared) |
    (table$misspecified_r_squared >= R_SQUARED_RANGE[1] &
       table$misspecified_r_squared <= R_SQUARED_RANGE[2])
  c1 <- table[table$config == "C1", ]
  if (nrow(c1) != 1L) stop("the C1 cell is required as the positivity and weight reference")
  table$truncation_ok <- table$inference_logit_truncation_fraction <= TRUNCATION_FRACTION_MAX &
    table$max_raw_weight <= MAX_RAW_WEIGHT_RELATIVE_TO_C1 * c1$max_raw_weight
  c1_clip_fraction <- c1$propensity_clip_fraction
  table$propensity_ok <- table$propensity_clip_fraction <= c1_clip_fraction + PROPENSITY_CLIP_FRACTION_MARGIN
  table$all_ok <- table$r_squared_ok & table$truncation_ok & table$propensity_ok
  # omega* = largest omega passing every criterion for both C2 and C3.
  candidates <- sort(unique(table$omega[table$config %in% c("C2", "C3")]), decreasing = TRUE)
  omega_star <- NA_real_
  for (omega in candidates) {
    rows <- table[table$config %in% c("C2", "C3") & table$omega == omega, ]
    if (nrow(rows) == 2L && all(rows$all_ok)) { omega_star <- omega; break }
  }
  c1_ok <- abs(c1$c1_reference_difference) <= C1_DIFFERENCE_MAX
  print(table[, c("config", "omega", "misspecified_r_squared", "inference_logit_truncation_fraction",
                  "max_raw_weight", "propensity_clip_fraction", "nonconverged_fits",
                  "tate_estimate", "tate_se", "target_only_estimate", "truth", "all_ok")],
        row.names = FALSE, digits = 4)
  cat(sprintf("omega_star=%s  c1_identity_ok=%s (difference %.6f)\n",
              format(omega_star), c1_ok, c1$c1_reference_difference))
  roce_write_atomic_directory(output, function(stage) {
    write.csv(table, file.path(stage, "omega_selection_table.csv"), row.names = FALSE)
    writeLines(c("task=dgp_common_basis", "history_entry=0002",
                 paste0("omega_star=", omega_star),
                 paste0("c1_identity_ok=", c1_ok),
                 paste0("cells=", nrow(table))),
               file.path(stage, "metadata.txt"))
    files <- list.files(stage, full.names = TRUE)
    writeLines(paste(vapply(files, roce_sha256_file, ""), basename(files), sep = "  "),
               file.path(stage, "sha256.txt"))
  }, caller = "common basis DGP omega selection")
}

if (sys.nframe() == 0L) main()
