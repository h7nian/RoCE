#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/manifest_cutoff_diagnostic.csv"
}
if (file.exists(output_path)) {
  stop("refusing to overwrite existing cutoff manifest: ", output_path,
       call. = FALSE)
}
n_sims <- if (length(args) >= 2L) as.integer(args[[2L]]) else 500L
if (is.na(n_sims) || n_sims < 1L) {
  stop("n_sims must be a positive integer.")
}
primary_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "ROCE_PRIMARY_CUTOFF", "1"
)))
cutoff_grid <- c(1, 1.5, 2, 2.5, 3)
if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
    primary_cutoff <= 0 ||
    !any(abs(cutoff_grid - primary_cutoff) <= 1e-12)) {
  stop("ROCE_PRIMARY_CUTOFF must be one positive member of the cutoff grid.")
}

# One task owns one generated dataset. All cutoffs reuse the same arm-specific
# nuisance fits and fold intermediates inside that task.
manifest <- expand.grid(
  experiment = "c3_cutoff_diagnostic",
  config = "C3",
  p = 100L,
  K = 4L,
  rho = c(0, 2.5),
  sim_id = seq_len(n_sims),
  stringsAsFactors = FALSE
)
manifest$n_site <- 1000L
manifest$n_folds <- 5L
manifest$nlambda_init <- 100L
manifest$M_tau <- 5
manifest$M_tau_inference <- 5
manifest$cutoffs <- paste(cutoff_grid, collapse = ",")
manifest$primary_cutoff <- primary_cutoff
# Keep all replications for one rho endpoint contiguous so bounded MSI batches
# accumulate an interpretable per-setting checkpoint instead of alternating
# between the two diagnostic endpoints.
manifest <- manifest[order(
  manifest$config, manifest$K, manifest$rho, manifest$sim_id
), , drop = FALSE]
rownames(manifest) <- NULL
manifest$task_id <- seq_len(nrow(manifest))
manifest <- manifest[
  ,
  c(
    "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
    "cutoffs", "primary_cutoff", "n_site", "n_folds", "nlambda_init",
    "M_tau", "M_tau_inference"
  )
]

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile(
  pattern = ".direct_tate_cutoff_manifest_", tmpdir = dirname(output_path),
  fileext = ".csv"
)
write.csv(manifest, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path, force = TRUE)
  stop("failed to atomically commit cutoff manifest: ", output_path,
       call. = FALSE)
}
message("wrote ", output_path, " with ", nrow(manifest), " tasks")
