#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
source(file.path("scripts", "slurm", "grouped_cutoff_diagnostic_helpers.R"))
output_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/grouped_cutoff_pilot/manifest.csv"
}
config <- if (length(args) >= 2L) args[[2L]] else "C1"
source_count <- if (length(args) >= 3L) as.integer(args[[3L]]) else 2L
n_sims <- if (length(args) >= 4L) as.integer(args[[4L]]) else 10L
rho_values <- if (length(args) >= 5L) args[[5L]] else "0;1;2.5"
cutoffs <- if (length(args) >= 6L) args[[6L]] else "1;1.5;2;2.5;3"
simulation_start <- if (length(args) >= 7L) as.integer(args[[7L]]) else 1L

if (!config %in% c("C1", "C2", "C3")) {
  stop("config must be C1, C2, or C3.", call. = FALSE)
}
if (is.na(source_count) || !source_count %in% c(2L, 4L, 8L)) {
  stop("K must be 2, 4, or 8.", call. = FALSE)
}
if (is.na(n_sims) || n_sims < 1L || n_sims > 25L) {
  stop("n_sims must be an integer in 1:25 for this diagnostic.", call. = FALSE)
}
if (is.na(simulation_start) || simulation_start < 1L ||
    simulation_start > 1000000L - n_sims + 1L) {
  stop(
    "simulation_start must define a positive contiguous seed range below 1e6.",
    call. = FALSE
  )
}
rho_grid <- roce_parse_semicolon_numeric_grid(
  rho_values, "rho_values", require_zero_first = TRUE
)
cutoff_grid <- roce_parse_semicolon_numeric_grid(
  cutoffs, "cutoffs", require_positive = TRUE
)
primary_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "ROCE_PRIMARY_CUTOFF", "1"
)))
if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
    primary_cutoff <= 0 ||
    !any(abs(cutoff_grid - primary_cutoff) <= 1e-12)) {
  stop("cutoffs must include the positive ROCE_PRIMARY_CUTOFF.",
       call. = FALSE)
}

manifest <- data.frame(
  task_id = seq_len(n_sims),
  experiment = "grouped_cutoff_diagnostic",
  sim_id = seq.int(simulation_start, length.out = n_sims),
  config = config,
  p = 100L,
  K = source_count,
  rho_values = paste(rho_grid, collapse = ";"),
  cutoffs = paste(cutoff_grid, collapse = ";"),
  primary_cutoff = primary_cutoff,
  n_site = 1000L,
  n_folds = 5L,
  nlambda_init = 100L,
  n_bootstrap = 5000L,
  M_tau = 5,
  M_tau_inference = 5,
  stringsAsFactors = FALSE
)

if (file.exists(output_path) || dir.exists(output_path)) {
  stop("refusing to overwrite existing grouped cutoff manifest: ",
       output_path, call. = FALSE)
}
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile(
  pattern = ".grouped_cutoff_manifest_", tmpdir = dirname(output_path),
  fileext = ".csv"
)
on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
write.csv(manifest, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  stop("failed to atomically commit grouped cutoff manifest.", call. = FALSE)
}
message(
  "wrote ", output_path, " with ", nrow(manifest),
  " seed-level jobs for seeds ", min(manifest$sim_id), "--",
  max(manifest$sim_id), "; each job evaluates rho={",
  paste(rho_grid, collapse = ","), "} and cutoff={",
  paste(cutoff_grid, collapse = ","), "}"
)
}

main()
