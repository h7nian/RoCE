#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/manifest.csv"
}
if (file.exists(output_path)) {
  stop("refusing to overwrite existing manifest: ", output_path,
       call. = FALSE)
}
mode <- if (length(args) >= 2L) args[[2L]] else "smoke"
n_sims <- if (length(args) >= 3L) as.integer(args[[3L]]) else 500L
nlambda_init <- if (length(args) >= 4L) as.integer(args[[4L]]) else 100L
if (!mode %in% c("smoke", "main", "shared_shift", "truncation",
                 "nlambda_validation")) {
  stop(paste0(
    "mode must be smoke, main, shared_shift, truncation, or ",
    "nlambda_validation. Use build_direct_tate_cutoff_manifest.R for cutoff ",
    "sensitivity so each dataset reuses one set of nuisance fits."
  ))
}
if (is.na(n_sims) || n_sims < 1L) {
  stop("n_sims must be a positive integer.")
}
if (is.na(nlambda_init) || nlambda_init < 2L) {
  stop("nlambda_init must be one integer of at least 2.")
}
if (!mode %in% c("smoke", "nlambda_validation") && nlambda_init != 100L) {
  stop(
    "primary and truncation manifests are locked to nlambda_init = 100."
  )
}
if (mode == "smoke") {
  n_sims <- 1L
}
if (mode == "nlambda_validation" && n_sims > 10L) {
  stop(
    "nlambda_validation is deliberately bounded to at most 10 paired seeds."
  )
}
primary_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "ROCE_PRIMARY_CUTOFF", "1"
)))
if (length(primary_cutoff) != 1L || !is.finite(primary_cutoff) ||
    primary_cutoff <= 0) {
  stop("ROCE_PRIMARY_CUTOFF must be one positive finite number.")
}

methods <- paste(
  c(
    "one_round_crossfit", "target_only", "sample_size",
    "inverse_variance", "federated_dr", "pooled_dr"
  ),
  collapse = ","
)
nlambda_validation_methods <- paste(
  c("one_round_crossfit", "target_only"),
  collapse = ","
)

default_truncation_settings <- data.frame(
  M_tau = 5,
  M_tau_inference = 5
)

make_grid <- function(experiment, configs, p_values, K_values, rho_values,
                      cutoffs,
                      truncation_settings = default_truncation_settings,
                      method_specification = methods,
                      deviation_mechanism = "treated_arm") {
  required_truncation_columns <- c("M_tau", "M_tau_inference")
  if (!is.data.frame(truncation_settings) ||
      nrow(truncation_settings) < 1L ||
      !all(required_truncation_columns %in% names(truncation_settings))) {
    stop(
      "truncation_settings must contain M_tau and M_tau_inference.",
      call. = FALSE
    )
  }
  grid <- expand.grid(
    experiment = experiment,
    config = configs,
    p = p_values,
    K = K_values,
    rho = rho_values,
    cutoff = cutoffs,
    sim_id = seq_len(n_sims),
    truncation_index = seq_len(nrow(truncation_settings)),
    stringsAsFactors = FALSE
  )
  grid$n_site <- 1000L
  grid$n_folds <- 5L
  grid$nlambda_init <- nlambda_init
  grid$n_bootstrap <- 5000L
  grid$M_tau <- truncation_settings$M_tau[grid$truncation_index]
  grid$M_tau_inference <-
    truncation_settings$M_tau_inference[grid$truncation_index]
  grid$truncation_index <- NULL
  grid$methods <- method_specification
  # How the deviated source deviates (R/data_generation_face.R): the treated
  # log-odds shift of the negative-transfer experiment or the shared shift of
  # both arms (main.tex sec:simulations, shared-shift experiment).
  grid$deviation_mechanism <- deviation_mechanism
  grid
}

parts <- list()
if (mode == "smoke") {
  parts[[length(parts) + 1L]] <- make_grid(
    experiment = "single_task_smoke",
    configs = "C3",
    p_values = 100L,
    K_values = 4L,
    rho_values = 0,
    cutoffs = primary_cutoff
  )
}
if (mode == "main") {
  parts[[length(parts) + 1L]] <- make_grid(
    experiment = "negative_transfer",
    configs = c("C1", "C2", "C3"),
    p_values = 100L,
    K_values = c(2L, 4L, 8L),
    rho_values = c(0, 0.5, 1, 1.5, 2, 2.5),
    cutoffs = primary_cutoff
  )
}
if (mode == "shared_shift") {
  # Pre-registered scope (HISTORY #0009): C1 only, all three source counts,
  # the same rho grid read as a shared log-odds shift of both arms of s1.
  parts[[length(parts) + 1L]] <- make_grid(
    experiment = "shared_shift",
    configs = "C1",
    p_values = 100L,
    K_values = c(2L, 4L, 8L),
    rho_values = c(0, 0.5, 1, 1.5, 2, 2.5),
    cutoffs = primary_cutoff,
    deviation_mechanism = "both_arms"
  )
}
if (mode == "truncation") {
  # Only fitting-radius changes require new nuisance fits. The primary C3/K4
  # rho-endpoint tasks write (5,4), (5,5), (5,6), and (5,Inf) sidecars from
  # their fitted (5,5) objects, so this separate manifest contains just the
  # genuinely new M_fit=4 and M_fit=6 fits.
  truncation_settings <- data.frame(
    M_tau = c(4, 6),
    M_tau_inference = c(5, 5)
  )
  parts[[length(parts) + 1L]] <- make_grid(
    experiment = "c3_truncation_diagnostic",
    configs = "C3",
    p_values = 100L,
    K_values = 4L,
    rho_values = c(0, 2.5),
    cutoffs = primary_cutoff,
    truncation_settings = truncation_settings,
    method_specification = nlambda_validation_methods
  )
}
if (mode == "nlambda_validation") {
  # A short paired-seed runtime/numerical validation. Comparison estimators do
  # not depend on the nuisance-path length and are therefore omitted here;
  # this keeps the diagnostic focused and avoids repeating their bootstrap.
  parts[[length(parts) + 1L]] <- make_grid(
    experiment = "nlambda_validation",
    configs = "C3",
    p_values = 100L,
    K_values = 4L,
    rho_values = 0,
    cutoffs = primary_cutoff,
    method_specification = nlambda_validation_methods
  )
}

manifest <- do.call(rbind, parts)
# Keep all replications from one statistical setting contiguous.  Bounded MSI
# submissions can then accumulate interpretable coverage/RMSE evidence for one
# setting before moving to the next, instead of scattering a tiny batch across
# many settings.  This changes scheduling order only, not any simulation seed
# or data-generating parameter.
manifest <- manifest[order(
  manifest$experiment,
  manifest$K,
  manifest$config,
  manifest$rho,
  manifest$cutoff,
  manifest$M_tau,
  manifest$M_tau_inference,
  manifest$sim_id
), , drop = FALSE]
rownames(manifest) <- NULL
manifest$task_id <- seq_len(nrow(manifest))
manifest <- manifest[
  ,
  c(
    "task_id", "experiment", "sim_id", "config", "p", "K",
    "rho", "cutoff", "n_site", "n_folds", "nlambda_init",
    "n_bootstrap", "M_tau", "M_tau_inference", "methods",
    "deviation_mechanism"
  )
]
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
temporary_path <- tempfile(
  pattern = ".direct_tate_manifest_", tmpdir = dirname(output_path),
  fileext = ".csv"
)
write.csv(manifest, temporary_path, row.names = FALSE)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path, force = TRUE)
  stop("failed to atomically commit manifest: ", output_path,
       call. = FALSE)
}
message("wrote ", output_path, " with ", nrow(manifest), " tasks")
