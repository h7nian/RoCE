#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
manifest_root <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000"
}
expected_replications <- if (length(args) >= 2L) {
  as.integer(args[[2L]])
} else {
  500L
}
expected_bootstrap <- if (length(args) >= 3L) {
  as.integer(args[[3L]])
} else {
  5000L
}
# Manifest family under audit (HISTORY #0009): the negative-transfer family
# (treated-arm deviation, C1-C3, with its diagnostic manifests) or the
# shared-shift family (both-arm deviation, C1 only, no diagnostic manifests).
family_name <- if (length(args) >= 4L) args[[4L]] else "main"
family <- switch(
  family_name,
  main = list(
    experiment = "negative_transfer", configs = c("C1", "C2", "C3"),
    deviation_mechanism = "treated_arm", audit_diagnostics = TRUE
  ),
  shared_shift = list(
    experiment = "shared_shift", configs = "C1",
    deviation_mechanism = "both_arms", audit_diagnostics = FALSE
  ),
  stop("family must be main or shared_shift.", call. = FALSE)
)
n_configs <- length(family$configs)
expected_primary_cutoff <- suppressWarnings(as.numeric(Sys.getenv(
  "ROCE_PRIMARY_CUTOFF", "1"
)))
expected_n_folds <- suppressWarnings(as.integer(Sys.getenv("ROCE_N_FOLDS", "5")))
cutoff_selection_fingerprint <- tolower(trimws(Sys.getenv(
  "ROCE_CUTOFF_SELECTION_FINGERPRINT", ""
)))

if (is.na(expected_replications) || expected_replications < 1L) {
  stop("expected_replications must be a positive integer.", call. = FALSE)
}
if (is.na(expected_n_folds) || expected_n_folds < 2L) {
  stop("ROCE_N_FOLDS must be one integer of at least 2.", call. = FALSE)
}
if (is.na(expected_bootstrap) || expected_bootstrap < 2L) {
  stop("expected_bootstrap must be an integer >= 2.", call. = FALSE)
}
if (length(expected_primary_cutoff) != 1L ||
    !is.finite(expected_primary_cutoff) || expected_primary_cutoff <= 0) {
  stop("ROCE_PRIMARY_CUTOFF must be one positive finite number.",
       call. = FALSE)
}
if (nzchar(cutoff_selection_fingerprint) &&
    !grepl("^[0-9a-f]{64}$", cutoff_selection_fingerprint)) {
  stop("ROCE_CUTOFF_SELECTION_FINGERPRINT must be SHA-256 when set.",
       call. = FALSE)
}
audit_output_paths <- file.path(
  manifest_root,
  c("manifest_audit_summary.csv", "manifest_audit_passed.txt")
)
if (any(file.exists(audit_output_paths))) {
  stop(
    "refusing to overwrite existing manifest audit artifact(s): ",
    paste(audit_output_paths[file.exists(audit_output_paths)], collapse = ", "),
    call. = FALSE
  )
}

read_manifest <- function(filename) {
  path <- file.path(manifest_root, filename)
  if (!file.exists(path)) {
    stop("manifest not found: ", path, call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

assert_replication_grid <- function(manifest, setting_columns, label) {
  setting_key <- interaction(
    manifest[setting_columns], drop = TRUE, lex.order = TRUE
  )
  by_setting <- split(manifest$sim_id, setting_key)
  expected_ids <- seq_len(expected_replications)
  valid <- vapply(by_setting, function(sim_ids) {
    identical(sort(as.integer(sim_ids)), expected_ids)
  }, logical(1L))
  assert_true(
    all(valid),
    sprintf("%s does not contain exactly sim_id 1:%d in every setting.",
            label, expected_replications)
  )
  length(by_setting)
}

main <- read_manifest("manifest_main.csv")
required_main <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
  "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
  "M_tau", "M_tau_inference", "methods", "deviation_mechanism"
)
assert_true(
  all(required_main %in% names(main)),
  "main manifest is missing required columns."
)
assert_true(
  nrow(main) == n_configs * 3L * 6L * expected_replications,
  "main manifest has the wrong row count."
)
assert_true(
  all(main$experiment == family$experiment) &&
    all(main$deviation_mechanism == family$deviation_mechanism),
  "main manifest experiment or deviation mechanism does not match the family."
)
assert_true(
  identical(as.integer(main$task_id), seq_len(nrow(main))) &&
    !anyDuplicated(main$task_id),
  "main manifest task_id must be the unique sequence 1:n."
)
assert_true(
  all(main$p == 100L) &&
    setequal(main$config, family$configs) &&
    setequal(as.integer(main$K), c(2L, 4L, 8L)) &&
    setequal(main$rho, c(0, 0.5, 1, 1.5, 2, 2.5)) &&
    all(abs(main$cutoff - expected_primary_cutoff) <= 1e-12) &&
    all(main$n_site == 1000L) &&
    all(main$n_folds == expected_n_folds) && all(main$nlambda_init == 100L) &&
    all(main$n_bootstrap == expected_bootstrap) &&
    all(main$M_tau == 5) && all(main$M_tau_inference == 5),
  "main manifest contains an unexpected simulation parameter."
)
main_setting_count <- assert_replication_grid(
  main, c("config", "p", "K", "rho", "cutoff"), "main manifest"
)
assert_true(
  main_setting_count == n_configs * 18L,
  sprintf("main manifest must contain %d settings.", n_configs * 18L)
)

split_parts <- lapply(c(2L, 4L, 8L), function(source_count) {
  split_manifest <- read_manifest(sprintf("manifest_main_K%d.csv", source_count))
  expected_part <- main[main$K == source_count, , drop = FALSE]
  rownames(split_manifest) <- NULL
  rownames(expected_part) <- NULL
  assert_true(
    identical(split_manifest, expected_part),
    sprintf("K=%d split manifest is not an exact ordered subset.", source_count)
  )
  split_manifest
})
split_task_ids <- unlist(lapply(split_parts, `[[`, "task_id"), use.names = FALSE)
assert_true(
  identical(sort(as.integer(split_task_ids)), seq_len(nrow(main))),
  "split manifests do not form an exact partition of main task IDs."
)

rho_groups <- read_manifest("manifest_main_rho_groups.csv")
required_group <- c(
  "group_task_id", "experiment", "sim_id", "config", "p", "K",
  "cutoff", "n_site", "n_folds", "nlambda_init", "n_bootstrap",
  "M_tau", "M_tau_inference", "methods", "deviation_mechanism",
  "rho_values", "primary_task_ids"
)
assert_true(
  all(required_group %in% names(rho_groups)) &&
    nrow(rho_groups) == n_configs * 3L * expected_replications &&
    all(rho_groups$experiment == family$experiment) &&
    all(rho_groups$deviation_mechanism == family$deviation_mechanism) &&
    identical(
      as.integer(rho_groups$group_task_id), seq_len(nrow(rho_groups))
    ) && all(rho_groups$p == 100L) &&
    all(rho_groups$nlambda_init == 100L) &&
    all(abs(rho_groups$cutoff - expected_primary_cutoff) <= 1e-12) &&
    all(rho_groups$n_bootstrap == expected_bootstrap) &&
    all(vapply(rho_groups$rho_values, function(value) {
      identical(
        as.numeric(strsplit(value, ";", fixed = TRUE)[[1L]]),
        c(0, 0.5, 1, 1.5, 2, 2.5)
      )
    }, logical(1L))),
  "rho-group manifest has an invalid schema, size, index, or parameter."
)
expanded_group_ids <- unlist(
  strsplit(rho_groups$primary_task_ids, ";", fixed = TRUE),
  use.names = FALSE
)
expanded_group_ids <- suppressWarnings(as.integer(expanded_group_ids))
assert_true(
  !anyNA(expanded_group_ids) && !anyDuplicated(expanded_group_ids) &&
    identical(sort(expanded_group_ids), seq_len(nrow(main))),
  "rho-group primary task IDs do not exactly partition the main manifest."
)
rho_grid <- c(0, 0.5, 1, 1.5, 2, 2.5)
for (group_index in seq_len(nrow(rho_groups))) {
  task_ids <- as.integer(strsplit(
    rho_groups$primary_task_ids[[group_index]], ";", fixed = TRUE
  )[[1L]])
  mapped <- main[match(task_ids, main$task_id), , drop = FALSE]
  group <- rho_groups[group_index, , drop = FALSE]
  assert_true(
    !anyNA(mapped$task_id) && identical(as.numeric(mapped$rho), rho_grid) &&
      length(unique(mapped$config)) == 1L && mapped$config[[1L]] == group$config &&
      length(unique(mapped$K)) == 1L && mapped$K[[1L]] == group$K &&
      length(unique(mapped$sim_id)) == 1L &&
      mapped$sim_id[[1L]] == group$sim_id,
    sprintf("rho-group row %d does not map exactly to main rows.", group_index)
  )
}

if (family$audit_diagnostics) {
smoke <- read_manifest("manifest_smoke_single.csv")
assert_true(
  nrow(smoke) == 1L && smoke$p == 100L && smoke$config == "C3" &&
    smoke$K == 4L && smoke$rho == 0 && smoke$nlambda_init == 100L &&
    abs(smoke$cutoff - expected_primary_cutoff) <= 1e-12 &&
    smoke$n_bootstrap == expected_bootstrap,
  "smoke manifest does not match the pre-production gate."
)

cutoff <- read_manifest("manifest_cutoff_diagnostic.csv")
required_cutoff <- c(
  "task_id", "experiment", "sim_id", "config", "p", "K", "rho",
  "cutoffs", "primary_cutoff", "n_site", "n_folds", "nlambda_init",
  "M_tau", "M_tau_inference"
)
assert_true(
  all(required_cutoff %in% names(cutoff)) &&
    nrow(cutoff) == 2L * expected_replications && all(cutoff$p == 100L) &&
    all(cutoff$config == "C3") && all(cutoff$K == 4L) &&
    setequal(cutoff$rho, c(0, 2.5)) &&
    all(cutoff$nlambda_init == 100L) &&
    all(abs(cutoff$primary_cutoff - expected_primary_cutoff) <= 1e-12) &&
    all(cutoff$cutoffs == "1,1.5,2,2.5,3"),
  "cutoff manifest contains an unexpected setting."
)
assert_true(
  assert_replication_grid(
    cutoff, c("config", "p", "K", "rho", "cutoffs"), "cutoff manifest"
  ) == 2L,
  "cutoff manifest must contain two settings."
)
assert_true(
  identical(
    as.numeric(cutoff$rho),
    rep(c(0, 2.5), each = expected_replications)
  ) && identical(
    as.integer(cutoff$sim_id),
    rep(seq_len(expected_replications), times = 2L)
  ),
  paste0(
    "cutoff manifest must keep each rho endpoint contiguous with simulation ",
    "IDs 1 through the expected replication count."
  )
)

truncation <- read_manifest("manifest_truncation_diagnostic.csv")
assert_true(
  nrow(truncation) == 2L * 2L * expected_replications &&
    all(truncation$p == 100L) && all(truncation$config == "C3") &&
    all(truncation$K == 4L) && setequal(truncation$rho, c(0, 2.5)) &&
    all(abs(truncation$cutoff - expected_primary_cutoff) <= 1e-12) &&
    all(truncation$nlambda_init == 100L) &&
    all(truncation$n_bootstrap == expected_bootstrap),
  "truncation manifest contains an unexpected setting."
)
assert_true(
  assert_replication_grid(
    truncation,
    c("config", "p", "K", "rho", "M_tau", "M_tau_inference"),
    "truncation manifest"
  ) == 4L,
  "truncation manifest must contain four refit-required settings."
)
assert_true(
  setequal(truncation$M_tau, c(4, 6)) &&
    all(truncation$M_tau_inference == 5) &&
    all(truncation$methods == "one_round_crossfit,target_only"),
  paste0(
    "truncation manifest must contain only the M_fit=4/6 nuisance-refit ",
    "settings; M_fit=5 sensitivities come from primary-task sidecars."
  )
)
}

summary <- data.frame(
  manifest = c("main", "main_K2", "main_K4", "main_K8", "rho_groups"),
  tasks = c(
    nrow(main), vapply(split_parts, nrow, integer(1L)), nrow(rho_groups)
  ),
  settings = n_configs * c(18L, 6L, 6L, 6L, 3L),
  replications_per_setting = rep(expected_replications, 5L),
  bootstrap_draws = rep(expected_bootstrap, 5L),
  primary_cutoff = rep(expected_primary_cutoff, 5L),
  stringsAsFactors = FALSE
)
if (family$audit_diagnostics) {
  summary <- rbind(summary, data.frame(
    manifest = c("smoke", "cutoff", "truncation"),
    tasks = c(nrow(smoke), nrow(cutoff), nrow(truncation)),
    settings = c(1L, 2L, 4L),
    replications_per_setting = c(
      1L, expected_replications, expected_replications
    ),
    bootstrap_draws = c(expected_bootstrap, NA_integer_, expected_bootstrap),
    primary_cutoff = rep(expected_primary_cutoff, 3L),
    stringsAsFactors = FALSE
  ))
}
summary$family <- family_name
summary$deviation_mechanism <- family$deviation_mechanism
utils::write.csv(
  summary, file.path(manifest_root, "manifest_audit_summary.csv"),
  row.names = FALSE
)

manifest_gate <- c(
  "manifest_audit=passed",
  paste0("family=", family_name),
  paste0("experiment=", family$experiment),
  paste0("deviation_mechanism=", family$deviation_mechanism),
  paste0("replications_per_setting=", expected_replications),
  paste0("bootstrap_draws=", expected_bootstrap),
  "p=100",
  paste0("primary_cutoff=", expected_primary_cutoff),
  paste0(
    "cutoff_selection_fingerprint=",
    if (nzchar(cutoff_selection_fingerprint)) {
      cutoff_selection_fingerprint
    } else {
      "not_recorded"
    }
  ),
  paste0(
    "audit_driver_md5=",
    unname(tools::md5sum(file.path(
      "scripts", "slurm", "audit_mc500_manifests.R"
    )))
  ),
  paste0(
    "rho_group_builder_md5=",
    unname(tools::md5sum(file.path(
      "scripts", "slurm", "build_direct_tate_rho_group_manifest.R"
    )))
  ),
  paste0(
    "manifest_main_md5=",
    unname(tools::md5sum(file.path(manifest_root, "manifest_main.csv")))
  ),
  paste0(
    "manifest_main_rho_groups_md5=",
    unname(tools::md5sum(file.path(
      manifest_root, "manifest_main_rho_groups.csv"
    )))
  )
)
for (source_count in c(2L, 4L, 8L)) {
  split_path <- file.path(
    manifest_root, sprintf("manifest_main_K%d.csv", source_count)
  )
  manifest_gate <- c(
    manifest_gate,
    sprintf(
      "manifest_main_K%d_md5=%s", source_count,
      unname(tools::md5sum(split_path))
    )
  )
}
manifest_gate_path <- file.path(manifest_root, "manifest_audit_passed.txt")
manifest_gate_temporary <- tempfile(
  pattern = ".manifest_audit_passed_", tmpdir = manifest_root,
  fileext = ".tmp"
)
writeLines(manifest_gate, manifest_gate_temporary)
if (!file.rename(manifest_gate_temporary, manifest_gate_path)) {
  unlink(manifest_gate_temporary)
  stop("failed to atomically write the manifest audit gate.", call. = FALSE)
}
print(summary, row.names = FALSE)
message("all MC500/B5000 manifest audits passed for family ", family_name)
