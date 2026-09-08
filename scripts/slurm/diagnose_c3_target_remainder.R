#!/usr/bin/env Rscript

main <- function(args = commandArgs(trailingOnly = TRUE)) {
if (length(args) > 3L) {
  stop(
    paste0(
      "usage: diagnose_c3_target_remainder.R [N_REPLICATIONS] ",
      "[OUTPUT_ROOT] [NUISANCE_LAMBDA_RULE]"
    ),
    call. = FALSE
  )
}
n_replications <- if (length(args) >= 1L) as.integer(args[[1L]]) else 500L
output_root <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "results/direct_tate_mc500_b5000/c3_target_remainder"
}
nuisance_lambda_rule <- if (length(args) >= 3L) args[[3L]] else "min"
if (!nuisance_lambda_rule %in% c("min", "1se")) {
  stop("NUISANCE_LAMBDA_RULE must be either 'min' or '1se'.", call. = FALSE)
}
if (length(n_replications) != 1L || is.na(n_replications) ||
    n_replications < 1L) {
  stop("N_REPLICATIONS must be one positive integer.", call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "atomic_output.R"))
source(file.path("scripts", "slurm", "direct_tate_task_helpers.R"))
package_provenance <- roce_runtime_package_provenance(project_library)
workflow_fingerprint <- roce_sha256_environment(
  "ROCE_C3_DIAGNOSTIC_WORKFLOW_FINGERPRINT"
)

if (dir.exists(output_root)) {
  stop(
    "diagnostic output directory already exists; use a fresh OUTPUT_ROOT to preserve it.",
    call. = FALSE
  )
}

.true_face_target_means <- function(data) {
  target <- data$R == "t"
  X <- data$X[target, , drop = FALSE]
  centered <- sweep(X, 2L, RoCE:::FACE_KAPPA, "-")
  eta <- as.numeric(
    centered %*% data$out_params$beta_linear +
      X^2 %*% data$out_params$beta_squared
  )
  calibration <- RoCE:::get_face_binary_calibration(
    data$p, config = data$config, kappa = RoCE:::FACE_KAPPA
  )
  baseline <- RoCE:::face_binary_logit(eta, calibration)
  list(
    mu0 = stats::plogis(baseline),
    mu1 = stats::plogis(baseline + RoCE:::FACE_BINARY_ATE_TARGET)
  )
}

.diagnose_one_replication <- function(sim_id) {
  set.seed(sim_id)
  data <- generate_simulation_data(
    n_total = 5000L,
    K = 4L,
    p = 100L,
    config = "C3",
    estimand_type = "superpopulation",
    outcome_type = "binary",
    dgp_type = "face",
    ate_deviation = 0,
    n_deviated_sites = 0L,
    warn_ignored = FALSE
  )
  data_split <- split_data_by_site(data)
  cell_diagnostics <- RoCE:::.binary_cell_diagnostics(
    data_split, "binomial"
  )
  folds <- RoCE:::build_crossfit_folds(
    data_split, n_folds = 5L
  )$target_folds
  target_rows <- which(data$R == "t")
  true_propensity <- data$p_treat_true[target_rows]
  true_means <- .true_face_target_means(data)
  n_target <- data_split$t$n
  propensity_cache <- new.env(hash = TRUE, parent = emptyenv())

  variants <- c(
    fitted_outcome_fitted_propensity = "fitted_fitted",
    true_outcome_fitted_propensity = "oracle_outcome",
    fitted_outcome_true_propensity = "oracle_propensity",
    true_outcome_true_propensity = "oracle_both"
  )
  arm_phi <- lapply(variants, function(unused) {
    list(mu0 = numeric(n_target), mu1 = numeric(n_target))
  })

  for (A_val in c(0L, 1L)) {
    arm_name <- paste0("mu", A_val)
    for (k1 in seq_len(5L)) {
      fit <- RoCE:::estimate_target_only_from_complement(
        target_folds = folds,
        k1 = k1,
        n_folds = 5L,
        family = "binomial",
        A_val = A_val,
        propensity_cache = propensity_cache,
        nuisance_lambda_rule = nuisance_lambda_rule
      )
      evaluation <- RoCE:::materialize_fold(folds, k1)
      index <- evaluation$original_idx
      fitted_propensity <- if (A_val == 1L) {
        fit$prop_scores
      } else {
        1 - fit$prop_scores
      }
      oracle_propensity <- if (A_val == 1L) {
        true_propensity[index]
      } else {
        1 - true_propensity[index]
      }
      oracle_outcome <- true_means[[arm_name]][index]

      predictions <- list(
        fitted_outcome_fitted_propensity = list(
          outcome = fit$m_pred, propensity = fitted_propensity
        ),
        true_outcome_fitted_propensity = list(
          outcome = oracle_outcome, propensity = fitted_propensity
        ),
        fitted_outcome_true_propensity = list(
          outcome = fit$m_pred, propensity = oracle_propensity
        ),
        true_outcome_true_propensity = list(
          outcome = oracle_outcome, propensity = oracle_propensity
        )
      )
      for (variant in names(predictions)) {
        arm_phi[[variant]][[arm_name]][index] <-
          RoCE:::calculate_aipw_pseudo_outcome(
            y = evaluation$Y,
            a = evaluation$A,
            m_hat = predictions[[variant]]$outcome,
            p_a = predictions[[variant]]$propensity,
            A_val = A_val
          )
      }
    }
  }

  truth <- data$mu1_true - data$mu0_true
  rows <- lapply(names(variants), function(variant) {
    phi <- arm_phi[[variant]]$mu1 - arm_phi[[variant]]$mu0
    if (length(phi) != n_target || any(!is.finite(phi))) {
      stop(
        sprintf("replication %d produced invalid %s pseudo-values.", sim_id, variant),
        call. = FALSE
      )
    }
    estimate <- mean(phi)
    se <- sqrt(mean((phi - estimate)^2) / n_target)
    ci_lower <- estimate - stats::qnorm(0.975) * se
    ci_upper <- estimate + stats::qnorm(0.975) * se
    row <- data.frame(
      sim_id = sim_id,
      method = unname(variants[[variant]]),
      estimate = estimate,
      se = se,
      bias = estimate - truth,
      coverage = truth >= ci_lower && truth <= ci_upper,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_width = ci_upper - ci_lower,
      truth = truth,
      experiment = "c3_target_remainder_diagnostic",
      dgp_type = "face",
      outcome_family = "binomial",
      heterogeneity_type = "none",
      estimand_type = "superpopulation",
      config = "C3",
      p = 100L,
      K = 4L,
      rho = 0,
      n_site = 1000L,
      n_folds = 5L,
      nlambda_init = 100L,
      nuisance_lambda_rule = nuisance_lambda_rule,
      estimand_scope = "tate",
      package_library = package_provenance$library,
      package_fingerprint = package_provenance$fingerprint,
      workflow_fingerprint = workflow_fingerprint,
      stringsAsFactors = FALSE
    )
    for (name in names(cell_diagnostics)) {
      row[[name]] <- cell_diagnostics[[name]]
    }
    row
  })
  do.call(rbind, rows)
}

allocated_cores <- suppressWarnings(as.integer(
  Sys.getenv("SLURM_CPUS_PER_TASK", "1")
))
if (length(allocated_cores) != 1L || is.na(allocated_cores) ||
    allocated_cores < 1L) {
  stop("SLURM_CPUS_PER_TASK must be one positive integer.", call. = FALSE)
}
workers <- min(allocated_cores, n_replications)
message(sprintf(
  paste0(
    "[C3 remainder] replications=%d workers=%d p=100 K=4 rho=0 ",
    "folds=5 nuisance_rule=%s"
  ),
  n_replications, workers, nuisance_lambda_rule
))

replication_ids <- seq_len(n_replications)
parts <- if (workers > 1L && .Platform$OS.type != "windows") {
  parallel::mclapply(
    replication_ids,
    .diagnose_one_replication,
    mc.cores = workers,
    mc.preschedule = FALSE
  )
} else {
  lapply(replication_ids, .diagnose_one_replication)
}
failed <- which(vapply(parts, inherits, logical(1L), what = "try-error"))
if (length(failed) > 0L) {
  stop(
    "C3 target-remainder worker failure(s) at replications: ",
    paste(failed, collapse = ","),
    call. = FALSE
  )
}
raw <- do.call(rbind, parts)
rownames(raw) <- NULL
scheduler_provenance <- roce_scheduler_provenance()
for (field in names(scheduler_provenance)) {
  raw[[field]] <- scheduler_provenance[[field]]
}

summary <- RoCE:::diagnose_simulation_results(
  raw,
  expected_replications = n_replications
)
summary <- summary[c(
  "experiment", "dgp_type", "outcome_family", "config",
  "heterogeneity_type", "estimand_type", "p", "K", "rho",
  "estimand_scope", "method", "nuisance_lambda_rule",
  "n_unique_replications",
  "bias", "bias_mcse", "empirical_sd", "bias_skewness",
  "bias_excess_kurtosis",
  "mean_se", "se_to_empirical_sd", "rmse", "coverage", "coverage_mcse",
  "truth_below_interval_fraction", "truth_above_interval_fraction",
  "normal_reference_coverage", "coverage_diagnosis", "diagnostic_status"
)]
summary$package_library <- package_provenance$library
summary$package_fingerprint <- package_provenance$fingerprint
summary$workflow_fingerprint <- workflow_fingerprint
for (field in names(scheduler_provenance)) {
  summary[[field]] <- scheduler_provenance[[field]]
}

paired_estimates <- reshape(
  raw[c("sim_id", "method", "estimate")],
  idvar = "sim_id", timevar = "method", direction = "wide"
)
paired_estimates <- paired_estimates[order(paired_estimates$sim_id), ]
required_paired_columns <- paste0(
  "estimate.",
  c("fitted_fitted", "oracle_outcome", "oracle_propensity", "oracle_both")
)
missing_paired_columns <- setdiff(
  required_paired_columns, names(paired_estimates)
)
if (length(missing_paired_columns) > 0L ||
    nrow(paired_estimates) != n_replications ||
    !identical(as.integer(paired_estimates$sim_id), seq_len(n_replications)) ||
    any(!is.finite(unlist(paired_estimates[required_paired_columns])))) {
  stop(
    "C3 paired remainder decomposition is incomplete: ",
    paste(missing_paired_columns, collapse = ","),
    call. = FALSE
  )
}
interaction_remainder <- with(
  paired_estimates,
  estimate.fitted_fitted - estimate.oracle_outcome -
    estimate.oracle_propensity + estimate.oracle_both
)
paired_remainder <- data.frame(
  sim_id = paired_estimates$sim_id,
  interaction_remainder = interaction_remainder,
  package_fingerprint = package_provenance$fingerprint,
  workflow_fingerprint = workflow_fingerprint,
  stringsAsFactors = FALSE
)
remainder_summary <- data.frame(
  contrast = paste0(
    "fitted_fitted-oracle_outcome-oracle_propensity+oracle_both"
  ),
  n_replications = length(interaction_remainder),
  mean = mean(interaction_remainder),
  empirical_sd = stats::sd(interaction_remainder),
  mcse = stats::sd(interaction_remainder) /
    sqrt(length(interaction_remainder)),
  rmse = sqrt(mean(interaction_remainder^2)),
  package_library = package_provenance$library,
  package_fingerprint = package_provenance$fingerprint,
  workflow_fingerprint = workflow_fingerprint,
  stringsAsFactors = FALSE
)

roce_write_atomic_directory(
  output_root,
  writer = function(staging_directory) {
    write.csv(
      raw,
      file.path(staging_directory, "c3_target_remainder_raw.csv"),
      row.names = FALSE
    )
    write.csv(
      summary,
      file.path(staging_directory, "c3_target_remainder_summary.csv"),
      row.names = FALSE
    )
    write.csv(
      paired_remainder,
      file.path(staging_directory, "c3_target_remainder_paired.csv"),
      row.names = FALSE
    )
    write.csv(
      remainder_summary,
      file.path(staging_directory, "c3_target_remainder_contrast.csv"),
      row.names = FALSE
    )
  },
  caller = "C3 target-remainder diagnostic"
)
print(summary, row.names = FALSE, digits = 5)
print(remainder_summary, row.names = FALSE, digits = 5)
message("[done] wrote versioned diagnostic directory ", output_root)
}

main()
