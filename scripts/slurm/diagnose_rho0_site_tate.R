#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_directory <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "results/direct_tate_mc500_b5000/rho0_site_tate"
}
n_reference <- if (length(args) >= 2L) as.integer(args[[2L]]) else 500000L
if (is.na(n_reference) || n_reference < 10000L) {
  stop("n_reference must be one integer of at least 10000.", call. = FALSE)
}
if (dir.exists(output_directory)) {
  stop("refusing to overwrite existing output: ", output_directory,
       call. = FALSE)
}

project_library <- Sys.getenv("ROCE_PROJECT_LIB", "")
if (nzchar(project_library)) {
  .libPaths(c(project_library, .libPaths()))
}
suppressPackageStartupMessages(library(RoCE))
source(file.path("scripts", "slurm", "result_provenance.R"))
source(file.path("scripts", "slurm", "atomic_output.R"))
package_provenance <- roce_runtime_package_provenance(project_library)

p <- 100L
kappa <- RoCE:::FACE_KAPPA
source_skewness_max <- RoCE:::FACE_NU_SOURCE_MAX
treatment_shift <- RoCE:::FACE_BINARY_ATE_TARGET
outcome_parameters <- get_face_outcome_parameters(p)
binary_calibration <- RoCE:::get_face_binary_calibration(p, config = "C1", kappa = kappa)
target_tate <- binary_calibration$mu1_superpop -
  binary_calibration$mu0_superpop

# Only the first four FACE outcome coefficients are nonzero. Generating these
# four coordinates directly gives the same marginal outcome law as the p=100
# DGP without allocating an unnecessary n_reference-by-100 matrix.
integrate_site_tate <- function(skewness, seed) {
  effect <- numeric(n_reference)
  set.seed(seed)
  for (j in seq_len(p)) {
    if (outcome_parameters$beta_linear[[j]] == 0 &&
        outcome_parameters$beta_squared[[j]] == 0) {
      next
    }
    x_j <- RoCE:::generate_skewed_normal(
      n_reference, kappa = kappa, phi = 1, nu = skewness
    )
    effect <- effect +
      (x_j - kappa) * outcome_parameters$beta_linear[[j]] +
      x_j^2 * outcome_parameters$beta_squared[[j]]
  }
  logit <- RoCE:::face_binary_logit(effect, binary_calibration)
  risk0 <- plogis(logit)
  risk1 <- plogis(logit + treatment_shift)
  individual_effect <- risk1 - risk0
  c(
    marginal_mu1 = mean(risk1),
    mu1_integration_mcse = stats::sd(risk1) / sqrt(n_reference),
    marginal_mu0 = mean(risk0),
    mu0_integration_mcse = stats::sd(risk0) / sqrt(n_reference),
    marginal_tate = mean(individual_effect),
    tate_integration_mcse = stats::sd(individual_effect) / sqrt(n_reference)
  )
}

parts <- lapply(c(2L, 4L, 8L), function(source_count) {
  skewness <- c(
    target = 0,
    stats::setNames(
      seq(
        source_skewness_max / source_count,
        source_skewness_max,
        length.out = source_count
      ),
      paste0("s", seq_len(source_count))
    )
  )
  rows <- lapply(seq_along(skewness), function(index) {
    integrated <- integrate_site_tate(
      skewness[[index]],
      seed = 420000L + 100L * source_count + index
    )
    data.frame(
      p = p,
      K = source_count,
      site = names(skewness)[[index]],
      skewness = unname(skewness[[index]]),
      conditional_log_odds_shift = treatment_shift,
      marginal_mu1 = unname(integrated[["marginal_mu1"]]),
      mu1_integration_mcse = unname(integrated[["mu1_integration_mcse"]]),
      target_estimand_mu1 = binary_calibration$mu1_superpop,
      difference_from_target_mu1 =
        unname(integrated[["marginal_mu1"]]) -
        binary_calibration$mu1_superpop,
      marginal_mu0 = unname(integrated[["marginal_mu0"]]),
      mu0_integration_mcse = unname(integrated[["mu0_integration_mcse"]]),
      target_estimand_mu0 = binary_calibration$mu0_superpop,
      difference_from_target_mu0 =
        unname(integrated[["marginal_mu0"]]) -
        binary_calibration$mu0_superpop,
      marginal_tate = unname(integrated[["marginal_tate"]]),
      tate_integration_mcse =
        unname(integrated[["tate_integration_mcse"]]),
      target_estimand_tate = target_tate,
      difference_from_target_tate =
        unname(integrated[["marginal_tate"]]) - target_tate,
      n_reference = n_reference,
      package_library = package_provenance$library,
      package_fingerprint = package_provenance$fingerprint,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
})
results <- do.call(rbind, parts)
rownames(results) <- NULL

summary_lines <- c(
  sprintf("p=100 rho=0 site-marginal TATE diagnostic; n=%d per row", n_reference),
  sprintf("Target estimand TATE: %.8f", target_tate),
  sprintf(
    "Largest absolute source-minus-target difference: %.8f",
    max(abs(results$difference_from_target_tate[results$site != "target"]))
  ),
  sprintf(
    "Largest absolute source-minus-target mu1 difference: %.8f",
    max(abs(results$difference_from_target_mu1[results$site != "target"]))
  ),
  sprintf(
    "Largest absolute source-minus-target mu0 difference: %.8f",
    max(abs(results$difference_from_target_mu0[results$site != "target"]))
  ),
  sprintf("Package fingerprint: %s", package_provenance$fingerprint)
)
roce_write_atomic_directory(
  output_directory,
  writer = function(staging_directory) {
    write.csv(
      results,
      file.path(staging_directory, "rho0_site_marginal_tate.csv"),
      row.names = FALSE
    )
    writeLines(
      summary_lines,
      file.path(staging_directory, "rho0_site_marginal_tate_report.txt")
    )
  },
  caller = "rho=0 site-marginal TATE diagnostic"
)
message(paste(summary_lines, collapse = "\n"))
