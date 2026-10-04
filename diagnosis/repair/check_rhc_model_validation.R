#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
template_path <- normalizePath(arguments[1L], mustWork = TRUE)
output <- arguments[2L]
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
source("diagnosis/repair/rhc_model_validation.R")
template <- readRDS(template_path)
design <- cbind(1, template$features)
diagnostics <- list()
for (name in names(template$populations)) {
  population <- template$populations[[name]]
  rhc_check_population(population, design)
  target <- population$sites$t
  # Independently enumerate the four treatment/outcome states for each profile,
  # checking both the oracle expectation and its full random-X variance.
  mean_score <- second_moment <- 0
  for (arm in 0:1) for (outcome in 0:1) {
    probability <- target$probability * if (arm == 1L) target$propensity else 1 - target$propensity
    mean_outcome <- target$outcome[, arm + 1L]
    probability <- probability * if (outcome == 1L) mean_outcome else 1 - mean_outcome
    score <- target$outcome[, 2L] - target$outcome[, 1L] +
      if (arm == 1L) (outcome - mean_outcome) / target$propensity else
        -(outcome - mean_outcome) / (1 - target$propensity)
    mean_score <- mean_score + sum(probability * score)
    second_moment <- second_moment + sum(probability * score^2)
  }
  sample <- rhc_generate_population_sample(template, name, 928101L)
  stopifnot(abs(mean_score - population$truth) < 1e-12,
    abs((second_moment - mean_score^2) / sample$data$t$n - sample$oracle["population_se"]^2) < 1e-12)
  for (site_name in setdiff(names(population$sites), "t")) {
    site <- population$sites[[site_name]]
    joint <- site$probability * cbind(1 - site$propensity, site$propensity)
    log_weight <- log(target$probability / joint)
    for (radius in c(2, 3, 5)) diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      scenario = name, site = site_name, radius = radius,
      max_abs_true_log_weight = max(abs(log_weight)),
      true_weight_clipping_probability = sum(joint * (abs(log_weight) > radius)),
      outcome_logit_max = max(abs(qlogis(site$outcome))))
  }
}
# Check the bounded construction on a simple design independently of the RHC fits.
toy_design <- cbind(1, c(-2, 0, 1, 3))
toy_probability <- c(.1, .2, .3, .4)
toy_coefficients <- c(-.7, 2)
bounded <- rhc_bounded_logistic(toy_coefficients, toy_design, toy_probability, 1.4)
stopifnot(max(abs(toy_design %*% bounded$coefficients)) <= 1.4 + 1e-10,
  abs(sum(toy_probability * plogis(toy_design %*% bounded$coefficients)) -
    sum(toy_probability * plogis(toy_design %*% toy_coefficients))) < 1e-10)
tilt <- rhc_joint_tilt(toy_design, toy_probability, matrix(c(-1, 1), nrow = 1L), c(.6, .4), 2.4)
stopifnot(max(abs(colSums(tilt$joint) - c(.6, .4))) < 1e-12,
  max(abs(tilt$joint * exp(tilt$log_weight) - toy_probability)) < 1e-12,
  max(abs(tilt$log_weight)) <= 2.4 + 1e-10)
stopifnot(identical(template$populations$W_tail$truth, template$populations$W_departure$truth),
  identical(template$populations$W_tail$sites$t, template$populations$W_departure$sites$t))
dir.create(output, recursive = TRUE)
write.csv(do.call(rbind, diagnostics), file.path(output, "population_diagnostics.csv"), row.names = FALSE)
paths <- c(template_path, "diagnosis/repair/rhc_model_validation.R", "diagnosis/repair/check_rhc_model_validation.R")
jsonlite::write_json(setNames(lapply(paths, function(path) digest::digest(file = path, algo = "sha256")),
  normalizePath(paths)), file.path(output, "checked_files.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines("Exact oracle expectation/variance, probability bounds and controlled departure passed.",
  file.path(output, "CHECKS_PASSED"))
cat("Independent population checks passed.\n")
