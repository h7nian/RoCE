#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly = TRUE)
stopifnot(length(arguments) == 2L)
output <- arguments[1L]; library_path <- normalizePath(arguments[2L], mustWork = TRUE)
stopifnot(startsWith(output, "/scratch.global/zhan9381/FACE-HD/real_data/rhc/"), !dir.exists(output))
dir.create(output, recursive = TRUE)
.libPaths(c(library_path, .libPaths()))
suppressPackageStartupMessages(library(RoCE, lib.loc = library_path))
source("diagnosis/repair/rhc_validation_helpers.R")
source("diagnosis/repair/rhc_model_validation.R")
observed <- rhc_validation_data()
features <- do.call(rbind, lapply(observed, `[[`, "W_outcome"))
raw_features <- do.call(rbind, lapply(observed, `[[`, "X"))
profile_hash <- vapply(seq_len(nrow(features)), function(i)
  digest::digest(as.numeric(features[i, ]), algo = "sha256"), character(1L))
unique_rows <- which(!duplicated(profile_hash))
profile <- match(profile_hash, profile_hash[unique_rows])
features <- features[unique_rows, , drop = FALSE]
raw_features <- raw_features[unique_rows, , drop = FALSE]
design <- cbind(1, features)
site_ids <- rep(names(observed), vapply(observed, `[[`, integer(1L), "n"))
probabilities <- lapply(names(observed), function(name) tabulate(profile[site_ids == name],
  nrow(features)) / observed[[name]]$n)
names(probabilities) <- names(observed)
pooled_probability <- tabulate(profile, nrow(features)) / length(profile)
target_probability <- .95 * probabilities$t + .05 * pooled_probability
set.seed(9282026L)
outcome_coefficients <- sapply(0:1, function(arm) as.numeric(fit_initial_outcome(
  observed$t$W_outcome, observed$t$Y, observed$t$A, A_val = arm, nlambda = 100L,
  lambda_rule = "min")))
outcome <- plogis(design %*% outcome_coefficients)
propensity_coefficients <- lapply(observed, function(site) as.numeric(fit_initial_outcome(
  site$W_outcome, site$A, rep(1L, site$n), A_val = 1L, nlambda = 100L, lambda_rule = "min")))
case_probability <- lapply(names(observed), function(name) if (name == "t") target_probability else
  .90 * probabilities[[name]] + .10 * pooled_probability)
names(case_probability) <- names(observed)
bounded_propensity <- lapply(names(observed), function(name) rhc_bounded_logistic(
  propensity_coefficients[[name]], design, case_probability[[name]], 2))
names(bounded_propensity) <- names(observed)
case_sites <- lapply(names(observed), function(name) list(probability = case_probability[[name]],
  propensity = plogis(drop(design %*% bounded_propensity[[name]]$coefficients)), outcome = outcome))
names(case_sites) <- names(observed)
truth <- function(sites) sum(sites$t$probability * (sites$t$outcome[, 2L] - sites$t$outcome[, 1L]))
populations <- list(O_case_mix = list(sites = case_sites, truth = truth(case_sites),
  description = "Correct common logistic OR; empirical site case mix with smoothed common support"))

# Preserve each arm's target outcome mean while bounding the linear predictor;
# add a prespecified omitted, bounded age-by-severity interaction.
bounded_outcome <- lapply(0:1, function(arm) rhc_bounded_logistic(
  outcome_coefficients[, arm + 1L], design, target_probability, 1.4))
stopifnot(all(c("age", "aps1") %in% colnames(features)))
interaction <- tanh(features[, "age"]) * tanh(features[, "aps1"])
interaction <- interaction - sum(target_probability * interaction)
nonlinear_outcome <- plogis(cbind(drop(design %*% bounded_outcome[[1L]]$coefficients) - .75 * interaction,
  drop(design %*% bounded_outcome[[2L]]$coefficients) + .75 * interaction))
source_slopes <- lapply(observed[-1L], function(site) sapply(0:1, function(arm) {
  coefficients <- fit_initial_density_ratio(site$Z_site, site$A,
    colMeans(cbind(1, observed$t$Z_site)), A_val = arm, M_tau = 5,
    nlambda = 100L, lambda_rule = "min", tol = 1e-10,
    nuisance_solver = "proximal_newton", nuisance_cv_certificate = TRUE)
  as.numeric(coefficients[-1L])
}))
tilt_scales <- list()
for (scenario in c("W_overlap", "W_tail", "W_departure")) {
  radius <- if (scenario == "W_overlap") 1.4 else 2.4
  sites <- list(t = list(probability = target_probability, propensity = case_sites$t$propensity,
    outcome = nonlinear_outcome))
  for (name in names(source_slopes)) {
    frequency <- mean(observed[[name]]$A)
    tilt <- rhc_joint_tilt(design, target_probability, source_slopes[[name]], c(1 - frequency, frequency), radius)
    probability <- rowSums(tilt$joint)
    source_outcome <- nonlinear_outcome
    if (scenario == "W_departure" && name == "s3")
      source_outcome[, 2L] <- plogis(qlogis(source_outcome[, 2L]) + 1)
    sites[[name]] <- list(probability = probability, propensity = tilt$joint[, 2L] / probability,
      outcome = source_outcome, weight_coefficients = tilt$coefficients)
    tilt_scales[[paste(scenario, name, sep = ":")]] <- tilt$scale
  }
  populations[[scenario]] <- list(sites = sites, truth = truth(sites), weight_radius = radius,
    description = paste("Correct exponential joint weights with log radius", radius,
      "; omitted outcome interaction; source s3 treated log-odds departure", as.integer(scenario == "W_departure")))
}
template <- list(features = features, raw_features = raw_features, observed_structure = observed,
  populations = populations, source_slopes = source_slopes, tilt_scales = tilt_scales,
  outcome_coefficients = outcome_coefficients, propensity_coefficients = propensity_coefficients,
  bounded_propensity = bounded_propensity, bounded_outcome = bounded_outcome,
  seed = 9282026L, library = library_path)
for (population in populations) rhc_check_population(population, design)
saveRDS(template, file.path(output, "template.rds"))
checks <- list()
for (scenario in names(populations)) {
  first <- rhc_generate_population_sample(template, scenario, 928001L)
  second <- rhc_generate_population_sample(template, scenario, 928001L)
  stopifnot(identical(first, second),
    !identical(first$data, rhc_generate_population_sample(template, scenario, 928002L)$data),
    all(vapply(first$data, function(site) identical(site$W_outcome, site$Z_site), logical(1L))))
  checks[[scenario]] <- data.frame(scenario = scenario, truth = first$truth,
    population_oracle_se = unname(first$oracle["population_se"]), unique_profiles = nrow(features))
}
write.csv(do.call(rbind, checks), file.path(output, "population_truth.csv"), row.names = FALSE)
paths <- c("diagnosis/repair/build_rhc_model_template.R", "diagnosis/repair/rhc_model_validation.R",
  "diagnosis/repair/rhc_validation_helpers.R", file.path(output, "template.rds"))
jsonlite::write_json(setNames(lapply(paths, function(path) digest::digest(file = path, algo = "sha256")),
  normalizePath(paths)), file.path(output, "checked_files.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines("Population laws, exact density identities, fixed truth and sample reproducibility checked.",
  file.path(output, "CHECKS_PASSED"))
print(do.call(rbind, checks))
